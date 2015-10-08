from __future__ import absolute_import, division, print_function

import os
import subprocess
import weakref

import numpy as np

from pymor.algorithms.timestepping import ExplicitEulerTimeStepper
from pymor.core.defaults import defaults
from pymor.core.interfaces import ImmutableInterface
from pymor.discretizations.basic import InstationaryDiscretization
from pymor.operators.basic import OperatorBase
from pymor.parameters.spaces import CubicParameterSpace
from pymor.tools import mpi
from pymor.vectorarrays.interfaces import VectorSpace
from pymor.vectorarrays.list import VectorInterface, ListVectorArray
from pymor.vectorarrays.numpy import NumpyVectorSpace

import libdune_burgers as dune_module


class DuneVector(VectorInterface):

    def __init__(self, impl):
        self._impl = impl

    @classmethod
    def make_zeros(cls, subtype):
        assert isinstance(subtype, int)
        impl = dune_module.Vector(subtype, 0.)
        return cls(impl)

    @property
    def subtype(self):
        return self._impl.dim

    @property
    def dim(self):
        return self._impl.dim

    @property
    def data(self):
        return np.frombuffer(self._impl.buffer())

    def copy(self):
        return DuneVector(self._impl.copy())

    def scal(self, alpha):
        self._impl.scal(alpha)

    def axpy(self, alpha, x):
        self._impl.axpy(alpha, x._impl)

    def dot(self, other):
        return self._impl.dot(other._impl)

    def l1_norm(self):
        return self._impl.l1Norm()

    def l2_norm(self):
        return self._impl.l2Norm()

    def sup_norm(self):
        return self._impl.supNorm()

    def components(self, component_indices):
        impl = self._impl
        return np.array([impl[int(i)] for i in component_indices])

    def amax(self):
        ind = self._impl.amaxInd()
        return ind, abs(self._impl[ind])

    def __getstate__(self):
        return (self.dim, self.data)

    def __setstate__(self, state):
        dim, data = state
        self._impl = dune_module.Vector(dim, 0.)
        self.data[:] = data


class DuneSpaceOperator(OperatorBase):

    linear = False

    def __init__(self, d):
        self._impl = impl = d.go
        self._d = weakref.ref(d)
        self.source = self.range = VectorSpace(ListVectorArray, (DuneVector, impl.dimSource))
        self.name = 'DuneBurgersSpaceOperator'
        self.build_parameter_type({'exponent': tuple()}, local_global=True)

    def apply(self, U, ind=None, mu=None):
        assert U in self.source
        mu = self.parse_parameter(mu)
        exponent = float(mu['exponent'])

        vectors = U._list if ind is None else [U._list[i] for i in ind]
        R = self.range.zeros(len(vectors))
        for r, u in zip(R._list, vectors):
            self._impl.apply(u._impl, r._impl, exponent, 1.)

        return R

    def restricted(self, components):
        if isinstance(components, np.ndarray):
            components = components.tolist()

        dims = mpi.comm.allgather(self.source.dim)
        offsets = np.cumsum(np.concatenate(([0], dims)))[:-1]
        offset = offsets[mpi.rank]

        diameter, local_centers, local_source_dofs = self._d().getSubgridData(components, offset)
        local_centers = np.array(local_centers, dtype=np.float64)
        local_source_dofs = np.array(local_source_dofs, dtype=np.int64)

        centers = np.empty((mpi.size,) + local_centers.shape, dtype=np.float64)
        source_dofs = np.empty((mpi.size,) + local_source_dofs.shape, dtype=np.int64) if mpi.rank0 else None
        mpi.comm.Gather(local_centers, centers, root=0)
        mpi.comm.Gather(local_source_dofs, source_dofs, root=0)

        if mpi.rank0:
            centers = np.sum(centers, axis=0)
            source_dofs = np.sum(source_dofs, axis=0)
            op, source_dofs, range_dofs = self._d().makeRestrictedSpaceOperator(diameter, centers.tolist(),
                                                                                source_dofs.tolist())
            source_dofs = np.array(source_dofs)
            range_dofs = np.array(range_dofs)
            return (RestrictedDuneSpaceOperator(op, range_dofs), source_dofs)


class RestrictedDuneSpaceOperator(OperatorBase):

    linear = False

    def __init__(self, impl, range_dofs):
        self._impl = impl
        self._range_dofs = range_dofs
        self.source = NumpyVectorSpace(impl.dimSource)
        self.range = NumpyVectorSpace(len(range_dofs))
        self.name = 'DuneBurgersSpaceOperator_restricted'

        self._source_vec = dune_module.Vector(impl.dimSource, 0.)
        self._range_vec = dune_module.Vector(impl.dimRange, 0.)
        self._source_array = np.frombuffer(self._source_vec.buffer())
        self._range_array = np.frombuffer(self._range_vec.buffer())
        self.build_parameter_type({'exponent': tuple()}, local_global=True)

    def apply(self, U, ind=None, mu=None):
        assert U in self.source
        mu = self.parse_parameter(mu)
        exponent = float(mu['exponent'])

        U = U.data if ind is None else \
            U.data[ind] if hasattr(ind, '__len__') else \
            U.data[ind:ind + 1]
        R = self.range.zeros(len(U))
        R_array = R.data

        for i, u in enumerate(U):
            self._source_array[:] = u
            self._range_array[:] = 0
            self._impl.apply(self._source_vec, self._range_vec, exponent, 1.)
            R_array[i] = self._range_array[:][self._range_dofs]

        return R


class DuneBurgersVisualizer(ImmutableInterface):

    def __init__(self, impl):
        self._impl = impl

    def visualize(self, U, d, filename=None):
        assert isinstance(U, ListVectorArray) \
            or (isinstance(U, tuple) and all(isinstance(u, ListVectorArray) for u in U)
                and all(len(u) == len(U[0]) for u in U))
        U = (U,) if isinstance(U, ListVectorArray) else U
        assert filename is None or len(U) == 1

        sources = []
        for i, u in enumerate(U):
            fn = filename or 'dune-burgers-visualization-{}'.format(i)
            sources.append(fn)

            if len(u) == 1:
                self._impl.visualize(u._list[0]._impl, fn)
            else:
                for i, uu in enumerate(u._list):
                    self._impl.visualize(uu._impl, '{}-{:05}'.format(fn, i))

        if filename is None:
            from pymor.tools import mpi
            if mpi.rank0:
                subprocess.call(['paraview', 's{:0>4}-dune-burgers-visualization-0-..pvtu'.format(mpi.size)])
                import glob
                for fn in glob.glob('s*-dune-burgers-visualization-*'):
                    os.remove(fn)


def discretize_dune_burgers(filename, exponent_range=(1., 2.), cache_region=None):

    impl = dune_module.Discretization(filename)

    operator = DuneSpaceOperator(impl)

    initial_data = operator.source.zeros()
    impl.initialProjection(initial_data._list[0]._impl)

    nt = impl.nt
    T = impl.T
    time_stepper = ExplicitEulerTimeStepper(nt)

    parameter_space = CubicParameterSpace({'exponent': tuple()}, *exponent_range)
    visualizer = DuneBurgersVisualizer(impl)

    d = InstationaryDiscretization(T, initial_data, operator, time_stepper=time_stepper,
                                   parameter_space=parameter_space,
                                   visualizer=visualizer, name='DuneBurgers',
                                   cache_region=cache_region)
    # d.generate_sid()
    return d
