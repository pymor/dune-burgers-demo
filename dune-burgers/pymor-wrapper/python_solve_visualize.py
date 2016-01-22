import sys

from pymor.discretizations.mpi import mpi_wrap_discretization
from pymor.tools import mpi
from pymor.vectorarrays.mpi import MPIVectorArrayAutoComm

from dune_burgers import discretize_dune_burgers

filename = sys.argv[1]
exponent = float(sys.argv[2])
times = map(int, sys.argv[3:])

if mpi.parallel:
    d = mpi_wrap_discretization(lambda: discretize_dune_burgers(filename),
                                use_with=True, with_apply2=False, array_type=MPIVectorArrayAutoComm)
else:
    d = discretize_dune_burgers(filename)

U = d.solve(exponent)
U_vis = U.empty()
U_vis.append(U, o_ind=times, remove_from_other=times)
d.visualize(U_vis, filename='out')
