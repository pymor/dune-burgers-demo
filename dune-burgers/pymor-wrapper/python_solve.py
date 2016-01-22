import sys

from pymor.discretizations.mpi import mpi_wrap_discretization
from pymor.vectorarrays.mpi import MPIVectorArrayAutoComm

from dune_burgers import discretize_dune_burgers

filename = sys.argv[1]
exponent = float(sys.argv[2])

d = mpi_wrap_discretization(lambda: discretize_dune_burgers(filename),
                            use_with=True, with_apply2=False, array_type=MPIVectorArrayAutoComm)

U = d.solve(exponent)
