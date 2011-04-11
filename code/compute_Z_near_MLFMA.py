import sys, os
from mpi4py import *
from MLFMA import compute_Z_near

if __name__=='__main__':
    #MPI.Init()
    my_id = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath('.'))
    # the simulation itself
    from simulation_parameters import *
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        compute_Z_near(params_simu)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    MPI.Finalize()

