import sys, os
from mpi4py import *
from MLFMA import setup_MLFMA, communicateParamsFile

if __name__=='__main__':
    MPI.Init()
    my_id = MPI.COMM_WORLD.Get_rank()
    if (my_id==0):
        if 'result' not in os.listdir('.'):
            os.mkdir('./result')
    # MLFMA parameters
    # we communicate to the others processors the params files
    communicateParamsFile("MLFMA_parameters.py")
    communicateParamsFile("simulation_parameters.py")
    sys.path.append(os.path.abspath('.'))
    # the simulation itself
    from simulation_parameters import *
    if (my_id==0):
        params_simu.display()
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_MLFMA(params_simu)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    MPI.Finalize()

