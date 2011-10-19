import sys, os, cPickle, time, argparse
from mpi4py import *
from FMM_Znear import Z_nearCRS_Computation
from FMM_precond import Mg_listsOfZnearBlocks_ToTransmitAndReceive


def compute_Z_near(params_simu, simuDirName):
    # computation of near field elements
    my_id = MPI.COMM_WORLD.Get_rank()
    ELEM_TYPE = 'F'
    Z_TMP_ELEM_TYPE = 'F'
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    Wall_t0 = time.time()
    CPU_t0 = time.clock()
    pathToSaveTo = os.path.join(tmpDirName, 'Z_tmp')
    Z_nearCRS_Computation(my_id, variables['processNumber_to_ChunksNumbers'], variables['chunkNumber_to_cubesNumbers'], variables['CFIE'], params_simu.MAX_BLOCK_SIZE, variables['w'], params_simu.eps_r, params_simu.mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, params_simu.TDS_APPROX, params_simu.Z_s, params_simu.MOM_FULL_PRECISION, pathToSaveTo)
    # we exchange the missing Z_near parts for each process
    pathToReadFrom = os.path.join(tmpDirName, 'Z_tmp')
    Mg_listsOfZnearBlocks_ToTransmitAndReceive(variables['chunkNumber_to_cubesNumbers'], variables['cubeNumber_to_chunkNumber'], variables['chunkNumber_to_processNumber'], variables['processNumber_to_ChunksNumbers'], pathToReadFrom, 'F')
    # we now dump-pickle the necessary variables
    CPU_time_Z_near_computation = time.clock() - CPU_t0
    Wall_time_Z_near_computation = time.time() - Wall_t0
    variables['Wall_time_Z_near_computation'] = Wall_time_Z_near_computation
    variables['CPU_time_Z_near_computation'] = CPU_time_Z_near_computation
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()


if __name__=='__main__':
    #MPI.Init()
    sys.path.append(os.path.abspath('.'))
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    if simuDirName==None:
        simuDirName = '.'

    # the simulation itself
    from simulation_parameters import *
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        compute_Z_near(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    MPI.Finalize()

