import sys, os, cPickle, time, argparse
from mpi4py import MPI
from FMM_precond import Mg_CSR
from FMM_Znear import Z_nearCRS_Assembling
from ReadWriteBlitzArray import writeASCIIBlitzArrayToDisk
from scipy import array

def compute_SAIpreconditioner(tmpDirName, C, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers, MAX_BLOCK_SIZE):
    my_id = MPI.COMM_WORLD.Get_rank()
    num_proc = MPI.COMM_WORLD.Get_size()
    # computation of near interactions matrix
    ELEM_TYPE = 'F'
    Z_TMP_ELEM_TYPE = 'F'
    # computation of the Frobenius preconditioner
    Wall_t0 = time.time()
    CPU_t0 = time.clock()
    pathToReadFrom = os.path.join(tmpDirName, 'Z_tmp')
    pathToSaveTo = os.path.join(tmpDirName, 'Mg_LeftFrob')
    # we look for the LIB_G2C type
    file = open('makefile.inc', 'r')
    content = file.readlines()
    file.close()
    for elem in content:
        if 'G2C' in elem:
            LIB_G2C = elem.split('=')[1].split()[0]
    if (my_id == 0):
        print "SAI preconditioner computation..."
    MPI.COMM_WORLD.Barrier()
    Mg_CSR(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()
    CPU_time_Mg_computation = time.clock() - CPU_t0 
    Wall_time_Mg_computation = time.time() - Wall_t0
    # assembling of near interactions matrix
    pathToReadFrom = os.path.join(tmpDirName, 'Z_tmp')
    pathToSaveTo = os.path.join(tmpDirName, 'Z_near')
    Z_nearCRS_Assembling(processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, MAX_BLOCK_SIZE, C, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()
    return Wall_time_Mg_computation, CPU_time_Mg_computation


def compute_SAI(params_simu, simuDirName):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    pathToReadFrom = os.path.join(tmpDirName, 'Z_tmp')
    pathToSaveTo = os.path.join(tmpDirName, 'Mg_LeftFrob')
    if params_simu.COMPUTE_Z_NEAR==1:
        pass
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    #if my_id==0:
        ##print variables['chunkNumber_to_cubesNumbers']
        #print variables['cubeNumber_to_chunkNumber']
    # writing the chunk numbers (per process)
    chunkNumbers = variables['processNumber_to_ChunksNumbers'][my_id]
    writeASCIIBlitzArrayToDisk(array(chunkNumbers).astype('i'), os.path.join(pathToSaveTo, 'chunkNumbers.txt'))
    # writing the cubes numbers (chunk-dependent)
    for chunk in chunkNumbers:
        writeASCIIBlitzArrayToDisk(array(variables['chunkNumber_to_cubesNumbers'][chunk]).astype('i'), os.path.join(tmpDirName, 'Mg_LeftFrob', "chunk" + str(chunk) + 'cubesNumbers.txt'))
    # writing the chunk numbers (cube-dependent)
    writeASCIIBlitzArrayToDisk(array(variables['cubeNumber_to_chunkNumber']).astype('i'), os.path.join(pathToSaveTo, 'cubeNumber_to_chunkNumber.txt')) 
    # SAI precond computation
    Wall_time_Mg_computation, CPU_time_Mg_computation = compute_SAIpreconditioner(tmpDirName, variables['C'], variables['chunkNumber_to_cubesNumbers'], variables['cubeNumber_to_chunkNumber'], variables['chunkNumber_to_processNumber'], variables['processNumber_to_ChunksNumbers'], params_simu.MAX_BLOCK_SIZE)
    variables['Wall_time_Mg_computation'] = Wall_time_Mg_computation
    variables['CPU_time_Mg_computation'] = CPU_time_Mg_computation
    if (my_id == 0) and (params_simu.VERBOSE == 1):
        print variables['CPU_time_Z_near_computation'], "CPU time (seconds) for constructing Z_CFIE_near"
        print variables['Wall_time_Z_near_computation'], "Wall time (seconds) for constructing Z_CFIE_near"
        print variables['CPU_time_Mg_computation'], "CPU time (seconds) for constructing SAI precond"
        print variables['Wall_time_Mg_computation'], "Wall time (seconds) for constructing SAI precond"
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

if __name__=='__main__':
    #MPI.Init()
    my_id = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath('.'))
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--simudir')
    parser.add_argument('--simuparams')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    simuParams = cmdline.simuparams
    if simuDirName==None:
        simuDirName = '.'
    if simuParams==None:
        simuParams = 'simulation_parameters'
    # the simulation itself
    exec 'from ' + simuParams + ' import *'
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        compute_SAI(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()
