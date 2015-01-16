import sys, os, time, argparse
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from mpi4py import MPI
from scipy import array, ones, compress
from ReadWriteBlitzArray import writeASCIIBlitzArrayToDisk, writeScalarToDisk
from meshClass import CubeClass

def reduceListRedundancy(listToReduce):
    if len(listToReduce) > 1:
        listToReduce.sort()
        tmp = array(listToReduce, 'i')
        tmp2 = ones(len(tmp), 'i')
        tmp2[1:] = tmp[1:] - tmp[:-1]
        newList = compress(tmp2!=0, listToReduce, axis=0).astype('i')
    else:
        newList = array(listToReduce, 'i')
    return newList

def Mg_listsOfZnearBlocks_ToTransmitAndReceive(ZnearChunkNumber_to_cubesNumbers, ZnearCubeNumber_to_chunkNumber, ZnearChunkNumber_to_processNumber, ZnearProcessNumber_to_ChunksNumbers, pathToReadFrom, Z_TMP_ELEM_TYPE):
    """this function creates 2 lists: Mg_listsOfZ_nearToTransmit and Mg_listsOfZ_nearToReceive"""
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    chunkNumbers = ZnearProcessNumber_to_ChunksNumbers[my_id]
    localPreconditionedCubesNumbers = []
    for i in chunkNumbers:
        localPreconditionedCubesNumbers.append(ZnearChunkNumber_to_cubesNumbers[i])
    listCubesNumbersToReceiveTmp, listCubesNumbersToSendTmp = [], []
    # initialization of the lists
    for i in range(num_proc):
        listCubesNumbersToReceiveTmp.append([])
        listCubesNumbersToSendTmp.append([])
    # we now fill the lists
    for elem in localPreconditionedCubesNumbers:
        for localCube in elem: # elem is a list of cubes Numbers
            chunkNumber = ZnearCubeNumber_to_chunkNumber[localCube]
            pathToReadCubeFrom = os.path.join(pathToReadFrom, "chunk" + str(chunkNumber))
            cube = CubeClass()
            cube.setIntArraysFromFile(pathToReadCubeFrom, localCube)
            for j in cube.cubeNeighborsIndexes:
                ZnearChunkNumber = ZnearCubeNumber_to_chunkNumber[j]
                ZnearProcessNumber = ZnearChunkNumber_to_processNumber[int(ZnearChunkNumber)]
                if not (my_id==ZnearProcessNumber):
                    listCubesNumbersToReceiveTmp[ZnearProcessNumber].append(j)
                    listCubesNumbersToSendTmp[ZnearProcessNumber].append(localCube)
    # we now reduce the redundancy of the lists
    listCubesNumbersToReceive, listCubesNumbersToSend = [], []
    for i in range(num_proc):
        tmp = reduceListRedundancy(listCubesNumbersToReceiveTmp[i])
        listCubesNumbersToReceive.append(tmp)
        tmp = reduceListRedundancy(listCubesNumbersToSendTmp[i])
        listCubesNumbersToSend.append(tmp)
    # now we construct the corresponding chunkNumbers and processNumbers lists
    listChunkNumbersToReceive, listChunkNumbersToSend = [], []
    for L in listCubesNumbersToReceive:
        listChunkNumbersToReceive.append([])
        for i in L:
            listChunkNumbersToReceive[-1].append(ZnearCubeNumber_to_chunkNumber[i])
    for L in listCubesNumbersToSend:
        listChunkNumbersToSend.append([])
        for i in L:
            listChunkNumbersToSend[-1].append(ZnearCubeNumber_to_chunkNumber[i])
    ## we create the missing directories
    for L in listChunkNumbersToReceive:
        for i in L:
            if 'chunk'+ str(i) not in os.listdir(pathToReadFrom):
                os.mkdir(os.path.join(pathToReadFrom, 'chunk'+ str(i)))
    ## now we write the data to be exchanged to disk
    for i in range(num_proc):
        if not (my_id==i):
            writeASCIIBlitzArrayToDisk(array(listCubesNumbersToSend[i]).astype('i'), os.path.join(pathToReadFrom, "CubesNumbersToSendToP" + str(i) + ".txt"))
            writeASCIIBlitzArrayToDisk(array(listChunkNumbersToSend[i]).astype('i'), os.path.join(pathToReadFrom, "ChunkNumbersToSendToP" + str(i) + ".txt"))
    #MPI.COMM_WORLD.Barrier()
    ## finally we write the format of the Near Field matrix elements
    NBytes = 8
    if Z_TMP_ELEM_TYPE=='D':
        NBytes = 16
        print("16 Bytes not supported yet in data transfer in communicateZnearBlocks. Exiting....")
        sys.exit(1)
    writeScalarToDisk(NBytes, os.path.join(pathToReadFrom, "itemsize.txt"))
    MPI.COMM_WORLD.Barrier()

def compute_Z_near(params_simu, simuDirName):
    # computation of near field elements
    my_id = MPI.COMM_WORLD.Get_rank()
    ELEM_TYPE = 'F'
    Z_TMP_ELEM_TYPE = 'F'
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'rb')
    variables = cPickle.load(file)
    file.close()
    Wall_t0 = time.time()
    CPU_t0 = time.clock()
    pathToSaveTo = os.path.join(tmpDirName, 'Z_tmp')
    #Z_nearCRS_Computation(my_id, variables['processNumber_to_ChunksNumbers'], variables['chunkNumber_to_cubesNumbers'], variables['CFIE'], params_simu.MAX_BLOCK_SIZE, variables['w'], params_simu.eps_r, params_simu.mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, params_simu.TDS_APPROX, params_simu.Z_s, params_simu.MOM_FULL_PRECISION, pathToSaveTo)
    # we exchange the missing Z_near parts for each process
    pathToReadFrom = os.path.join(tmpDirName, 'Z_tmp')
    Mg_listsOfZnearBlocks_ToTransmitAndReceive(variables['chunkNumber_to_cubesNumbers'], variables['cubeNumber_to_chunkNumber'], variables['chunkNumber_to_processNumber'], variables['processNumber_to_ChunksNumbers'], pathToReadFrom, 'F')

    # we now dump-pickle the necessary variables
    CPU_time_Z_near_computation = time.clock() - CPU_t0
    Wall_time_Z_near_computation = time.time() - Wall_t0
    variables['Wall_time_Z_near_computation'] = Wall_time_Z_near_computation
    variables['CPU_time_Z_near_computation'] = CPU_time_Z_near_computation
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'wb')
    cPickle.dump(variables, file)
    file.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--inputdir')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    inputDirName = cmdline.inputdir
    simuParams = 'simulation_parameters'

    # the simulation itself
    my_id = MPI.COMM_WORLD.Get_rank()
    if (my_id==0):
        sys.path.append(os.path.abspath(inputDirName))
        exec('from ' + simuParams + ' import *')
    else:
        params_simu = ['blabla']
    params_simu = MPI.COMM_WORLD.bcast(params_simu)
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        compute_Z_near(params_simu, simuDirName)
    else:
        print("you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings.")
        sys.exit(1)
    #MPI.Finalize()

