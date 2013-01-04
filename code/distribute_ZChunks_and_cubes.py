import sys, os, argparse
import time, copy, commands, pickle, cPickle
from mpi4py import MPI
from scipy import zeros, floor, array
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk, readASCIIBlitzIntArray2DFromDisk, writeBlitzArrayToDisk

def Z_near_size_computation(cubes_lists_edges_numbers, cubes_lists_NeighborsIndexes):
    C = len(cubes_lists_edges_numbers)
    N_nearPerCube = zeros(C, 'd')
    for i in range(C):
	N_nearPerCubeTmp = 0
        for j in cubes_lists_NeighborsIndexes[i]:
            tmp = cubes_lists_edges_numbers[int(i)].shape[0] * cubes_lists_edges_numbers[int(j)].shape[0]
            N_nearPerCubeTmp += tmp
        N_nearPerCube[i] = N_nearPerCubeTmp
    return N_nearPerCube

def Z_nearChunksDistribution(MAX_BLOCK_SIZE, N_nearPerCube, C, pathToWriteTo):
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    if ( (MAX_BLOCK_SIZE<0.1) | (MAX_BLOCK_SIZE>10000.) ):
        print "Error: Z_nearChunksDistribution: MAX_BLOCK_SIZE too big or too small"
        sys.exit(1)

    chunkNumber_to_cubesIndexes, cubeIndex_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = ['blabla'], ['blabla'], ['blabla'], ['blabla']
    if my_id==0:
        Total_Z_near_size = sum(N_nearPerCube)*(2.0*4.0/(1024.**2))
        Forecasted_number_of_chunks = max(int(floor(Total_Z_near_size/MAX_BLOCK_SIZE) + 1), num_procs*2)
        Average_N_cubes = C * 1./Forecasted_number_of_chunks
        print "Total size of Z_near matrix =", Total_Z_near_size, "MBytes"
        print "Number of leaf cubes =", C
        print "Forecasted number of chunks =", Forecasted_number_of_chunks
        cubesIndexAndNumberToProcessNumber = readASCIIBlitzIntArray2DFromDisk(os.path.join(pathToWriteTo, 'octtree_data/cubesIndexAndNumberToProcessNumber_FOR_Z_NEAR.txt') )
        # processNumber_to_cubesIndexes
        processNumber_to_cubesIndexes = {}
        for i in range(num_procs):
            processNumber_to_cubesIndexes[i] = []
        for elem in cubesIndexAndNumberToProcessNumber:
            cubeIndex, ProcNum = elem[0], elem[2]
            processNumber_to_cubesIndexes[ProcNum].append(cubeIndex)
        for i in range(num_procs):
            processNumber_to_cubesIndexes[i].sort() # we sort by cube index
        # processNumber_to_ChunksNumbers, chunkNumber_to_cubesIndexes
        processNumber_to_ChunksNumbers, chunkNumber_to_cubesIndexes, chunkNumber_to_processNumber= [], [], []
        # cubeIndex_to_chunkNumber
        cubeIndex_to_chunkNumber = zeros(C, 'i')
        for i in range(num_procs):
            processNumber_to_ChunksNumbers.append([])
        chunkNumber = 0
        for i in range(num_procs):
            elem = processNumber_to_cubesIndexes[i]
            NC = len(elem) # number of cubes per process
            startIndex, stopIndex = 0, int(min(floor(Average_N_cubes), NC))
            while (startIndex<NC):
                processNumber_to_ChunksNumbers[i].append(chunkNumber)
                chunkNumber_to_cubesIndexes.append(elem[startIndex:stopIndex])
                # cubeIndex_to_chunkNumber
                cubeIndex_to_chunkNumber[elem[startIndex:stopIndex]] = chunkNumber
                chunkNumber_to_processNumber.append(i)
                startIndex = stopIndex
                stopIndex += int(floor(Average_N_cubes))
                stopIndex = min(stopIndex, NC)
                chunkNumber += 1
    return chunkNumber_to_cubesIndexes, cubeIndex_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers

def createChunkDirs(ZprocessNumber_to_ChunksNumbers, tmpDirName, my_id):
    pathToSaveTo = os.path.join(tmpDirName, 'Z_tmp')
    for chunkNumber in ZprocessNumber_to_ChunksNumbers[my_id]:
        pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
        os.mkdir(pathToSaveToChunk)

def distribute_Chunks(params_simu, simuDirName):
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    if my_id==0:
        file = open(os.path.join(tmpDirName, "mesh", 'cubes_lists_RWGsNumbers.txt'), 'r')
        cubes_lists_RWGsNumbers = cPickle.load(file)
        file.close()
        file = open(os.path.join(tmpDirName, "mesh", 'cubes_lists_NeighborsIndexes.txt'), 'r')
        cubes_lists_NeighborsIndexes = cPickle.load(file)
        file.close()
        N_nearPerCube = Z_near_size_computation(cubes_lists_RWGsNumbers, cubes_lists_NeighborsIndexes)
    else:
        N_nearPerCube = ['blabla']
    N_nearPerCube = MPI.COMM_WORLD.bcast(N_nearPerCube)
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = Z_nearChunksDistribution(params_simu.MAX_BLOCK_SIZE, N_nearPerCube, variables['C'], tmpDirName)
    chunkNumber_to_cubesNumbers = MPI.COMM_WORLD.bcast(chunkNumber_to_cubesNumbers)
    cubeNumber_to_chunkNumber = MPI.COMM_WORLD.bcast(cubeNumber_to_chunkNumber)
    chunkNumber_to_processNumber = MPI.COMM_WORLD.bcast(chunkNumber_to_processNumber)
    processNumber_to_ChunksNumbers = MPI.COMM_WORLD.bcast(processNumber_to_ChunksNumbers)
    # distributing chunks of the mesh
    local_ChunksNumbers = processNumber_to_ChunksNumbers[my_id]
    local_chunkNumber_to_cubesNumbers = []
    local_chunkNumber_N_cubesNumbers = []
    startIndex = 0
    for chunkNumber in local_ChunksNumbers:
        list_cubes_tmp = chunkNumber_to_cubesNumbers[chunkNumber]
        local_chunkNumber_to_cubesNumbers += list_cubes_tmp
        length = len(list_cubes_tmp)
        local_chunkNumber_N_cubesNumbers += [length]
    writeBlitzArrayToDisk(array(local_ChunksNumbers, 'i'), os.path.join(tmpDirName, 'Z_tmp', 'local_ChunksNumbers.txt'))
    writeBlitzArrayToDisk(array(local_chunkNumber_to_cubesNumbers, 'i'), os.path.join(tmpDirName, 'Z_tmp', 'local_chunkNumber_to_cubesNumbers.txt'))
    writeBlitzArrayToDisk(array(local_chunkNumber_N_cubesNumbers, 'i'), os.path.join(tmpDirName, 'Z_tmp', 'local_chunkNumber_N_cubesNumbers.txt'))
    writeScalarToDisk(len(local_ChunksNumbers), os.path.join(tmpDirName, 'Z_tmp', "N_local_Chunks.txt"))
    writeScalarToDisk(len(local_chunkNumber_to_cubesNumbers), os.path.join(tmpDirName, 'Z_tmp', "N_local_cubes.txt"))

    createChunkDirs(processNumber_to_ChunksNumbers, tmpDirName, my_id)

    variables['chunkNumber_to_cubesNumbers'] = chunkNumber_to_cubesNumbers
    variables['cubeNumber_to_chunkNumber'] = cubeNumber_to_chunkNumber
    variables['chunkNumber_to_processNumber'] = chunkNumber_to_processNumber
    variables['processNumber_to_ChunksNumbers'] = processNumber_to_ChunksNumbers
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

if __name__=='__main__':
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
    exec 'from ' + simuParams + ' import *'
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        my_id = MPI.COMM_WORLD.Get_rank()
        CPU_time, Wall_time = time.clock(), time.time()
        distribute_Chunks(params_simu, simuDirName)
        CPU_time, Wall_time = time.clock() - CPU_time, time.time() - Wall_time
        print "Process", my_id, "chunks numbers distribution/folders creation: CPU time =", CPU_time, "sec"
        print "Process", my_id, "chunks numbers distribution/folders creation: Wall time =", Wall_time, "sec"
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()

