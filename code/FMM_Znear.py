from mpi4py import MPI
import os, sys, commands
from scipy import zeros, array, floor, compress, reshape
from scipy import sort, argsort, take, arange
from Z_MoM import Z_MoM, Z_MoM_triangles_arraysFromCube
from ReadWriteBlitzArray import *
from meshClass import MeshClass, CubeClass
import copy
from runMPIsystemCommand import runMPIsystemCommand

def Z_near_size_computation(cubes_lists_edges_numbers, cubesNeighborsIndexes):
    N_nearBlockDiag = 0.0
    N_near = 0.0
    C = len(cubes_lists_edges_numbers)
    N_nearPerCube = zeros(C, 'd')
    for i in range(C):
	N_tmp = cubes_lists_edges_numbers[i].shape[0]
        N_nearBlockDiag += 1.0 * N_tmp**2
	N_nearPerCubeTmp = 0
        for j in cubesNeighborsIndexes[i]:
            tmp = cubes_lists_edges_numbers[int(i)].shape[0] * cubes_lists_edges_numbers[int(j)].shape[0]
            N_nearPerCubeTmp += tmp
            N_near += tmp
        N_nearPerCube[i] = N_nearPerCubeTmp
    return N_nearBlockDiag, N_near, N_nearPerCube

def Z_nearPerCube(path, cube, CFIE, cubeNumber, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION):
    """this function computes the part of the MoM Z matrix corresponding to
    the cube 'cubeNumber' as the observer cube and all its neighbors cubes,
    including itself, as the source cubes."""
    # call of the C++ routine for computing the interaction matrix
    Z_CFIE_near_local = Z_MoM_triangles_arraysFromCube(cube, CFIE, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
    # we create a 1-D array out of Z_CFIE_near_local
    Z_CFIE_nearPerCube = array(Z_CFIE_near_local.astype(ELEM_TYPE).flat)
    writeBlitzArrayToDisk(Z_CFIE_near_local.astype(Z_TMP_ELEM_TYPE), os.path.join(path, str(cubeNumber)))
    return Z_CFIE_nearPerCube

def read_Z_perCube_fromFile(path, cubeNumber, cube, Z_TMP_ELEM_TYPE):
    Z_tmp = readBlitzArrayFromDisk(os.path.join(path, str(cubeNumber)), cube.N_RWG_test, cube.N_RWG_src, Z_TMP_ELEM_TYPE)
    return Z_tmp

def compute_list_cubes(cubesNumbers, pathToReadFrom):
    """returns a list of cubes"""
    list_cubes = {}
    for i in range(len(cubesNumbers)):
        cubeNumber = cubesNumbers[i]
        pathToReadCubeFrom = pathToReadFrom
        cube = CubeClass()
        cube.setIntDoubleArraysFromFile(pathToReadCubeFrom, cubeNumber)
        list_cubes[cubeNumber] = copy.copy(cube)
    return list_cubes

def chunk_of_Z_nearCRS_Computation(CFIE, cubesNumbers, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION, pathToSaveTo):
    """this function computes a chunk of the near non-diagonal part of the MoM matrix,
    but saves all the atomic elements on the disk. These elementsz will be later on used 
    by chunk_of_Z_nearCRS_Assembling and MgPrecondition"""
    pathToReadCubeFrom = pathToSaveTo
    list_cubes = compute_list_cubes(cubesNumbers, pathToReadCubeFrom)
    for cubeNumber, cube in list_cubes.iteritems():
        Z_CFIE_near_tmp = Z_nearPerCube(pathToSaveTo, cube, CFIE, cubeNumber, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION)

def chunk_of_Z_nearCRS_Assembling(cubesNumbers, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom):
    """this function computes a chunk of the near non-diagonal part of the MoM matrix
    and uses a scheme similar to CRS for holding the sparse matrix in memory, leading to a
    ~25% decrease in memory consumption of the Z_near matrix and its associated indexes.

    Well, even better than that, since for each cube all RWG functions happen to be test
    functions with exactly the same RWG source functions, we can regroup the column indexes
    per cube and not per RWG function anymore.

    All in all, this leads to a ~50% decrease in memory consumption of the Z_near matrix 
    and its associated indexes with regards to a coordinate scheme (Z_ij, i, j). This new
    scheme is called the Compressed Block Scheme, or CBS."""
    list_cubes = compute_list_cubes(cubesNumbers, pathToReadFrom)
    N_RWG, N_near, N_srcFromNeighbors = 0, 0, 0
    for key, val in list_cubes.iteritems():
        Nl, Nc = val.N_RWG_test, val.N_RWG_src
        N_RWG += Nl
        N_near += Nl * Nc
        N_srcFromNeighbors += Nc
    RWG_numbers = zeros(N_RWG, 'i')
    # number of elements in the Z_near chunk
    Z_CFIE_near = zeros(N_near, ELEM_TYPE) # matrix elements
    # for the q_array, each src function for all the testing functions of a cube appears only once
    # instead of once per testing function. This allows a dramatic reduction in q_array.shape.
    # See the help for this function.
    q_array = zeros(N_srcFromNeighbors, 'i') # column indexes
    rowIndexToColumnIndexes = zeros((N_RWG, 2), 'i') # start and end indexes
    startIndex, startIndexInRWGNumbers, startIndexInQArray = 0, 0, 0
    index_in_rowIndexToColumnIndexes = 0
    for cubeNumber, cube in list_cubes.iteritems():
        # reading the sparse matrix
        Z_CFIE_near_tmp1 = read_Z_perCube_fromFile(pathToReadFrom, cubeNumber, cube, Z_TMP_ELEM_TYPE)
        Z_CFIE_near_tmp = array(Z_CFIE_near_tmp1.astype(ELEM_TYPE).flat)
        Z_CFIE_near[startIndex:startIndex + len(Z_CFIE_near_tmp)] = Z_CFIE_near_tmp.tolist()
        startIndex += Z_CFIE_near_tmp.shape[0] # index update
        # q_array gives the column indexes. This is the second column of pq_array
        q_array[startIndexInQArray:startIndexInQArray + cube.N_RWG_src] = cube.testSrc_RWGsNumbers
        # finding the RWGs numbers
        RWG_numbers[startIndexInRWGNumbers:startIndexInRWGNumbers + cube.N_RWG_test] = cube.testSrc_RWGsNumbers[:cube.N_RWG_test]
        startIndexInRWGNumbers += cube.N_RWG_test
        # now we have to construct rowIndexToColumnIndexes
        for j in range(cube.N_RWG_test):
            rowIndexToColumnIndexes[index_in_rowIndexToColumnIndexes, 0] = startIndexInQArray
            rowIndexToColumnIndexes[index_in_rowIndexToColumnIndexes, 1] = startIndexInQArray + cube.N_RWG_src
            index_in_rowIndexToColumnIndexes += 1
        # update startIndexInQArray
        startIndexInQArray += cube.N_RWG_src
    return Z_CFIE_near, q_array, rowIndexToColumnIndexes, RWG_numbers

def Z_nearChunksDistribution(MAX_BLOCK_SIZE, N_nearPerCube, C, pathToWriteTo, params_simu):
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    if ( (MAX_BLOCK_SIZE<0.1) | (MAX_BLOCK_SIZE>250.) ):
        print "Error: Z_nearChunksDistribution: MAX_BLOCK_SIZE too big or too small"
        sys.exit(1)
    if (my_id==0):
        Total_Z_near_size = sum(N_nearPerCube)*(2.0*4.0/(1024.**2))
        Forecasted_number_of_chunks = max(int(floor(Total_Z_near_size/MAX_BLOCK_SIZE) + 1), num_procs*2)
        Average_N_cubes = C * 1./Forecasted_number_of_chunks
        print "Total size of Z_near matrix =", Total_Z_near_size, "MBytes"
        print "Number of leaf cubes =", C
        print "Forecasted number of chunks =", Forecasted_number_of_chunks
    CUBES_DISTRIBUTION = 1
    writeScalarToDisk(CUBES_DISTRIBUTION, os.path.join('.',pathToWriteTo,'octtree_data/CUBES_DISTRIBUTION.txt') )
    MPI.COMM_WORLD.Barrier()
    # we use the octtree C++ algorithm for repartition of the leaf cubes between the processes
    if (my_id == 0):
        runMPIsystemCommand("./code/MoM", "mpi_mlfma", num_procs)
    MPI.COMM_WORLD.Barrier()

    chunkNumber_to_cubesIndexes, cubeIndex_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = ['blabla'], ['blabla'], ['blabla'], ['blabla']
    if my_id==0:
        cubesIndexAndNumberToProcessNumberTmp = readASCIIBlitzIntArray2DFromDisk(os.path.join('.',pathToWriteTo, 'octtree_data/cubesIndexAndNumberToProcessNumber_FOR_Z_NEAR.txt') )
        # and now new repartition
        # we first sort cubesIndexAndNumberToProcessNumberTmp by the second column (cube number)
        sortedIndexes = argsort(cubesIndexAndNumberToProcessNumberTmp[:,1], axis=-1, kind='mergesort')
        # finally to each cube index we will assign a process number
        cubesIndexToProcessNumberTmp = zeros((C, 2), 'i')
        cubesIndexToProcessNumberTmp[:,0] = take(cubesIndexAndNumberToProcessNumberTmp[:,0], sortedIndexes, axis=0)
        cubesIndexToProcessNumberTmp[:,1] = take(cubesIndexAndNumberToProcessNumberTmp[:,2], sortedIndexes, axis=0)
        # we sort cubesIndexToProcessNumberTmp according to the process number
        sortedIndexes = argsort(cubesIndexToProcessNumberTmp[:,1], axis=-1, kind='mergesort')
        cubesIndexToProcessNumber = take(cubesIndexToProcessNumberTmp, sortedIndexes, axis=0)
        # processNumber_to_cubesIndexes
        processNumber_to_cubesIndexes = {}
        for i in range(num_procs):
            processNumber_to_cubesIndexes[i] = []
        for elem in cubesIndexToProcessNumber:
            cubeIndex, ProcNum = elem[0], elem[1]
            processNumber_to_cubesIndexes[ProcNum].append(cubeIndex)
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

def Z_nearCRS_Computation(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, CFIE, MAX_BLOCK_SIZE, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION, pathToSaveTo):
    """this function computes Z_CFIE_near by slices and stores them on the disk."""
    index, percentage = 0, 0
    for chunkNumber in processNumber_to_ChunksNumbers[my_id]:
        if my_id==0:
            newPercentage = int(index * 100.0/len(processNumber_to_ChunksNumbers[my_id]))
            if (newPercentage - percentage)>=5:
                print "Process", my_id, ": computing Z_CFIE_near chunk.", newPercentage, "% completed"
                sys.stdout.flush()
                percentage = newPercentage
        pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
        cubesNumbers = chunkNumber_to_cubesNumbers[chunkNumber]
        chunk_of_Z_nearCRS_Computation(CFIE, cubesNumbers, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION, pathToSaveToChunk)
        index += 1

def Z_nearCRS_Assembling(processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, MAX_BLOCK_SIZE, C, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom, pathToSaveTo):
    """this function computes Z_CFIE_near by slices and stores them on the disk.
    The maximum size of a block is given by the variable MAX_BLOCK_SIZE in MegaBytes"""
    # test on MAX_BLOCK_SIZE
    if ( (MAX_BLOCK_SIZE<0.1) | (MAX_BLOCK_SIZE>250.) ):
        print "Error: MAX_BLOCK_SIZE too big or too small"
        sys.exit(1)
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    NAME = "Z_CFIE_near"
    if (my_id==0):
        print "Number of leaf cubes =", C
        print "assembling Z_CFIE_near chunks..."
    chunkNumbers = processNumber_to_ChunksNumbers[my_id]
    for chunkNumber in chunkNumbers:
        cubesNumbers = chunkNumber_to_cubesNumbers[chunkNumber]
        pathToReadFromChunk = os.path.join(pathToReadFrom, "chunk" + str(chunkNumber))
        Z_CFIE_near, q_array, rowIndexToColumnIndexes, RWG_numbers = chunk_of_Z_nearCRS_Assembling(cubesNumbers, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFromChunk)
        writeToDisk_chunk_of_Z_sparse(pathToSaveTo, NAME, Z_CFIE_near, q_array, rowIndexToColumnIndexes, RWG_numbers, chunkNumber)
        del Z_CFIE_near, q_array, rowIndexToColumnIndexes, RWG_numbers
        commands.getoutput("rm -rf " + os.path.join(pathToReadFromChunk))
    # we write the chunks numbers of the process
    writeASCIIBlitzArrayToDisk(array(chunkNumbers).astype('i'), os.path.join(pathToSaveTo, 'chunkNumbers.txt'))
    #commands.getoutput("rm -rf " + os.path.join(pathToReadFrom))
