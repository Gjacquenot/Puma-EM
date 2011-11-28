from mpi4py import MPI
import sys, os, copy
from scipy import zeros, array, floor, compress, arange, take, product, dot, sort, argsort, ones, weave
from scipy.weave import converters
from scipy import prod, rand, eye, put, transpose #, linalg
from myPseudoInv import computeMyPinvCC #, computeMyPinv
from FMM_Znear import read_Z_perCube_fromFile, compute_list_cubes, writeToDisk_chunk_of_Z_sparse
from ReadWriteBlitzArray import writeBlitzArrayToDisk, readBlitzArrayFromDisk, writeScalarToDisk, writeASCIIBlitzArrayToDisk
from meshClass import CubeClass

def findEdgesInRadiusAroundCube(RWGNumber_nodes, nodesCoord, rCubeCenter, R_NORM_TYPE_1):
    """this function yields a 1 to an edge number if it is within
    a certain 1-Norm around the center of the cube of interest"""
    isEdgeInCartesianRadius = ones(RWGNumber_nodes.shape[0], 'i')
    wrapping_code = """
    int N_edges = isEdgeInCartesianRadius.size();
    for (int i=0 ; i<N_edges ; ++i) {
      const int node1 = RWGNumber_nodes(i, 1), node2 = RWGNumber_nodes(i, 2);
      for (int j=0 ; j<3 ; ++j) {
        const double center_coord = (nodesCoord(node1, j) + nodesCoord(node2, j))/2.0;
        if (std::abs(center_coord - rCubeCenter(j)) > R_NORM_TYPE_1) {
          isEdgeInCartesianRadius(i) = 0;
          break;
        }
      }
    }
    """
    weave.inline(wrapping_code,
             ['isEdgeInCartesianRadius','RWGNumber_nodes','nodesCoord','rCubeCenter','R_NORM_TYPE_1'],
             type_converters = converters.blitz,
             compiler = 'gcc',
             extra_compile_args = ['-O3', '-pthread', '-w'])
    # old code
    #diff = (abs(edges_centroids - rCubeCenter) < R_NORM_TYPE_1).astype('i')
    #isEdgeInCartesianRadius = product(diff, axis=1).astype('i')
    return isEdgeInCartesianRadius

def computePreconditionerColumnsPerCube(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber, R_NORM_TYPE_1, Z_TMP_ELEM_TYPE):
    """this function computes the number of columns per cube for the Frobenius preconditioner"""
    NumberOfColumnsPerCube = zeros(len(list_cubes), 'i')
    list_Z_tmp = {}
    i = 0
    for cubeNumber, cube in list_cubes.iteritems():
        areEdgesInRadiusAroundCube = findEdgesInRadiusAroundCube(cube.localTestSrcRWGNumber_nodes, cube.nodesCoord, cube.rCubeCenter, R_NORM_TYPE_1)
        cube.localEdgesInRadiusAroundCube = (compress(areEdgesInRadiusAroundCube==1, arange(len(areEdgesInRadiusAroundCube)))).astype('i')
        chunkNumber = cubeNumber_to_chunkNumber[cubeNumber]
        pathToReadFromChunk = os.path.join( pathToReadFrom, "chunk" + str(chunkNumber) )
        list_Z_tmp[cubeNumber] = read_Z_perCube_fromFile(pathToReadFromChunk, cubeNumber, cube, Z_TMP_ELEM_TYPE)
        NumberOfColumnsPerCube[i] = sum(areEdgesInRadiusAroundCube)
        i += 1
    return NumberOfColumnsPerCube, list_Z_tmp

def numberOfElemsInPrecond(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber, R_NORM_TYPE_1, Z_TMP_ELEM_TYPE):
    """computes the number of elements in the preconditioner"""
    preconditionerColumnsPerCube, list_Z_tmp = computePreconditionerColumnsPerCube(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber, R_NORM_TYPE_1, Z_TMP_ELEM_TYPE)
    N_precond = 0
    i = 0
    for cubeNumber, cube in list_cubes.iteritems():
        N_precond += preconditionerColumnsPerCube[i] * cube.N_RWG_test
        i += 1
    return N_precond, preconditionerColumnsPerCube, list_Z_tmp

def assignValuesToMatrix(M_toFill, M_toCopy, indexes_lines, indexes_columns):
    wrapping_code = """
    using namespace blitz;
    int Nl = indexes_lines.extent(0), Nc = indexes_columns.extent(0);
    for (int i=0 ; i<Nl ; ++i) {
      for (int j=0 ; j<Nc ; ++j) {
        M_toFill(static_cast<int>(indexes_lines(i)), static_cast<int>(indexes_columns(j))) = M_toCopy(i, j);
      }
    }
    """
    weave.inline(wrapping_code,
             ['M_toFill', 'M_toCopy', 'indexes_lines', 'indexes_columns'],
             type_converters = converters.blitz,
             compiler = 'gcc',
             extra_compile_args = ['-O3', '-pthread', '-w'])

def MgPreconditionerComputationPerCube(cube, list_cubes_with_neighbors, list_Z_tmp, pathToReadFrom, cubeNumber_to_chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C):
    Z_local = eye(cube.N_RWG_src, cube.N_RWG_src).astype(Z_TMP_ELEM_TYPE)
    src_edges_numbers_local_src_edges_numbers = {}
    index = 0
    for RWGnumber in cube.testSrc_RWGsNumbers:
        src_edges_numbers_local_src_edges_numbers[RWGnumber] = index
        index += 1
    ## construction of the local matrix to be inverted
    for i in range(cube.N_neighbors):
        neighborCubeNumber = cube.cubeNeighborsIndexes[int(i)]
        chunkNumber = cubeNumber_to_chunkNumber[neighborCubeNumber]
        pathToReadFromChunk = os.path.join( pathToReadFrom, "chunk" + str(chunkNumber) )
        if not list_cubes_with_neighbors.has_key(neighborCubeNumber):
            neighborCube = CubeClass()
            neighborCube.setIntDoubleArraysFromFile(pathToReadFromChunk, neighborCubeNumber)
            list_cubes_with_neighbors[neighborCubeNumber] = copy.copy(neighborCube)
            list_Z_tmp[neighborCubeNumber] = read_Z_perCube_fromFile(pathToReadFromChunk, neighborCubeNumber, neighborCube, Z_TMP_ELEM_TYPE)
        else:
            neighborCube = list_cubes_with_neighbors[neighborCubeNumber]
        Z_tmp = list_Z_tmp[neighborCubeNumber]
        if i==0:
            ## we first fill in the first lines of Z_local
            indexes_lines = arange(Z_tmp.shape[0]).astype('i')
            indexes_columns = arange(Z_tmp.shape[1]).astype('i')
            assignValuesToMatrix(Z_local, Z_tmp, indexes_lines, indexes_columns)
        else:
            ## we then fill in the remaining lines
            Z_local_lines_indexes = zeros(neighborCube.N_RWG_test, 'i')
            index = 0
            for RWGnumber in neighborCube.testSrc_RWGsNumbers[:neighborCube.N_RWG_test]:
                Z_local_lines_indexes[index] = src_edges_numbers_local_src_edges_numbers[RWGnumber]
                index += 1
            columnsOfNeighborCubeToBeConsidered = findEdgesInRadiusAroundCube(neighborCube.localTestSrcRWGNumber_nodes, neighborCube.nodesCoord, cube.rCubeCenter, a * 1.5)
            srcEdgesNumbersOfNeighborCubeToBeConsidered = compress(columnsOfNeighborCubeToBeConsidered, neighborCube.testSrc_RWGsNumbers, axis=0)
            Z_local_columns_indexes = 0 * srcEdgesNumbersOfNeighborCubeToBeConsidered
            index = 0
            for RWGnumber in srcEdgesNumbersOfNeighborCubeToBeConsidered:
                Z_local_columns_indexes[index] = src_edges_numbers_local_src_edges_numbers[RWGnumber]
                index += 1
            Z_tmp2 = compress(columnsOfNeighborCubeToBeConsidered, Z_tmp, axis=1)
            assignValuesToMatrix(Z_local, Z_tmp2, Z_local_lines_indexes, Z_local_columns_indexes)
    Z_local_2 = take(Z_local, cube.localEdgesInRadiusAroundCube, axis=0).astype('D')
    src_edges_numbers_2 = take(cube.testSrc_RWGsNumbers, cube.localEdgesInRadiusAroundCube, axis=0)
    #Y_CFIE_near_local = linalg.pinv(Z_local_2)
    #Y_CFIE_near_local = computeMyPinv(Z_local_2)
    Y_CFIE_near_local = computeMyPinvCC(Z_local_2, LIB_G2C)
    # we take the first N_RWG_test lines
    Mg_tmp = take(Y_CFIE_near_local, arange(cube.N_RWG_test), axis=0)
    # we now "flatten" the preconditioner, in order to make a sparse matrix
    return (Mg_tmp.flat[:]).astype(ELEM_TYPE), src_edges_numbers_2

def MgPreconditionerComputation(target_mesh, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom):
    """This function computes the Left Frobenius preconditioner"""
    #C = target_mesh.cubes_centroids.shape[0]
    #cubeNumber_to_chunkNumber = zeros(C, 'i')
    #N_precond = numberOfElemsInPrecond(0, C, R_NORM_TYPE_1)
    #Mg, Mg_pq_array = zeros(N_precond, ELEM_TYPE), zeros((N_precond, 2), 'i')
    #startIndex2 = 0
    #for i in range(C):
        #Mg_tmp, Mg_pq_array_tmp = MgPreconditionerComputationPerCube(pathToReadFrom, cubeNumber_to_chunkNumber, target_mesh, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE)
        #Mg[startIndex2:startIndex2 + Mg_tmp.shape[0]], Mg_pq_array[startIndex2:startIndex2 + Mg_tmp.shape[0], :] = Mg_tmp, Mg_pq_array_tmp
        #startIndex2 += Mg_tmp.shape[0]
        #sys.stdout.write("\r" + "Left Frobenius preconditioner computation. Percentage = %.4s" %str(i*100./C))
        #sys.stdout.flush()
    #print
    #return Mg, Mg_pq_array

def chunk_of_Mg_CSR(cubesNumbers, chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, cubeNumber_to_chunkNumber):
    """same as for chunk_of_Z_near_CRS, but for the preconditioner"""
    pathToReadChunkFrom = os.path.join(pathToReadFrom, "chunk" + str(chunkNumber))
    list_cubes = compute_list_cubes(cubesNumbers, pathToReadChunkFrom)
    list_cubes_with_neighbors = copy.copy(list_cubes)
    N_RWG = 0
    for cubeNumber, cube in list_cubes.iteritems():
        Nl, Nc = cube.N_RWG_test, cube.N_RWG_src
        N_RWG += Nl
    test_RWG_numbers = zeros(N_RWG, 'i')
    # number of elements in the preconditioner chunk
    N_precond, N_ColumnsPerCube, list_Z_tmp = numberOfElemsInPrecond(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber, R_NORM_TYPE_1, Z_TMP_ELEM_TYPE)
    Mg = zeros(N_precond, ELEM_TYPE)
    # for the q_array, each src function for all the testing functions of a cube appears only once
    # instead of once per testing function. This allows a dramatic reduction in q_array.size
    q_array = zeros(sum(N_ColumnsPerCube), 'i') # column indexes
    rowIndexToColumnIndexes = zeros((N_RWG, 2), 'i') # start and end indexes
    startIndex, startIndexInRWGNumbers, startIndexInQArray = 0, 0, 0
    indexN_ColumnsPerCube, index_in_rowIndexToColumnIndexes = 0, 0
    for cubeNumber, cube in list_cubes.iteritems():
        # finding the RWGs numbers
        chunkNumber = cubeNumber_to_chunkNumber[cubeNumber]
        test_RWG_numbers[startIndexInRWGNumbers:startIndexInRWGNumbers + cube.N_RWG_test] = cube.testSrc_RWGsNumbers[:cube.N_RWG_test]
        startIndexInRWGNumbers += cube.N_RWG_test
        # computing the sparse matrix
        Mg_tmp, Mg_q_array = MgPreconditionerComputationPerCube(cube, list_cubes_with_neighbors, list_Z_tmp, pathToReadFrom, cubeNumber_to_chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C)
        Mg[startIndex:startIndex + Mg_tmp.shape[0]] = Mg_tmp
        startIndex += Mg_tmp.shape[0] # index update
        # q_array gives the column indexes.
        q_array[startIndexInQArray:startIndexInQArray + N_ColumnsPerCube[indexN_ColumnsPerCube]] = Mg_q_array
        # now we have to construct rowIndexToColumnIndexes
        indInf, indSup = index_in_rowIndexToColumnIndexes, index_in_rowIndexToColumnIndexes+cube.N_RWG_test
        rowIndexToColumnIndexes[indInf:indSup, 0] = startIndexInQArray
        rowIndexToColumnIndexes[indInf:indSup, 1] = startIndexInQArray + N_ColumnsPerCube[indexN_ColumnsPerCube]
        index_in_rowIndexToColumnIndexes = indSup
        # update startIndexInQArray
        startIndexInQArray += N_ColumnsPerCube[indexN_ColumnsPerCube]
        indexN_ColumnsPerCube += 1
    return Mg, q_array, rowIndexToColumnIndexes, test_RWG_numbers

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
            cube.setIntDoubleArraysFromFile(pathToReadCubeFrom, localCube)
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
    ## we find the dimensions of the cubes to exchange
    dimensionsOfCubesToSend = []
    for L in listCubesNumbersToSend:
        dimensionsOfCubesToSend.append([])
        for cubeNumber in L:
            chunkNumber = ZnearCubeNumber_to_chunkNumber[cubeNumber]
            pathToReadCubeFrom = os.path.join(pathToReadFrom, "chunk" + str(chunkNumber))
            cube = CubeClass()
            cube.setIntDoubleArraysFromFile(pathToReadCubeFrom, cubeNumber)
            dimensionsOfCubesToSend[-1].append([cube.N_RWG_test, cube.N_RWG_src])
    ## now we write the data to be exchanged to disk
    for i in range(num_proc):
        if not (my_id==i):
            writeASCIIBlitzArrayToDisk(array(listCubesNumbersToReceive[i]).astype('i'), os.path.join(pathToReadFrom,  "CubesNumbersToReceiveFromP" + str(i) + ".txt"))
            writeASCIIBlitzArrayToDisk(array(listCubesNumbersToSend[i]).astype('i'), os.path.join(pathToReadFrom, "CubesNumbersToSendToP" + str(i) + ".txt"))
            writeASCIIBlitzArrayToDisk(array(dimensionsOfCubesToSend[i]).astype('i'), os.path.join(pathToReadFrom, "dimensionsOfCubesToSendToP" + str(i) + ".txt"))
            writeASCIIBlitzArrayToDisk(array(listChunkNumbersToReceive[i]).astype('i'), os.path.join(pathToReadFrom, "ChunkNumbersToReceiveFromP" + str(i) + ".txt"))
            writeASCIIBlitzArrayToDisk(array(listChunkNumbersToSend[i]).astype('i'), os.path.join(pathToReadFrom, "ChunkNumbersToSendToP" + str(i) + ".txt"))
    #MPI.COMM_WORLD.Barrier()
    ## finally we write the format of the Near Field matrix elements
    NBytes = 8
    if Z_TMP_ELEM_TYPE=='D':
        NBytes = 16
        print "16 Bytes not supported yet in data transfer in communicateZnearBlocks. Exiting...."
        sys.exit(1)
    writeScalarToDisk(NBytes, os.path.join(pathToReadFrom, "itemsize.txt"))
    MPI.COMM_WORLD.Barrier()

def Mg_CSR(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, pathToSaveTo):
    """this function computes Mg by slices and stores them on the disk."""
    NAME = "Mg_LeftFrob"
    chunkNumbers = processNumber_to_ChunksNumbers[my_id]
    index, percentage = 0, 0
    for chunkNumber in chunkNumbers:
        if my_id==0:
            newPercentage = int(index * 100.0/len(processNumber_to_ChunksNumbers[my_id]))
            if (newPercentage - percentage)>=5:
                print "Process", my_id, ": computing SAI precond chunk.", newPercentage, "% completed"
                sys.stdout.flush()
                percentage = newPercentage
        pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
        os.mkdir(pathToSaveToChunk)
        cubesNumbers = chunkNumber_to_cubesNumbers[chunkNumber]
        Mg, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers = chunk_of_Mg_CSR(cubesNumbers, chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, cubeNumber_to_chunkNumber)
        writeToDisk_chunk_of_Z_sparse(pathToSaveTo, NAME, Mg, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers, chunkNumber)
        index += 1
    # we write the chunks numbers of the process
    writeASCIIBlitzArrayToDisk(array(chunkNumbers).astype('i'), os.path.join(pathToSaveTo, 'chunkNumbers.txt'))
