from mpi4py import MPI
import sys, os, copy
from scipy import zeros, array, compress, arange, take, product, dot, sort, argsort, ones, weave
from scipy.weave import converters
from scipy import prod, rand, eye, put, transpose #, linalg
from myPseudoInv import computeMyPinvCC #, computeMyPinv
from FMM_Znear import read_Z_perCube_fromFile, writeToDisk_chunk_of_Z_sparse
from ReadWriteBlitzArray import writeBlitzArrayToDisk, readBlitzArrayFromDisk, writeScalarToDisk, writeASCIIBlitzArrayToDisk
from meshClass import CubeClass


def compute_list_cubes(cubesNumbers, pathToReadCubesFrom, Z_TMP_ELEM_TYPE):
    """returns a list of cubes"""
    list_cubes = {}
    list_Z_tmp = {}
    for i in range(len(cubesNumbers)):
        cubeNumber = cubesNumbers[i]
        cube = CubeClass()
        cube.setIntArraysFromFile(pathToReadCubesFrom, cubeNumber)
        list_cubes[cubeNumber] = copy.copy(cube)
        list_Z_tmp[cubeNumber] = read_Z_perCube_fromFile(pathToReadCubesFrom, cubeNumber, cube, Z_TMP_ELEM_TYPE)
    return list_cubes, list_Z_tmp

def computePreconditionerColumnsPerCube(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber):
    """this function computes the number of columns per cube for the Frobenius preconditioner"""
    NumberOfColumnsPerCube = zeros(len(list_cubes), 'i')
    i = 0
    for cubeNumber, cube in list_cubes.iteritems():
        cube.localEdgesInRadiusAroundCube = (compress(cube.isEdgeInCartesianRadius==1, arange(len(cube.isEdgeInCartesianRadius)))).astype('i')
        chunkNumber = cubeNumber_to_chunkNumber[cubeNumber]
        pathToReadFromChunk = os.path.join( pathToReadFrom, "chunk" + str(chunkNumber) )
        NumberOfColumnsPerCube[i] = sum(cube.isEdgeInCartesianRadius)
        i += 1
    return NumberOfColumnsPerCube

def numberOfElemsInPrecond(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber):
    """computes the number of elements in the preconditioner"""
    preconditionerColumnsPerCube = computePreconditionerColumnsPerCube(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber)
    N_precond = 0
    i = 0
    for cubeNumber, cube in list_cubes.iteritems():
        N_precond += preconditionerColumnsPerCube[i] * cube.N_RWG_test
        i += 1
    return N_precond, preconditionerColumnsPerCube

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

def MgPreconditionerComputationPerCube(cube, list_cubes_with_neighbors, list_Z_tmp, pathToReadFrom, cubeNumber_to_chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C):
    Z_local = eye(cube.N_RWG_src, cube.N_RWG_src).astype(Z_TMP_ELEM_TYPE)
    src_edges_numbers_local_src_edges_numbers = {}
    index = 0
    for RWGnumber in cube.testSrc_RWGsNumbers:
        src_edges_numbers_local_src_edges_numbers[RWGnumber] = index
        index += 1
    ## construction of the local matrix to be inverted
    for i in range(cube.N_neighbors):
        neighborCubeNumber = cube.cubeNeighborsIndexes[int(i)]
        if not list_cubes_with_neighbors.has_key(neighborCubeNumber):
            chunkNumber = cubeNumber_to_chunkNumber[neighborCubeNumber]
            pathToReadFromChunk = os.path.join( pathToReadFrom, "chunk" + str(chunkNumber) )
            neighborCube = CubeClass()
            neighborCube.setIntArraysFromFile(pathToReadFromChunk, neighborCubeNumber)
            list_cubes_with_neighbors[neighborCubeNumber] = copy.copy(neighborCube)
            list_Z_tmp[neighborCubeNumber] = read_Z_perCube_fromFile(pathToReadFromChunk, neighborCubeNumber, neighborCube, Z_TMP_ELEM_TYPE)

    set_cubeNeighborsIndexes = set(cube.cubeNeighborsIndexes)

    for i in range(cube.N_neighbors):
        neighborCubeNumber = cube.cubeNeighborsIndexes[int(i)]
        neighborCube = list_cubes_with_neighbors[neighborCubeNumber]
        Z_neighbor = list_Z_tmp[neighborCubeNumber]
        if i==0:
            ## we first fill in the first lines of Z_local
            indexes_lines = arange(Z_neighbor.shape[0]).astype('i')
            indexes_columns = arange(Z_neighbor.shape[1]).astype('i')
            assignValuesToMatrix(Z_local, Z_neighbor, indexes_lines, indexes_columns)
        else:
            ## we then fill in the remaining lines
            Z_local_lines_indexes = zeros(neighborCube.N_RWG_test, 'i')
            index = 0
            for RWGnumber in neighborCube.testSrc_RWGsNumbers[:neighborCube.N_RWG_test]:
                Z_local_lines_indexes[index] = src_edges_numbers_local_src_edges_numbers[RWGnumber]
                index += 1

            common_neighborsNumbers = [val for val in neighborCube.cubeNeighborsIndexes if val in set_cubeNeighborsIndexes]
            srcEdgesNumbersOfNeighborCubeToBeConsidered_tmp = []
            for common_neighbor in common_neighborsNumbers:
                common_neighborCube = list_cubes_with_neighbors[common_neighbor]
                testRWGs = common_neighborCube.testSrc_RWGsNumbers[:common_neighborCube.N_RWG_test]
                for RWG in testRWGs:
                    srcEdgesNumbersOfNeighborCubeToBeConsidered_tmp.append(RWG)
            set_srcEdgesNumbersOfNeighborCubeToBeConsidered = set(srcEdgesNumbersOfNeighborCubeToBeConsidered_tmp)
            columnsOfNeighborCubeToBeConsidered = array([val for val in range(neighborCube.N_RWG_src) if neighborCube.testSrc_RWGsNumbers[val] in set_srcEdgesNumbersOfNeighborCubeToBeConsidered], 'i')

            Z_local_columns_indexes = zeros(len(srcEdgesNumbersOfNeighborCubeToBeConsidered_tmp), 'i')
            index = 0
            for RWGnumber in srcEdgesNumbersOfNeighborCubeToBeConsidered_tmp:
                Z_local_columns_indexes[index] = src_edges_numbers_local_src_edges_numbers[RWGnumber]
                index += 1
            Z_tmp2 = take(Z_neighbor, columnsOfNeighborCubeToBeConsidered, axis=1)
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

def chunk_of_Mg_CSR(cubesNumbers, chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, cubeNumber_to_chunkNumber):
    """same as for chunk_of_Z_near_CRS, but for the preconditioner"""
    pathToReadChunkFrom = os.path.join(pathToReadFrom, "chunk" + str(chunkNumber))
    list_cubes, list_Z_tmp = compute_list_cubes(cubesNumbers, pathToReadChunkFrom, Z_TMP_ELEM_TYPE)
    list_cubes_with_neighbors = copy.copy(list_cubes)
    N_RWG = 0
    for cubeNumber, cube in list_cubes.iteritems():
        Nl, Nc = cube.N_RWG_test, cube.N_RWG_src
        N_RWG += Nl
    test_RWG_numbers = zeros(N_RWG, 'i')
    # number of elements in the preconditioner chunk
    N_precond, N_ColumnsPerCube = numberOfElemsInPrecond(list_cubes, pathToReadFrom, cubeNumber_to_chunkNumber)
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
        Mg_tmp, Mg_q_array = MgPreconditionerComputationPerCube(cube, list_cubes_with_neighbors, list_Z_tmp, pathToReadFrom, cubeNumber_to_chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C)
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

def Mg_CSR(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, pathToSaveTo):
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
        Mg, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers = chunk_of_Mg_CSR(cubesNumbers, chunkNumber, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, cubeNumber_to_chunkNumber)
        writeToDisk_chunk_of_Z_sparse(pathToSaveTo, NAME, Mg, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers, chunkNumber)
        index += 1
    # we write the chunks numbers of the process
    writeASCIIBlitzArrayToDisk(array(chunkNumbers).astype('i'), os.path.join(pathToSaveTo, 'chunkNumbers.txt'))

