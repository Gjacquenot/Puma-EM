from mpi4py import MPI
import os, sys, commands
from scipy import zeros, array
from Z_MoM import Z_MoM, Z_MoM_triangles_arraysFromCube
from ReadWriteBlitzArray import writeBlitzArrayToDisk, readBlitzArrayFromDisk, writeScalarToDisk, writeASCIIBlitzArrayToDisk
from meshClass import MeshClass, CubeClass
import copy

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

def chunk_of_Z_nearCRS_Computation(CFIE, cubesNumbers, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION, pathToSaveTo):
    """this function computes a chunk of the near non-diagonal part of the MoM matrix,
    but saves all the atomic elements on the disk. These elements will later on be used 
    by chunk_of_Z_nearCRS_Assembling and MgPrecondition"""
    pathToReadCubeFrom = pathToSaveTo
    list_cubes = compute_list_cubes(cubesNumbers, pathToReadCubeFrom)
    for cubeNumber, cube in list_cubes.iteritems():
        Z_CFIE_near_tmp = Z_nearPerCube(pathToSaveTo, cube, CFIE, cubeNumber, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, TDS_APPROX, Z_s, MOM_FULL_PRECISION)

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


# Assembling

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
    test_RWG_numbers = zeros(N_RWG, 'i')
    # number of elements in the Z_near chunk
    Z_CFIE_near = zeros(N_near, ELEM_TYPE) # matrix elements
    # for the src_RWG_numbers, each src function for all the testing functions of a cube appears only once
    # instead of once per testing function. This allows a dramatic reduction in src_RWG_numbers.shape.
    # See the help for this function.
    src_RWG_numbers = zeros(N_srcFromNeighbors, 'i') # column indexes
    rowIndexToColumnIndexes = zeros((N_RWG, 2), 'i') # start and end indexes
    startIndex, startIndexInTestRWGNumbers, startIndexInSrcRWGNumbers = 0, 0, 0
    index_in_rowIndexToColumnIndexes = 0
    for cubeNumber, cube in list_cubes.iteritems():
        # reading the sparse matrix
        Z_CFIE_near_tmp1 = read_Z_perCube_fromFile(pathToReadFrom, cubeNumber, cube, Z_TMP_ELEM_TYPE)
        Z_CFIE_near_tmp = array(Z_CFIE_near_tmp1.astype(ELEM_TYPE).flat)
        Z_CFIE_near[startIndex:startIndex + len(Z_CFIE_near_tmp)] = Z_CFIE_near_tmp.tolist()
        startIndex += Z_CFIE_near_tmp.shape[0] # index update
        # src_RWG_numbers gives the column indexes. This is the second column of pq_array
        src_RWG_numbers[startIndexInSrcRWGNumbers:startIndexInSrcRWGNumbers + cube.N_RWG_src] = cube.testSrc_RWGsNumbers
        # finding the RWGs numbers
        test_RWG_numbers[startIndexInTestRWGNumbers:startIndexInTestRWGNumbers + cube.N_RWG_test] = cube.testSrc_RWGsNumbers[:cube.N_RWG_test]
        startIndexInTestRWGNumbers += cube.N_RWG_test
        # now we have to construct rowIndexToColumnIndexes
        for j in range(cube.N_RWG_test):
            rowIndexToColumnIndexes[index_in_rowIndexToColumnIndexes, 0] = startIndexInSrcRWGNumbers
            rowIndexToColumnIndexes[index_in_rowIndexToColumnIndexes, 1] = startIndexInSrcRWGNumbers + cube.N_RWG_src
            index_in_rowIndexToColumnIndexes += 1
        # update startIndexInSrcArray
        startIndexInSrcRWGNumbers += cube.N_RWG_src
    return Z_CFIE_near, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers

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
        Z_CFIE_near, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers = chunk_of_Z_nearCRS_Assembling(cubesNumbers, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFromChunk)
        writeToDisk_chunk_of_Z_sparse(pathToSaveTo, NAME, Z_CFIE_near, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers, chunkNumber)
        del Z_CFIE_near, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers
        commands.getoutput("rm -rf " + os.path.join(pathToReadFromChunk))
    # we write the chunks numbers of the process
    writeASCIIBlitzArrayToDisk(array(chunkNumbers).astype('i'), os.path.join(pathToSaveTo, 'chunkNumbers.txt'))

def writeToDisk_chunk_of_Z_sparse(path, name, Z, src_RWG_numbers, rowIndexToColumnIndexes, test_RWG_numbers, chunkNumber):
    """this function writes to disk the chunks of Z sparse and the corresponding indexes arrays, each with a number"""
    chunkNumberString = str(chunkNumber)
    writeBlitzArrayToDisk(Z, os.path.join(path, name) + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(src_RWG_numbers, os.path.join(path, 'src_RWG_numbers') + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(rowIndexToColumnIndexes, os.path.join(path, 'rowIndexToColumnIndexes') + str(chunkNumber) + '.txt')
    writeBlitzArrayToDisk(test_RWG_numbers, os.path.join(path, 'test_RWG_numbers') + str(chunkNumber) + '.txt')
    # now we write the scalar values
    N_test_RWG_File = os.path.join(path, 'N_test_RWG') + str(chunkNumber) + '.txt'
    writeScalarToDisk(test_RWG_numbers.shape[0], N_test_RWG_File)
    N_near_File = os.path.join(path, 'N_near') + str(chunkNumber) + '.txt'
    writeScalarToDisk(Z.shape[0], N_near_File)
    N_src_RWG_File = os.path.join(path, 'N_src_RWG') + str(chunkNumber) + '.txt'
    writeScalarToDisk(src_RWG_numbers.shape[0], N_src_RWG_File)
