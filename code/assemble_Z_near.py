import sys, os, cPickle, time, argparse, commands
from mpi4py import MPI
from scipy import array, zeros
from ReadWriteBlitzArray import readBlitzArrayFromDisk, writeBlitzArrayToDisk, writeScalarToDisk, writeASCIIBlitzArrayToDisk
from meshClass import CubeClass
import copy

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
        cube.setIntArraysFromFile(pathToReadCubeFrom, cubeNumber)
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
    if ( (MAX_BLOCK_SIZE<0.1) | (MAX_BLOCK_SIZE>10000.) ):
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


def assemble_Z_near(params_simu, simuDirName):
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    ELEM_TYPE, Z_TMP_ELEM_TYPE = 'F', 'F'
    # assembling of near interactions matrix
    pathToReadFrom, pathToSaveTo = os.path.join(tmpDirName, 'Z_tmp'), os.path.join(tmpDirName, 'Z_near')
    Z_nearCRS_Assembling(variables['processNumber_to_ChunksNumbers'], variables['chunkNumber_to_cubesNumbers'], params_simu.MAX_BLOCK_SIZE, variables['C'], ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()

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
        assemble_Z_near(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()
