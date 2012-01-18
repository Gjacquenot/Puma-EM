#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <blitz/array.h>
#include <mpi.h>

using namespace std;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"
#include "Z_EJ_Z_HJ.h"
#include "scatter_mesh_per_cube.h"

void compute_cube_local_arrays(int & N_RWG_test,
                               int & N_RWG_src,
                               int & N_neighbors,
                               int & N_nodes,
                               int & S,
                               blitz::Array<int, 1>& test_RWGsNumbers,
                               blitz::Array<int, 1>& src_RWGsNumbers,
                               blitz::Array<int, 1>& localTestRWGNumber_CFIE_OK,
                               blitz::Array<int, 1>& localSrcRWGNumber_M_CURRENT_OK,
                               blitz::Array<int, 1>& neighbors_cubes,
                               blitz::Array<int, 2>& localTestSrcRWGNumber_signedTriangles,
                               blitz::Array<int, 2>& localTestSrcRWGNumber_nodes,
                               blitz::Array<double, 2>& nodesCoord,
                               blitz::Array<double, 1>& rCubeCenter,
                               const blitz::Array<int, 1>& cubeIntArrays,
                               const blitz::Array<double, 1>& cubeDoubleArrays)
{
  // first the ints
  N_RWG_test = cubeIntArrays(0);
  N_RWG_src = cubeIntArrays(1);
  N_neighbors = cubeIntArrays(2);
  N_nodes = cubeIntArrays(3);
  S = cubeIntArrays(4);
  
  // test_RWGsNumbers, src_RWGsNumbers
  int startIndex = 5;
  int stopIndex = startIndex + N_RWG_src;
  test_RWGsNumbers.resize(N_RWG_test);
  src_RWGsNumbers.resize(N_RWG_src);
  for (int i=0; i<N_RWG_test; i++) test_RWGsNumbers(i) = i; // cubeIntArrays(startIndex + i);
  for (int i=0; i<N_RWG_src; i++) src_RWGsNumbers(i) = i; //cubeIntArrays(startIndex + i);

  // signed triangles
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src * 2;
  localTestSrcRWGNumber_signedTriangles.resize(N_RWG_src, 2);
  for (int i=0; i<N_RWG_src; i++) {
    localTestSrcRWGNumber_signedTriangles(i, 0) = cubeIntArrays(startIndex + i*2);
    localTestSrcRWGNumber_signedTriangles(i, 1) = cubeIntArrays(startIndex + i*2 + 1);
  }

  // the local nodes
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src * 4;
  localTestSrcRWGNumber_nodes.resize(N_RWG_src, 4);
  for (int i=0; i<N_RWG_src; i++) {
    localTestSrcRWGNumber_nodes(i, 0) = cubeIntArrays(startIndex + i*4);
    localTestSrcRWGNumber_nodes(i, 1) = cubeIntArrays(startIndex + i*4 + 1);
    localTestSrcRWGNumber_nodes(i, 2) = cubeIntArrays(startIndex + i*4 + 2);
    localTestSrcRWGNumber_nodes(i, 3) = cubeIntArrays(startIndex + i*4 + 3);
  }

  // the local CFIE OK  
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_test;
  localTestRWGNumber_CFIE_OK.resize(N_RWG_test);
  for (int i=0; i<N_RWG_test; i++) localTestRWGNumber_CFIE_OK(i) = cubeIntArrays(startIndex + i);

  // local M CURRENTS OK
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src;
  localSrcRWGNumber_M_CURRENT_OK.resize(N_RWG_src);
  for (int i=0; i<N_RWG_src; i++) localSrcRWGNumber_M_CURRENT_OK(i) = cubeIntArrays(startIndex + i);

  // neighbors cubes
  startIndex = stopIndex;
  stopIndex = startIndex + N_neighbors;
  neighbors_cubes.resize(N_neighbors);
  for (int i=0; i<N_neighbors; i++) neighbors_cubes(i) = cubeIntArrays(startIndex + i);

  // now on to the double arrays
  startIndex = 0;
  nodesCoord.resize(N_nodes, 3);
  for (int i=0; i<N_nodes; i++) {
    nodesCoord(i, 0) = cubeDoubleArrays(3*i);
    nodesCoord(i, 1) = cubeDoubleArrays(3*i + 1);
    nodesCoord(i, 2) = cubeDoubleArrays(3*i + 2);
  }
  
  // rCubeCenter
  rCubeCenter(0) = cubeDoubleArrays(3*N_nodes);
  rCubeCenter(1) = cubeDoubleArrays(3*N_nodes + 1);
  rCubeCenter(2) = cubeDoubleArrays(3*N_nodes + 2);
}

void findRWGsInRadiusAroundCube(blitz::Array<int, 1>& isRWGInCartesianRadius,
                                const int N_RWG_src,
                                const blitz::Array<int, 2>& RWGNumber_nodes,
                                const blitz::Array<double, 2>& nodesCoord,
                                const blitz::Array<double, 1>& rCubeCenter,
                                const double R_NORM_TYPE_1)
{
  // this function yields a 1 to an edge number if it is within
  // a certain 1-Norm around the center of the cube to which it belongs
  isRWGInCartesianRadius.resize(N_RWG_src);
  isRWGInCartesianRadius = 1;
  for (int i=0 ; i<N_RWG_src ; ++i) {
    const int node1 = RWGNumber_nodes(i, 1), node2 = RWGNumber_nodes(i, 2);
    for (int j=0 ; j<3 ; ++j) {
      const double center_coord = (nodesCoord(node1, j) + nodesCoord(node2, j))/2.0;
      if (std::abs(center_coord - rCubeCenter(j)) > R_NORM_TYPE_1) {
        isRWGInCartesianRadius(i) = 0;
        break;
      }
    }
  }
}


void compute_cubeSmallIntArray(blitz::Array<int, 1>& cubeSmallIntArray,
                               const blitz::Array<int, 1>& cubeIntArrays,
                               const blitz::Array<int, 1>& neighbors_cubes,
                               const blitz::Array<int, 2>& localTestSrcRWGNumber_nodes,
                               const blitz::Array<double, 2>& nodesCoord, 
                               const blitz::Array<double, 1>& rCubeCenter,
                               const double R_NORM_TYPE_1)
{
  const int N_RWG_src = cubeIntArrays(1);
  const int N_neighbors = cubeIntArrays(2);
  // we will only save the 5 first elements, the RWG numbers, the nodes numbers, and the cube neighbors
  const int size_small_array = 5 + N_RWG_src + N_RWG_src + N_neighbors;
  cubeSmallIntArray.resize(size_small_array);
  cubeSmallIntArray(0) = cubeIntArrays(0);
  cubeSmallIntArray(1) = cubeIntArrays(1);
  cubeSmallIntArray(2) = cubeIntArrays(2);
  cubeSmallIntArray(3) = cubeIntArrays(3);
  cubeSmallIntArray(4) = cubeIntArrays(4);

  // src_RWGsNumbers
  int startIndex = 5;
  int stopIndex = startIndex + N_RWG_src;
  for (int i=0; i<N_RWG_src; i++) cubeSmallIntArray(startIndex + i) = cubeIntArrays(startIndex + i);

  blitz::Array<int, 1> isRWGInCartesianRadius;
  findRWGsInRadiusAroundCube(isRWGInCartesianRadius, N_RWG_src, localTestSrcRWGNumber_nodes, nodesCoord, rCubeCenter, R_NORM_TYPE_1);
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src;
  for (int i=0; i<N_RWG_src; i++) cubeSmallIntArray(startIndex + i) = isRWGInCartesianRadius(i);

  // neighbors cubes
  startIndex = stopIndex;
  stopIndex = startIndex + N_neighbors;
  for (int i=0; i<N_neighbors; i++) cubeSmallIntArray(startIndex + i) = neighbors_cubes(i);
}


int main(int argc, char* argv[]) {

  MPI::Init();
  int ierror;
  const int num_procs = MPI::COMM_WORLD.Get_size();
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0;
  MPI_Status status;

  string simuDir = ".";
  if ( argc > 2 ) {
     if( string(argv[1]) == "--simudir" ) simuDir = argv[2];
  }

  // general variables
  const string SIMU_DIR = simuDir;
  const string TMP = SIMU_DIR + "/tmp" + intToString(my_id);
  const string Z_TMP_DATA_PATH = TMP + "/Z_tmp/";
  const string OCTTREE_DATA_PATH = TMP + "/octtree_data/";
  string filename;

  // reading the simulation variables
  int TDS_APPROX, MOM_FULL_PRECISION;
  readIntFromASCIIFile(OCTTREE_DATA_PATH  + "TDS_APPROX.txt", TDS_APPROX);
  readIntFromASCIIFile(OCTTREE_DATA_PATH  + "MOM_FULL_PRECISION.txt", MOM_FULL_PRECISION);
  blitz::Array<std::complex<double>, 1> CFIEcoeffs;
  readComplexDoubleBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "CFIEcoeffs.txt", CFIEcoeffs);
  float w;
  readFloatFromASCIIFile(OCTTREE_DATA_PATH  + "w.txt", w);
  std::complex<float> eps_r, mu_r, Z_s;
  readComplexFloatFromASCIIFile(OCTTREE_DATA_PATH  + "eps_r.txt", eps_r);
  readComplexFloatFromASCIIFile(OCTTREE_DATA_PATH  + "mu_r.txt", mu_r);
  readComplexFloatFromASCIIFile(OCTTREE_DATA_PATH  + "Z_s.txt", Z_s);
  double R_NORM_TYPE_1;
  readDoubleFromASCIIFile(OCTTREE_DATA_PATH  + "leaf_side_length.txt", R_NORM_TYPE_1);
  
  // reading mesh data
  int N_local_Chunks, N_local_cubes;

  filename = Z_TMP_DATA_PATH + "N_local_Chunks.txt";
  readIntFromASCIIFile(filename, N_local_Chunks);

  filename = Z_TMP_DATA_PATH + "N_local_cubes.txt";
  readIntFromASCIIFile(filename, N_local_cubes);
  
  blitz::Array<int, 1> local_ChunksNumbers(N_local_Chunks);
  blitz::Array<int, 1> local_chunkNumber_N_cubesNumbers(N_local_Chunks);
  blitz::Array<int, 1> local_chunkNumber_to_cubesNumbers(N_local_cubes);

  filename = Z_TMP_DATA_PATH + "local_ChunksNumbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_ChunksNumbers);
  
  filename = Z_TMP_DATA_PATH + "local_chunkNumber_N_cubesNumbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_chunkNumber_N_cubesNumbers);
  
  filename = Z_TMP_DATA_PATH + "local_chunkNumber_to_cubesNumbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_chunkNumber_to_cubesNumbers);

  // creating a cube_to_chunk correspondance. Needed for reading and saving the files.
  blitz::Array<int, 1> local_cubeNumber_to_chunkNumbers(N_local_cubes);
  int startIndex = 0;
  for (int i=0; i<N_local_Chunks; i++) {
    const int chunkNumber = local_ChunksNumbers(i);
    const int N_cubes_in_chunk = local_chunkNumber_N_cubesNumbers(i);
    for (int j=0; j<N_cubes_in_chunk; j++) local_cubeNumber_to_chunkNumbers(startIndex + j) = chunkNumber;
    startIndex += N_cubes_in_chunk;
  }

  blitz::Array<blitz::Array<int, 1>, 1> allCubeIntArrays;
  blitz::Array<blitz::Array<double, 1>, 1> allCubeDoubleArrays;
  scatter_mesh_per_cube(allCubeIntArrays, allCubeDoubleArrays, SIMU_DIR, local_ChunksNumbers, local_chunkNumber_N_cubesNumbers, local_chunkNumber_to_cubesNumbers, local_cubeNumber_to_chunkNumbers);

  // now we compute the Z near and write them to disk
  float percentage = 0.0;
  for (int i=0; i<N_local_cubes; i++) {
    if (my_id==master) {
      float newPercentage = i * 100.0/N_local_cubes;
      if ((newPercentage - percentage)>=10.0) {
        std::cout << "Process " <<  my_id << " : computing Z_CFIE_near chunk. " << std::floor(newPercentage) <<  " % completed." << endl;
        flush(std::cout);
        percentage = newPercentage;
      }
    }
    // computing the local arrays
    int N_RWG_test, N_RWG_src, N_neighbors, N_nodes, S;
    blitz::Array<int, 1> test_RWGsNumbers, src_RWGsNumbers;
    blitz::Array<int, 1> localTestRWGNumber_CFIE_OK, localSrcRWGNumber_M_CURRENT_OK, neighbors_cubes;
    blitz::Array<int, 2> localTestSrcRWGNumber_signedTriangles, localTestSrcRWGNumber_nodes;
    blitz::Array<double, 2> nodesCoord;
    blitz::Array<double, 1> rCubeCenter(3);
    compute_cube_local_arrays(N_RWG_test, N_RWG_src, N_neighbors, N_nodes, S, test_RWGsNumbers, src_RWGsNumbers, localTestRWGNumber_CFIE_OK, localSrcRWGNumber_M_CURRENT_OK, neighbors_cubes, localTestSrcRWGNumber_signedTriangles, localTestSrcRWGNumber_nodes, nodesCoord, rCubeCenter, allCubeIntArrays(i), allCubeDoubleArrays(i));

    // computing the Z_CFIE
    blitz::Array<std::complex<double>, 2> Z_CFIE_J(N_RWG_test, N_RWG_src), Z_CFIE_M(1, 1);
    Z_CFIE_J = 0.0;
    Z_CFIE_M = 0.0;
    localSrcRWGNumber_M_CURRENT_OK *= 0; // no dielectric in MLFMA yet
    const double signSurfObs = 1.0, signSurfSrc = 1.0; // no dielectric in MLFMA yet
    Z_CFIE_J_computation(Z_CFIE_J, Z_CFIE_M, CFIEcoeffs, signSurfObs, signSurfSrc, test_RWGsNumbers, src_RWGsNumbers, localTestRWGNumber_CFIE_OK, localSrcRWGNumber_M_CURRENT_OK, localTestSrcRWGNumber_signedTriangles, localTestSrcRWGNumber_nodes, nodesCoord, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION);

    // transforming and writing the matrix to the disk
    blitz::Array<std::complex<float>, 1> Z_CFIE_J_linear(N_RWG_test * N_RWG_src);
    for (int ii=0; ii<N_RWG_test; ii++) {
      for (int jj=0; jj<N_RWG_src; jj++) Z_CFIE_J_linear(ii*N_RWG_src + jj) = Z_CFIE_J(ii, jj);
    }
    // now we write the arrays to disk
    const int chunkNumber = local_cubeNumber_to_chunkNumbers(i);
    const int cubeNumber = local_chunkNumber_to_cubesNumbers(i);
    const string filenameIntArray = Z_TMP_DATA_PATH + "chunk" + intToString(chunkNumber) + "/" + intToString(cubeNumber) + "_IntArrays.txt";
    const string filenameZnear = Z_TMP_DATA_PATH + "chunk" + intToString(chunkNumber) + "/" + intToString(cubeNumber);
    writeComplexFloatBlitzArray1DToBinaryFile(filenameZnear, Z_CFIE_J_linear);
    // we now create a smaller version of the int array to save
    blitz::Array<int, 1> cubeSmallIntArray;
    compute_cubeSmallIntArray(cubeSmallIntArray, allCubeIntArrays(i), neighbors_cubes, localTestSrcRWGNumber_nodes, nodesCoord, rCubeCenter, R_NORM_TYPE_1);
    writeIntBlitzArray1DToBinaryFile(filenameIntArray, cubeSmallIntArray);
  }
  
  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;
}
