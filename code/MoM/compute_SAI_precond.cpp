#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <blitz/array.h>
#include <mpi.h>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"
#include "GetMemUsage.h"
#include "./lapack/zgels_interface.h"

void computeLwork(int & lwork, const int M, const int N, const int nrhs) {
  // computes the space necessary for WORK to.....work!!
  // returns the dimension of the array WORK.
  // LWORK >= max( 1, MN + max( MN, NRHS ) ).
  // For optimal performance, LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
  // where MN = min(M,N) and NB is the optimum block size.
  int NB = 20;
  int mn = min(M, N);
  lwork = max( 1, mn + max( mn, nrhs )*NB );
}

void computeMyPinvCC(blitz::Array<std::complex<double>, 2>& Y, const blitz::Array<std::complex<double>, 2>& A) {
  /* this routine computes the pseudo inverse of A
  and is a wrapper to the fortran function zgels.f
  By the way, to understand this wrapping structure, you
  better check the comments at the beginning of zgels.f */
  int m = A.extent(0);
  int n = A.extent(1);
  int lda = m;
  int ldb = max(n, m);
  int nrhs = m;
  const int N = min(ldb, nrhs);
  blitz::Array<std::complex<double>, 2> B(ldb, nrhs, blitz::fortranArray), AA(ldb, nrhs, blitz::fortranArray);
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) AA(i, j) = A(i, j);
  }
  B = 0.0;
  for (int i=0 ; i<N ; ++i) B(i, i) = 1.0;
  int lwork;
  computeLwork(lwork, m, n, nrhs);
  blitz::Array<std::complex<double>, 1> work(lwork);
  work = 0.0;
  int info = 0;
  char trans = 'N';
  zgels(trans, m, n, nrhs, AA, lda, B, ldb, work, lwork, info);
  if (info==0) {
    if (m<=n) {
      Y.resize(B.extent(0), B.extent(1));
      for (int i=0; i<B.extent(0); i++) {
        for (int j=0; j<B.extent(1); j++) Y(i, j) = B(i, j);
      }
    }
    else {
      Y.resize(n, B.extent(1));
      for (int i=0; i<n; i++) {
        for (int j=0; j<B.extent(1); j++) Y(i, j) = B(i, j);
      }
    }
  }
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
  const string SAI_PRECOND_DATA_PATH = TMP + "/Mg_LeftFrob/";
  const string OCTTREE_DATA_PATH = TMP + "/octtree_data/";
  string filename;

  blitz::Array<int, 1> chunkNumbers, cubeNumber_to_chunkNumber;
  readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "chunkNumbers.txt", chunkNumbers);
  readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "cubeNumber_to_chunkNumber.txt", cubeNumber_to_chunkNumber);

  const int N_chunks = chunkNumbers.size();
  for (int i=0; i<N_chunks; i++) {
    const int chunk = chunkNumbers(i);
    blitz::Array<int, 1> cubesNumbers;
    readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "chunk" + intToString(chunk) + "cubesNumbers.txt", cubesNumbers);
    const int N_cubes = cubesNumbers.size();
    for (int j=0; j<N_cubes; j++) {
      const int cubeNumber = cubesNumbers(j);
      const string pathToCubeIntArrays = Z_TMP_DATA_PATH + "chunk" + intToString(chunk) + "/" + intToString(cubeNumber) + "_IntArrays.txt";
      // reading cubeIntArrays (cube information)
      blitz::Array<int, 1> cubeIntArrays;
      blitz::ifstream ifs(pathToCubeIntArrays.c_str(), blitz::ios::binary);
      ifs.seekg (0, blitz::ios::end);
      int length = ifs.tellg();
      ifs.close();
      int N_cubeIntArrays = length/4;
      cubeIntArrays.resize(N_cubeIntArrays);
      readIntBlitzArray1DFromBinaryFile(pathToCubeIntArrays, cubeIntArrays);
      int N_RWG_test = cubeIntArrays(0);
      int N_RWG_src = cubeIntArrays(1);
      int N_neighbors = cubeIntArrays(2);
      int N_nodes = cubeIntArrays(3);
      int startIndex = 5, stopIndex = startIndex + N_RWG_src;
      std::vector<int> testSrc_RWGsNumbers;
      testSrc_RWGsNumbers.reserve(N_RWG_src);
      for (int index=startIndex; index<stopIndex; index++) testSrc_RWGsNumbers.push_back(cubeIntArrays(index));

      std::vector<int> isEdgeInCartesianRadius;
      isEdgeInCartesianRadius.reserve(N_RWG_src);
      startIndex = stopIndex;
      stopIndex = startIndex + N_RWG_src;
      for (int index=startIndex; index<stopIndex; index++) isEdgeInCartesianRadius.push_back(cubeIntArrays(index));
    
      std::vector<int> cubeNeighborsIndexes;
      cubeNeighborsIndexes.reserve(N_neighbors);
      startIndex = stopIndex;
      stopIndex = startIndex + N_neighbors;
      for (int index=startIndex; index<stopIndex; index++) cubeNeighborsIndexes.push_back(cubeIntArrays(index));

      blitz::Array<std::complex<float>, 1> Z_CFIE_J_linear(N_RWG_test * N_RWG_src);
      blitz::Array<std::complex<float>, 2> Z_CFIE_J(N_RWG_test, N_RWG_src);
      const string pathToCube_Z = Z_TMP_DATA_PATH + "chunk" + intToString(chunk) + "/" + intToString(cubeNumber);
      readComplexFloatBlitzArray1DFromBinaryFile(pathToCube_Z, Z_CFIE_J_linear);
      for (int ii=0; ii<N_RWG_test; ii++) {
        for (int jj=0; jj<N_RWG_src; jj++) Z_CFIE_J(ii, jj) = Z_CFIE_J_linear(ii*N_RWG_src + jj);
      }
    }
  }
  
  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;
}
