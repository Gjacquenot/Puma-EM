#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <blitz/array.h>
#include <mpi.h>
#include <map>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"
#include "GetMemUsage.h"
#include "./lapack/zgels_interface.h"

class CubeArrays {
  public:
    int number;
    int N_RWG_test;
    int N_RWG_src;
    int N_neighbors;
    std::vector<int> testSrc_RWGsNumbers;
    std::vector<int> isEdgeInCartesianRadius;
    std::vector<int> neighborsIndexes;
    blitz::Array<std::complex<float>, 2> Z_CFIE_J;
    // constructors
    CubeArrays(void){};
    CubeArrays(const int /*cubeNumber*/, const string /*pathToReadFrom*/); 
    void copyCubeArrays (const CubeArrays& cubeArraysToCopy);
    CubeArrays(const CubeArrays&); // copy constructor
    CubeArrays& operator=(const CubeArrays&); // copy assignment operator
    ~CubeArrays();
};

CubeArrays::CubeArrays(const int cubeNumber, const string pathToReadFrom)
{
  number = cubeNumber;
  const string pathToCubeIntArrays = pathToReadFrom + intToString(cubeNumber) + "_IntArrays.txt";  
  const string pathToCube_Z = pathToReadFrom + intToString(cubeNumber);
  // reading cubeIntArrays (cube information)
  blitz::Array<int, 1> cubeIntArrays;
  blitz::ifstream ifs(pathToCubeIntArrays.c_str(), blitz::ios::binary);
  ifs.seekg (0, blitz::ios::end);
  int length = ifs.tellg();
  ifs.close();
  int N_cubeIntArrays = length/4;
  cubeIntArrays.resize(N_cubeIntArrays);
  readIntBlitzArray1DFromBinaryFile(pathToCubeIntArrays, cubeIntArrays);

  N_RWG_test = cubeIntArrays(0);
  N_RWG_src = cubeIntArrays(1);
  N_neighbors = cubeIntArrays(2);

  int startIndex = 5, stopIndex = startIndex + N_RWG_src;
  testSrc_RWGsNumbers.reserve(N_RWG_src);
  for (int index=startIndex; index<stopIndex; index++) testSrc_RWGsNumbers.push_back(cubeIntArrays(index));

  isEdgeInCartesianRadius.reserve(N_RWG_src);
  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src;
  for (int index=startIndex; index<stopIndex; index++) isEdgeInCartesianRadius.push_back(cubeIntArrays(index));
    
  neighborsIndexes.reserve(N_neighbors);
  startIndex = stopIndex;
  stopIndex = startIndex + N_neighbors;
  for (int index=startIndex; index<stopIndex; index++) neighborsIndexes.push_back(cubeIntArrays(index));

   blitz::Array<std::complex<float>, 1> Z_CFIE_J_linear(N_RWG_test * N_RWG_src);
   readComplexFloatBlitzArray1DFromBinaryFile(pathToCube_Z, Z_CFIE_J_linear);
   Z_CFIE_J.resize(N_RWG_test, N_RWG_src);
   for (int ii=0; ii<N_RWG_test; ii++) {
     for (int jj=0; jj<N_RWG_src; jj++) Z_CFIE_J(ii, jj) = Z_CFIE_J_linear(ii*N_RWG_src + jj);
   }
}

void CubeArrays::copyCubeArrays(const CubeArrays& cubeArraysToCopy) // copy member function
{
  number = cubeArraysToCopy.number;
  N_RWG_test = cubeArraysToCopy.N_RWG_test;
  N_RWG_src = cubeArraysToCopy.N_RWG_src;
  N_neighbors = cubeArraysToCopy.N_neighbors;
  testSrc_RWGsNumbers.resize(cubeArraysToCopy.testSrc_RWGsNumbers.size());
  testSrc_RWGsNumbers = cubeArraysToCopy.testSrc_RWGsNumbers;
  isEdgeInCartesianRadius.resize(cubeArraysToCopy.isEdgeInCartesianRadius.size());
  isEdgeInCartesianRadius = cubeArraysToCopy.isEdgeInCartesianRadius;
  neighborsIndexes.resize(cubeArraysToCopy.neighborsIndexes.size());
  neighborsIndexes = cubeArraysToCopy.neighborsIndexes;
  Z_CFIE_J.resize(N_RWG_test, N_RWG_src);
  Z_CFIE_J = cubeArraysToCopy.Z_CFIE_J;
}

CubeArrays::CubeArrays(const CubeArrays& cubeArraysToCopy) // copy constructor
{
  copyCubeArrays(cubeArraysToCopy);
}

CubeArrays& CubeArrays::operator=(const CubeArrays& cubeArraysToCopy) { // copy assignment
  copyCubeArrays(cubeArraysToCopy);
  return *this;
}

CubeArrays::~CubeArrays() {
  testSrc_RWGsNumbers.clear();
  isEdgeInCartesianRadius.clear();
  neighborsIndexes.clear();
  Z_CFIE_J.free();
}

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
  typedef std::map<int, CubeArrays> CubeArraysMap;

  blitz::Array<int, 1> chunkNumbers, cubeNumber_to_chunkNumber;
  readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "chunkNumbers.txt", chunkNumbers);
  readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "cubeNumber_to_chunkNumber.txt", cubeNumber_to_chunkNumber);

  const int N_chunks = chunkNumbers.size();
  for (int i=0; i<N_chunks; i++) {
    const int chunk = chunkNumbers(i);
    blitz::Array<int, 1> cubesNumbers;
    readIntBlitzArray1DFromASCIIFile(SAI_PRECOND_DATA_PATH + "chunk" + intToString(chunk) + "cubesNumbers.txt", cubesNumbers);
    const int N_cubes = cubesNumbers.size();
    CubeArraysMap List_1, List_2;
    for (int j=0; j<N_cubes; j++) {
      const int cubeNumber = cubesNumbers(j);
      const string pathToCube = Z_TMP_DATA_PATH + "chunk" + intToString(chunk) + "/";
      List_1.insert(CubeArraysMap::value_type(cubeNumber, CubeArrays(cubeNumber, pathToCube)));
      
    }
  }
  
  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;
}
