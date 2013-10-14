#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <blitz/array.h>
#include <mpi.h>
#include <map>
#include <set>

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

typedef std::map<int, CubeArrays> CubeArraysMap;
typedef std::map<int, CubeArrays>::const_iterator CubeArraysMapIterator;


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

void computeMyPinvCC(blitz::Array<std::complex<double>, 2>& Y, blitz::Array<std::complex<double>, 1>& A, int m, int n) {
  /* this routine computes the pseudo inverse of A
  and is a wrapper to the fortran function zgels.f
  By the way, to understand this wrapping structure, you
  better check the comments at the beginning of zgels.f */
  int lda = m;
  int ldb = max(n, m);
  int nrhs = m;
  const int N = min(ldb, nrhs);
  blitz::Array<std::complex<double>, 1> B(ldb * nrhs); // B(ldb, nrhs)
  B = 0.0;
  for (int i=0 ; i<N ; i++) B(i + i*ldb) = 1.0;
  int lwork;
  computeLwork(lwork, m, n, nrhs);
  blitz::Array<std::complex<double>, 1> work(lwork);
  work = 0.0;
  int info = 0;
  char trans = 'N';
  zgels2(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
  if (info==0) {
    if (m<=n) {
      Y.resize(ldb, nrhs);
      for (int j=0; j<nrhs; j++) {
        for (int i=0; i<ldb; i++) Y(i, j) = B(i + j*ldb);
      }
    }
    else {
      Y.resize(n, nrhs);
      for (int i=0; i<n; i++) {
        for (int j=0; j<nrhs; j++) Y(i, j) = B(i + j*n);
      }
    }
  }
}

void MgPreconditionerComputationPerCube(const int cubeNumber,
                                        const CubeArraysMap & ListCubes)
{
  CubeArraysMapIterator it = ListCubes.find(cubeNumber);
  blitz::Array<std::complex<float>, 2> Z_local((*it).second.N_RWG_src, (*it).second.N_RWG_src);
  Z_local = std::complex<float>(0.0, 0.0);
  for (int i=0; i<(*it).second.N_RWG_src; i++) Z_local(i, i) = std::complex<float>(1.0, 0.0);

  std::map<int, int> src_edges_numbers_local_src_edges_numbers;
  for (int index = 0; index<(*it).second.N_RWG_src; index++) {
    const int RWG_number = (*it).second.testSrc_RWGsNumbers[index];
    src_edges_numbers_local_src_edges_numbers[RWG_number] = index;
  }
  std::set<int> set_cubeNeighborsIndexes((*it).second.neighborsIndexes.begin(), (*it).second.neighborsIndexes.end());

  for (int index=0; index<(*it).second.N_neighbors; index++) {
    const int neighborCubeNumber = (*it).second.neighborsIndexes[index];
    CubeArraysMapIterator it_neighbor = ListCubes.find(neighborCubeNumber);
    // we first fill in the first lines of Z_local
    if (index==0) {
      for (int i=0; i<(*it).second.N_RWG_test; i++) {
        for (int j=0; j<(*it).second.N_RWG_src; j++) Z_local(i, j) = (*it).second.Z_CFIE_J(i, j);
      }
    }
    // we then fill in the remaining lines
    else {
      // first we find the line indexes
      std::vector<int> Z_local_lines_indexes;
      Z_local_lines_indexes.resize((*it_neighbor).second.N_RWG_test);
      for (int i=0; i<(*it_neighbor).second.N_RWG_test; i++) {
        const int RWG_number = (*it_neighbor).second.testSrc_RWGsNumbers[i];  
        Z_local_lines_indexes[i] = src_edges_numbers_local_src_edges_numbers[RWG_number];
      }
      // then we find the column indexes. A little more complicated
      std::vector<int> common_neighborsNumbers;
      for (int i=0; i<(*it_neighbor).second.N_neighbors; i++) {
        const int val = (*it_neighbor).second.neighborsIndexes[i];
        if (set_cubeNeighborsIndexes.find(val) != set_cubeNeighborsIndexes.end()) common_neighborsNumbers.push_back(val);
      }
      std::vector<int> src_tmp;
      for (int i=0; i<common_neighborsNumbers.size(); i++) {
        const int commonNeighbor = common_neighborsNumbers[i];
        CubeArraysMapIterator it_commonNeighbor = ListCubes.find(commonNeighbor);
        for (int kk=0; kk<(*it_commonNeighbor).second.N_RWG_test; kk++) {
          src_tmp.push_back((*it_commonNeighbor).second.testSrc_RWGsNumbers[kk]);
        }
      }
      std::set<int> set_srcEdges(src_tmp.begin(), src_tmp.end());
      std::vector<int> columnsOfNeighborCubeToBeConsidered;
      for (int i=0; i<(*it_neighbor).second.N_RWG_src; i++) {
        if (set_srcEdges.find((*it_neighbor).second.testSrc_RWGsNumbers[i]) != set_srcEdges.end()) {
          columnsOfNeighborCubeToBeConsidered.push_back(i);
        }
      }
      std::vector<int> Z_local_columns_indexes;
      Z_local_columns_indexes.resize(src_tmp.size());
      for (int i=0; i<src_tmp.size(); i++) {
        const int RWG_number = src_tmp[i];
        Z_local_columns_indexes[i] = src_edges_numbers_local_src_edges_numbers[RWG_number];
      }
      // we construct Z_local
      for (int i=0; i<Z_local_lines_indexes.size(); i++) {
        for (int j=0; j<Z_local_columns_indexes.size(); j++) {
          const int index_line = Z_local_lines_indexes[i];
          const int index_column = Z_local_columns_indexes[j];
          Z_local(index_line, index_column) = (*it_neighbor).second.Z_CFIE_J(i, columnsOfNeighborCubeToBeConsidered[j]);
        }
      }
    } // end else
  } // end for 

  // further reduce the matrix size
  int N_lines = 0;
  for (int i=0; i<(*it).second.N_RWG_src; i++) N_lines += (*it).second.isEdgeInCartesianRadius[i];
  blitz::Array<int, 1> src_edges_numbers_2(N_lines);
  blitz::Array<std::complex<double>, 1> Z_local_2(N_lines * (*it).second.N_RWG_src);
  Z_local_2 = 0.0;
  int line_index = 0;
  for (int i=0; i<(*it).second.N_RWG_src; i++) {
    if ((*it).second.isEdgeInCartesianRadius[i]==1) {
      for (int j=0; j<(*it).second.N_RWG_src; j++) Z_local_2(line_index + j*N_lines) = Z_local (i, j);
      src_edges_numbers_2(line_index) = (*it).second.testSrc_RWGsNumbers[i];
      line_index++;
    }
  }
  // compute the SAI matrix
  blitz::Array<std::complex<double>, 2> Y_CFIE;
  computeMyPinvCC(Y_CFIE, Z_local_2, N_lines, (*it).second.N_RWG_src);
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
    CubeArraysMap ListCubes;
    // variables needed later
    int N_RWG = 0, N_precond = 0, N_q_array = 0;
    std::vector<int> N_ColumnsPerCube;
    N_ColumnsPerCube.resize(N_cubes);
    // we construct a list of cubes
    for (int j=0; j<N_cubes; j++) {
      const int cubeNumber = cubesNumbers(j);
      // we read the cube data and construct the cube
      const string pathToCube = Z_TMP_DATA_PATH + "chunk" + intToString(chunk) + "/";
      const CubeArrays cube(cubeNumber, pathToCube);
      // we add the cube only if it is not in the list
      if (ListCubes.find(cubeNumber) == ListCubes.end()) ListCubes.insert(CubeArraysMap::value_type(cubeNumber, cube));
      // we then add the neighbor cubes to the list
      for (int kk=0; kk<cube.N_neighbors; kk++) {
        const int neighborCubeNumber = cube.neighborsIndexes[kk];
        if (ListCubes.find(neighborCubeNumber) == ListCubes.end()) {
          const int chunkNumberNeighbor = cubeNumber_to_chunkNumber(neighborCubeNumber);
          const string pathToCubeNeighbor = Z_TMP_DATA_PATH + "chunk" + intToString(chunkNumberNeighbor) + "/";
          const CubeArrays cubeNeighbor(neighborCubeNumber, pathToCubeNeighbor);
          ListCubes.insert(CubeArraysMap::value_type(neighborCubeNumber, cubeNeighbor));
        }
      }
      // then we need to construct a series of indexes and offsets
      N_RWG += cube.N_RWG_test;
      int N_isEdgeInCartesianRadius = 0;
      for (int kk=0; kk<cube.isEdgeInCartesianRadius.size(); kk++) N_isEdgeInCartesianRadius += cube.isEdgeInCartesianRadius[kk];
      N_ColumnsPerCube[j] = N_isEdgeInCartesianRadius;
      N_precond += N_isEdgeInCartesianRadius * cube.N_RWG_test;
      N_q_array += N_isEdgeInCartesianRadius;      
    }
    // for the q_array, each src function for all the testing functions of a cube appears only once
    // instead of once per testing function. This allows a dramatic reduction in q_array.size
    std::vector<int> test_RWG_numbers, q_array;
    test_RWG_numbers.resize(N_RWG);
    q_array.resize(N_q_array);
    // N_precond = number of elements in the preconditioner chunk
    blitz::Array<std::complex<float>, 1> Mg(N_precond);
    blitz::Array<int, 2> rowIndexToColumnIndexes(N_RWG, 2);
    int startIndex = 0, startIndexInRWGNumbers = 0, startIndexInQArray = 0;
    int indexN_ColumnsPerCube = 0, index_in_rowIndexToColumnIndexes = 0;
    // constructing test_RWG_numbers
    for (int j=0; j<N_cubes; j++) {
      const int cubeNumber = cubesNumbers(j);
      CubeArraysMapIterator it = ListCubes.find(cubeNumber);
      const CubeArrays cube((*it).second);
      for (int kk=0; kk<cube.N_RWG_test; kk++) test_RWG_numbers[kk + startIndexInRWGNumbers] = cube.testSrc_RWGsNumbers[kk];
      startIndexInRWGNumbers += cube.N_RWG_test;
    }
    // calculation of the SAI preconditioner
    for (int j=0; j<N_cubes; j++) {
      const int cubeNumber = cubesNumbers(j);
      MgPreconditionerComputationPerCube(cubeNumber, ListCubes);
    }
  }
  
  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;
}
