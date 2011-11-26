#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <blitz/array.h>
#include <map>
#include <mpi.h>

using namespace std;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"

void RWG_renumber(const string pathToReadFrom)
{
  blitz::Array<int, 1> chunkNumbers;
  string filename = pathToReadFrom + "chunkNumbers.txt";
  readIntBlitzArray1DFromASCIIFile(filename, chunkNumbers);
  
  multimap <int, int> src_RWG_number_To_Index;
  multimap <int, int>::iterator itr;
  int local_index = 0;
  // first pass to fill in the multimap
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int N_src_RWG;
    const int chunkNumber = chunkNumbers(i);
    const string chunkNumberString = intToString(chunkNumber);
    // reading the number of src RWGs
    filename = pathToReadFrom + "N_src_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_src_RWG);
    
    // reading src_RWG_numbers
    blitz::Array<int, 1> src_RWG_numbers(N_src_RWG);
    filename = pathToReadFrom + "src_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, src_RWG_numbers);

    // creating the multimap
    for (int j=0; j<N_src_RWG; j++) {
      const int src_RWG = src_RWG_numbers(j);
      if ( (itr=src_RWG_number_To_Index.find(src_RWG)) == src_RWG_number_To_Index.end() ) {
        src_RWG_number_To_Index.insert(pair<int, int> (src_RWG, local_index));
        local_index++;
      }
    }
  }

  // we now re-index the src_RWGs of the multimap by order of magnitude of their key. 
  // Simply traverse the multimap and assign an increasing counter to its second field. 
  // We also create an array of all the locally needed src_RWGs
  const int N_src_RWG_local = local_index;
  blitz::Array<int, 1> local_src_RWG_numbers(N_src_RWG_local);
  local_index = 0;
  for (itr=src_RWG_number_To_Index.begin(); itr!=src_RWG_number_To_Index.end(); ++itr) {
    (*itr).second = local_index;
    local_src_RWG_numbers(local_index) = (*itr).first;
    local_index++;
  }
  // writing to the file the new total_src_RWG_numbers to a file
  filename = pathToReadFrom + "local_src_RWG_numbers.txt";
  writeIntBlitzArray1DToBinaryFile(filename, local_src_RWG_numbers);
  filename = pathToReadFrom + "local_N_src_RWG.txt";
  writeIntToASCIIFile(filename, N_src_RWG_local);

  // second pass to re-index the chunk arrays
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int N_src_RWG;
    const int chunkNumber = chunkNumbers(i);
    const string chunkNumberString = intToString(chunkNumber);
    // reading the number of src RWGs
    filename = pathToReadFrom + "N_src_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_src_RWG);

    // reading src_RWG_numbers
    blitz::Array<int, 1> src_RWG_numbers(N_src_RWG);
    filename = pathToReadFrom + "src_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, src_RWG_numbers);

    // replacing the absolute src_RWG_numbers by local indexes
    for (int j=0; j<N_src_RWG; j++) {
      const int src_RWG = src_RWG_numbers(j);
      if ( (itr=src_RWG_number_To_Index.find(src_RWG)) != src_RWG_number_To_Index.end() ) src_RWG_numbers(j) = (*itr).second;
    }
    // writing to the file the new src_RWG_numbers to a file
    filename = pathToReadFrom + "src_RWG_numbers" + chunkNumberString + ".txt";
//    writeIntBlitzArray1DToBinaryFile(filename, src_RWG_numbers);
  }
  
  const int my_id = MPI::COMM_WORLD.Get_rank();
  std::cout << "Process " << my_id << " : local_index = " << local_index << endl;
  flush(std::cout);
}

int main(int argc, char* argv[]) {

  MPI::Init();
  int ierror;
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

  // first Znear. The code for preconditioner will be the same.
  const string pathToReadZnearFrom = TMP + "/Z_near/";
  RWG_renumber(pathToReadZnearFrom);
  
  // then Mg
  const string pathToReadMgFrom = TMP + "/Mg_LeftFrob/";
  RWG_renumber(pathToReadMgFrom);
  
  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;  
}
