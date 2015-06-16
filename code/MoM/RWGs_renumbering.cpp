#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <blitz/array.h>
#include <map>
#include <mpi.h>

using namespace std;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"

void RWG_renumber(const string pathToReadFrom, const string RWG_type)
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  blitz::Array<int, 1> chunkNumbers;
  string filename = pathToReadFrom + "chunkNumbers.txt";
  readIntBlitzArray1DFromASCIIFile(filename, chunkNumbers);
  
  multimap <int, int> RWG_number_To_Index;
  multimap <int, int>::iterator itr;
  int local_index = 0;
  // first pass to fill in the multimap
  for (unsigned int i=0 ; i<chunkNumbers.size() ; i++) {
    int N_RWG;
    const int chunkNumber = chunkNumbers(i);
    const string chunkNumberString = intToString(chunkNumber);
    // reading the number of RWGs
    filename = pathToReadFrom + "N_" + RWG_type + "_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_RWG);
    
    // reading RWG_numbers
    blitz::Array<int, 1> RWG_numbers(N_RWG);
    filename = pathToReadFrom + RWG_type + "_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWG_numbers);

    // creating the multimap
    for (int j=0; j<N_RWG; j++) {
      const int RWG = RWG_numbers(j);
      if ( (itr=RWG_number_To_Index.find(RWG)) == RWG_number_To_Index.end() ) {
        RWG_number_To_Index.insert(pair<int, int> (RWG, local_index));
        local_index++;
      }
    }
  }
  // we now re-index the RWGs of the multimap by order of magnitude of their key. 
  // Simply traverse the multimap and assign an increasing counter to its second field. 
  // We also create an array of all the locally needed RWGs
  const int N_RWG_local = local_index;
  blitz::Array<int, 1> local_RWG_numbers(N_RWG_local);
  local_index = 0;
  for (itr=RWG_number_To_Index.begin(); itr!=RWG_number_To_Index.end(); ++itr) {
    (*itr).second = local_index;
    local_RWG_numbers(local_index) = (*itr).first;
    local_index++;
  }
  // writing to the file the new total_RWG_numbers to a file
  filename = pathToReadFrom + "local_" + RWG_type + "_RWG_numbers.txt";
  writeIntBlitzArray1DToBinaryFile(filename, local_RWG_numbers);
  filename = pathToReadFrom + "local_N_" + RWG_type + "_RWG.txt";
  writeIntToASCIIFile(filename, N_RWG_local);

  // second pass to re-index the chunk arrays
  for (unsigned int i=0 ; i<chunkNumbers.size() ; i++) {
    int N_RWG;
    const int chunkNumber = chunkNumbers(i);
    const string chunkNumberString = intToString(chunkNumber);
    // reading the number of RWGs
    filename = pathToReadFrom + "N_" + RWG_type + "_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_RWG);

    // reading RWG_numbers
    blitz::Array<int, 1> RWG_numbers(N_RWG);
    filename = pathToReadFrom + RWG_type + "_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWG_numbers);

    // replacing the absolute RWG_numbers by local indexes
    for (int j=0; j<N_RWG; j++) {
      const int RWG = RWG_numbers(j);
      if ( (itr=RWG_number_To_Index.find(RWG)) != RWG_number_To_Index.end() ) RWG_numbers(j) = (*itr).second;
    }
    // writing to the file the new RWG_numbers to a file
    filename = pathToReadFrom + RWG_type + "_RWG_numbers" + chunkNumberString + ".txt";
    writeIntBlitzArray1DToBinaryFile(filename, RWG_numbers);
  }
}

int main(int argc, char* argv[]) {

  MPI::Init();
  const int my_id = MPI::COMM_WORLD.Get_rank(), num_procs = MPI::COMM_WORLD.Get_size();

  string simuDir = ".";
  if ( argc > 2 ) {
     if( string(argv[1]) == "--simudir" ) simuDir = argv[2];
  }

  // general variables
  const string SIMU_DIR = simuDir;
  const string TMP = SIMU_DIR + "/tmp" + intToString(my_id);

  // first Znear. The code for preconditioner will be the same.
  const string pathToReadZnearFrom = TMP + "/Z_near/";
  string RWG_type = "src";
  RWG_renumber(pathToReadZnearFrom, RWG_type);
  RWG_type = "test";
  RWG_renumber(pathToReadZnearFrom, RWG_type);

  // then Mg
  const string pathToReadMgFrom = TMP + "/Mg_LeftFrob/";
  RWG_type = "src";
  RWG_renumber(pathToReadMgFrom, RWG_type);
  RWG_type = "test";
  RWG_renumber(pathToReadMgFrom, RWG_type);

  // Get peak memory usage of each rank
  float memusage_local_MB = static_cast<float>(MemoryUsageGetPeak())/(1024.0*1024.0);
  float tot_memusage_MB, max_memusage_MB;
  MPI_Reduce(&memusage_local_MB, &max_memusage_MB, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&memusage_local_MB, &tot_memusage_MB, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_id==0) {
    std::cout << "MEMINFO " << argv[0] << " : max mem = " << max_memusage_MB << " MB, tot mem = " << tot_memusage_MB << " MB, avg mem = " << tot_memusage_MB/num_procs << " MB" << std::endl;
  }
  MPI::Finalize();
  return 0;  
}
