#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <blitz/array.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace blitz;

#include "GetMemUsage.h"
#include "octtree.h"
#include "EMConstants.h"

int main(int argc, char* argv[]) {

  MPI::Init();
  int num_procs = MPI::COMM_WORLD.Get_size();
  int my_id = MPI::COMM_WORLD.Get_rank();

  string simuDir = ".";
  if ( argc > 2 ) {
     if( string(argv[1]) == "--simudir" ) simuDir = argv[2];
  }

  // general variables
  const string TMP = simuDir + "/tmp" + intToString(my_id);
  const string OCTTREE_DATA_PATH = TMP + "/octtree_data/";
  const string MESH_DATA_PATH = TMP + "/mesh/";
  writeIntToASCIIFile(OCTTREE_DATA_PATH + "CUBES_DISTRIBUTION.txt", 1);
  if (my_id==0) {
    int C;
    string filename = MESH_DATA_PATH + "C.txt";
    readIntFromASCIIFile(filename, C);
    blitz::Array<double, 2> cubes_centroids(C, 3);
    filename = MESH_DATA_PATH + "cubes_centroids.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, cubes_centroids);

    Octtree octtree(OCTTREE_DATA_PATH, cubes_centroids, my_id, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  // we write 0 on the file
  writeIntToASCIIFile(OCTTREE_DATA_PATH + "CUBES_DISTRIBUTION.txt", 0);

  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  MPI::Finalize();
  return 0;
}
