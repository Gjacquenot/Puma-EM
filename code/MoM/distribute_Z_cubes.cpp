#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace blitz;

#include "octtree.h"
#include "mesh.h"
#include "EMConstants.h"

int main(void) {

  MPI::Init();
  int ierror, num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();

  // general variables
  const string TMP = "./tmp" + intToString(my_id), OCTTREE_DATA_PATH = TMP + "/octtree_data/", MESH_DATA_PATH = TMP + "/mesh/";
  writeIntToASCIIFile(OCTTREE_DATA_PATH + "CUBES_DISTRIBUTION.txt", 1);
  if (my_id==0) {
      Mesh target_mesh(MESH_DATA_PATH);
      Octtree octtree(OCTTREE_DATA_PATH, target_mesh.cubes_centroids, my_id, num_procs);
    }
  ierror = MPI_Barrier(MPI::COMM_WORLD);
  // we write 0 on the file
  writeIntToASCIIFile(OCTTREE_DATA_PATH + "CUBES_DISTRIBUTION.txt", 0);
  MPI::Finalize();
  return 0;
}
