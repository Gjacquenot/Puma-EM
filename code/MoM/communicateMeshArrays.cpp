
#include <fstream>
#include <iostream>
#include <string>
#include <blitz/array.h>
#include <mpi.h>

using namespace blitz;

#include "readWriteBlitzArrayFromFile.h"
/****************************************************************************/
/******************************* main ***************************************/
/****************************************************************************/

int main(void) {
  int my_id;
  int num_procs;

  MPI::Init();

  num_procs = MPI::COMM_WORLD.Get_size();
//
//  Get the individual process ID.
//
  my_id = MPI::COMM_WORLD.Get_rank();
//  general variables
  if (my_id==0) cout << "Exchanging mesh arrays in C++ (more memory-economic than in Python)......" << endl;
  const string MESH_PATH = "./tmp" + intToString(my_id) + "/mesh/";
  int MAX_N_RWG_per_cube, C, V, T, E, N_levels, S, ierror;
  if (my_id==0) {
    {
      string filename = MESH_PATH + "MAX_N_RWG_per_cube.txt";
      readIntFromASCIIFile(filename, MAX_N_RWG_per_cube);
    }
    {
      string filename = MESH_PATH + "C.txt";
      readIntFromASCIIFile(filename, C);
    }
    {
      string filename = MESH_PATH + "V.txt";
      readIntFromASCIIFile(filename, V);
    }
    {
      string filename = MESH_PATH + "S.txt";
      readIntFromASCIIFile(filename, S);
    }
    {
      string filename = MESH_PATH + "T.txt";
      readIntFromASCIIFile(filename, T);
    }
    {
      string filename = MESH_PATH + "E.txt";
      readIntFromASCIIFile(filename, E);
    }
    {
      string filename = MESH_PATH + "N_levels.txt";
      readIntFromASCIIFile(filename, N_levels);
    }
  }
  // broadcasting the integers
  MPI_Bcast(&MAX_N_RWG_per_cube, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&C, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&V, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&S, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&T, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&E, 1, MPI::INT, 0, MPI::COMM_WORLD);
  MPI_Bcast(&N_levels, 1, MPI::INT, 0, MPI::COMM_WORLD);
  // saving the integers
  if (my_id != 0) {
    string filename = MESH_PATH + "MAX_N_RWG_per_cube.txt";
    writeIntToASCIIFile(filename, MAX_N_RWG_per_cube);
    filename = MESH_PATH + "C.txt";
    writeIntToASCIIFile(filename, C);
    filename = MESH_PATH + "V.txt";
    writeIntToASCIIFile(filename, V);
    filename = MESH_PATH + "S.txt";
    writeIntToASCIIFile(filename, S);
    filename = MESH_PATH + "T.txt";
    writeIntToASCIIFile(filename, T);
    filename = MESH_PATH + "E.txt";
    writeIntToASCIIFile(filename, E);
    filename = MESH_PATH + "N_levels.txt";
    writeIntToASCIIFile(filename, N_levels);
  }
  ierror = MPI_Barrier(MPI::COMM_WORLD);

  // now broadcasting the mesh arrays
  {// broadcasting big_cube_lower_coord
    blitz::Array<double, 1> big_cube_lower_coord(3);
    string filename = MESH_PATH + "big_cube_lower_coord.txt";
    if (my_id==0) readDoubleBlitzArray1DFromBinaryFile(filename, big_cube_lower_coord);
    // broadcasting
    MPI_Bcast(big_cube_lower_coord.data(), big_cube_lower_coord.size(), MPI::DOUBLE, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeDoubleBlitzArray1DToBinaryFile(filename, big_cube_lower_coord);
  }
  {// broadcasting vertexes_coord
    blitz::Array<double, 2> vertexes_coord(V, 3);
    string filename = MESH_PATH + "vertexes_coord.txt";
    if (my_id==0) readDoubleBlitzArray2DFromBinaryFile(filename, vertexes_coord);
    // broadcasting
    MPI_Bcast(vertexes_coord.data(), vertexes_coord.size(), MPI::DOUBLE, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeDoubleBlitzArray2DToBinaryFile(filename, vertexes_coord);
  }
  {// broadcasting triangles_vertexes
    blitz::Array<int, 2> triangle_vertexes(T, 3);
    string filename = MESH_PATH + "triangle_vertexes.txt";
    if (my_id==0) readIntBlitzArray2DFromBinaryFile(filename, triangle_vertexes);
    // broadcasting
    MPI_Bcast(triangle_vertexes.data(), triangle_vertexes.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray2DToBinaryFile(filename, triangle_vertexes);
  }
  {// broadcasting cubes_centroids
    blitz::Array<double, 2> cubes_centroids(C, 3);
    string filename = MESH_PATH + "cubes_centroids.txt";
    if (my_id==0) readDoubleBlitzArray2DFromBinaryFile(filename, cubes_centroids);
    // broadcasting
    MPI_Bcast(cubes_centroids.data(), cubes_centroids.size(), MPI::DOUBLE, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeDoubleBlitzArray2DToBinaryFile(filename, cubes_centroids);
  }
  { // broadcasting cubes_RWGsNumbers
    blitz::Array<int, 2> cubes_RWGsNumbers(C, MAX_N_RWG_per_cube);
    string filename = MESH_PATH + "cubes_RWGsNumbers.txt";
    if (my_id==0) readIntBlitzArray2DFromBinaryFile(filename, cubes_RWGsNumbers);
    // broadcasting
    MPI_Bcast(cubes_RWGsNumbers.data(), cubes_RWGsNumbers.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray2DToBinaryFile(filename, cubes_RWGsNumbers);
  }
  { // broadcasting RWGNumber_signedTriangles
    blitz::Array<int, 2> RWGNumber_signedTriangles(E, 2);
    string filename = MESH_PATH + "RWGNumber_signedTriangles.txt";
    if (my_id==0) readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_signedTriangles);
    // broadcasting
    MPI_Bcast(RWGNumber_signedTriangles.data(), RWGNumber_signedTriangles.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray2DToBinaryFile(filename, RWGNumber_signedTriangles);
  }
  { // broadcasting RWGNumber_edgeVertexes
    blitz::Array<int, 2> RWGNumber_edgeVertexes(E, 2);
    string filename = MESH_PATH + "RWGNumber_edgeVertexes.txt";
    if (my_id==0) readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_edgeVertexes);
    // broadcasting
    MPI_Bcast(RWGNumber_edgeVertexes.data(), RWGNumber_edgeVertexes.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray2DToBinaryFile(filename, RWGNumber_edgeVertexes);
  }
  { // broadcasting RWGNumber_oppVertexes
    blitz::Array<int, 2> RWGNumber_oppVertexes(E, 2);
    string filename = MESH_PATH + "RWGNumber_oppVertexes.txt";
    if (my_id==0) readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_oppVertexes);
    // broadcasting
    MPI_Bcast(RWGNumber_oppVertexes.data(), RWGNumber_oppVertexes.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray2DToBinaryFile(filename, RWGNumber_oppVertexes);
  }
  { // broadcasting triangles_surfaces
    blitz::Array<int, 1> triangles_surfaces(T);
    string filename = MESH_PATH + "triangles_surfaces.txt";
    if (my_id==0) readIntBlitzArray1DFromBinaryFile(filename, triangles_surfaces);
    // broadcasting
    MPI_Bcast(triangles_surfaces.data(), triangles_surfaces.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray1DToBinaryFile(filename, triangles_surfaces);
  }
  { // broadcasting isClosedSurface
    blitz::Array<int, 1> isClosedSurface(S);
    string filename = MESH_PATH + "isClosedSurface.txt";
    if (my_id==0) readIntBlitzArray1DFromBinaryFile(filename, isClosedSurface);
    // broadcasting
    MPI_Bcast(isClosedSurface.data(), isClosedSurface.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray1DToBinaryFile(filename, isClosedSurface);
  }
  { // broadcasting RWGNumber_CFIE_OK
    blitz::Array<int, 1> RWGNumber_CFIE_OK(E);
    string filename = MESH_PATH + "RWGNumber_CFIE_OK.txt";
    if (my_id==0) readIntBlitzArray1DFromBinaryFile(filename, RWGNumber_CFIE_OK);
    // broadcasting
    MPI_Bcast(RWGNumber_CFIE_OK.data(), RWGNumber_CFIE_OK.size(), MPI::INT, 0, MPI::COMM_WORLD);
    // saving
    if (my_id != 0) writeIntBlitzArray1DToBinaryFile(filename, RWGNumber_CFIE_OK);
  }
  ierror = MPI_Barrier(MPI::COMM_WORLD);
  MPI::Finalize();
  return 0;
}

