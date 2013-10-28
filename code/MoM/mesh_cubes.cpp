#include <iostream>
#include <string>
#include <blitz/array.h>
#include <vector>
#include <algorithm>
#include <math.h> 

using namespace std;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"

void cube_lower_coord_computation(int & N_levels,
                                  int & max_N_cubes_1D,
                                  std::vector<double> & big_cube_center_coord,
                                  std::vector<double> & big_cube_lower_coord,
                                  const double a, 
                                  const blitz::Array<double, 2>& vertexes_coord, 
                                  const int V)
{
  // This function computes the coordinates of the 1st cube,
  // and the number of small cubes in each direction x, y, z.
  // a is the length of the side of a cube.
  double x_min, x_max, y_min, y_max, z_min, z_max;
  x_min = vertexes_coord(0, 0);
  x_max = vertexes_coord(0, 0);
  y_min = vertexes_coord(0, 1);
  y_max = vertexes_coord(0, 1);
  z_min = vertexes_coord(0, 2);
  z_max = vertexes_coord(0, 2);
  for (int i=1; i<V; i++) {
    x_min = min(vertexes_coord(i, 0), x_min);
    x_max = max(vertexes_coord(i, 0), x_max);
    y_min = min(vertexes_coord(i, 1), y_min);
    y_max = max(vertexes_coord(i, 1), y_max);
    z_min = min(vertexes_coord(i, 2), z_min);
    z_max = max(vertexes_coord(i, 2), z_max);
  }
  const double Delta_x = x_max - x_min + (sqrt(2.0)-1.0) * a;
  const double Delta_y = y_max - y_min + (sqrt(2.0)-1.0) * a;
  const double Delta_z = z_max - z_min + (sqrt(2.0)-1.0) * a;
  const double length_side_of_big_cube = max(Delta_x, max(Delta_y, Delta_z));
  N_levels = max(2, static_cast<int>(ceil(log2(length_side_of_big_cube/a))));
  big_cube_center_coord.resize(3);
  big_cube_center_coord[0] = (x_min + x_max)/2.0;
  big_cube_center_coord[1] = (y_min + y_max)/2.0;
  big_cube_center_coord[2] = (z_min + z_max)/2.0;
  max_N_cubes_1D = pow(2, N_levels);
  big_cube_lower_coord.resize(3);
  for (int i=0; i<3; i++) big_cube_lower_coord[i] = big_cube_center_coord[i] - length_side_of_big_cube/2.0;
}

void compute_RWGNumber_edgeCentroidCoord(blitz::Array<double, 2>& edgeCentroidCoord,
                                         const blitz::Array<double, 2>& vertexes_coord,
                                         const blitz::Array<int, 2>& RWGNumber_edgeVertexes,
                                         const int N_RWG)
{
  for (int i=0; i<N_RWG; i++) {
    const int index_1 = RWGNumber_edgeVertexes(i, 0), index_2 = RWGNumber_edgeVertexes(i, 1);
    for (int j=0; j<3; j++) edgeCentroidCoord(i, j) = (vertexes_coord(index_1, j) + vertexes_coord(index_2, j))/2.0;
  }
}

void RWGNumber_cubeNumber_computation(blitz::Array<int, 1>& RWGNumber_cubeNumber,
                                      blitz::Array<double, 2>& RWGNumber_cubeCentroidCoord,
                                      const double a,
                                      const int max_N_cubes_1D, 
                                      const std::vector<double> & cube_lower_coord, 
                                      const blitz::Array<double, 2>& RWGNumber_edgeCentroidCoord,
                                      const int N_RWG)
{
  // This function finds for each edge the cube to which it belongs.
  // a is the length of the side of a cube
  RWGNumber_cubeNumber.resize(N_RWG);
  RWGNumber_cubeCentroidCoord.resize(N_RWG, 3);
  for (int i=0; i<N_RWG; i++) {
    const int RWGNumber_cube0 = static_cast<int>(floor((RWGNumber_edgeCentroidCoord(i, 0) - cube_lower_coord[0])/a));
    const int RWGNumber_cube1 = static_cast<int>(floor((RWGNumber_edgeCentroidCoord(i, 1) - cube_lower_coord[1])/a));
    const int RWGNumber_cube2 = static_cast<int>(floor((RWGNumber_edgeCentroidCoord(i, 2) - cube_lower_coord[2])/a));    
    RWGNumber_cubeNumber(i) = RWGNumber_cube0 * max_N_cubes_1D*max_N_cubes_1D;
    RWGNumber_cubeNumber(i) += RWGNumber_cube1 * max_N_cubes_1D + RWGNumber_cube2;
    RWGNumber_cubeCentroidCoord(i, 0) = cube_lower_coord[0] + a * RWGNumber_cube0 + a/2.0;
    RWGNumber_cubeCentroidCoord(i, 1) = cube_lower_coord[1] + a * RWGNumber_cube1 + a/2.0;
    RWGNumber_cubeCentroidCoord(i, 2) = cube_lower_coord[2] + a * RWGNumber_cube2 + a/2.0;
  }
}

int main(int argc, char *argv[]) {
  
  int V, T, N_RWG;
  double a;
  const string READING_PATH = argv[1];
  cout << "mesh_cubes.cpp: reading in " << READING_PATH << endl;
  string filename;

  readIntFromASCIIFile(READING_PATH + "T.txt", T);
  readIntFromASCIIFile(READING_PATH + "V.txt", V);
  readIntFromASCIIFile(READING_PATH + "N_RWG.txt", N_RWG);
  readDoubleFromASCIIFile(READING_PATH + "a.txt", a);
    
  blitz::Array<double, 2> vertexes_coord(V, 3);
  readDoubleBlitzArray2DFromBinaryFile(READING_PATH + "vertexes_coord.txt", vertexes_coord);

  int N_levels, max_N_cubes_1D;
  std::vector<double> big_cube_center_coord, big_cube_lower_coord;
  cube_lower_coord_computation(N_levels, max_N_cubes_1D, big_cube_center_coord, big_cube_lower_coord, a, vertexes_coord, V);
  
  cout << "mesh_cubes.cpp: N_levels  = " << N_levels << endl;
  cout << "mesh_cubes.cpp: max_N_cubes_1D = " << max_N_cubes_1D << endl;
  cout << "mesh_cubes.cpp: big_cube_center_coord = " << big_cube_center_coord[0] << ", " << big_cube_center_coord[1] << ", " << big_cube_center_coord[2] << endl;
  cout << "mesh_cubes.cpp: big_cube_lower_coord = " << big_cube_lower_coord[0] << ", " << big_cube_lower_coord[1] << ", " << big_cube_lower_coord[2] << endl;
  
  blitz::Array<int, 2> RWGNumber_edgeVertexes(N_RWG, 2);
  readIntBlitzArray2DFromBinaryFile(READING_PATH + "RWGNumber_edgeVertexes.txt", RWGNumber_edgeVertexes);

  blitz::Array<double, 2> RWGNumber_edgeCentroidCoord(N_RWG, 3);
  compute_RWGNumber_edgeCentroidCoord(RWGNumber_edgeCentroidCoord, vertexes_coord, RWGNumber_edgeVertexes, N_RWG);
  
  blitz::Array<int, 1> RWGNumber_cubeNumber(N_RWG);
  blitz::Array<double, 2> RWGNumber_cubeCentroidCoord(N_RWG, 3);
  RWGNumber_cubeNumber_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord, a, max_N_cubes_1D, big_cube_lower_coord, RWGNumber_edgeCentroidCoord, N_RWG);

  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  return 0;
}
