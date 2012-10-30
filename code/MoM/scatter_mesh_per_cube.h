#ifndef SCATTER_MESH_PER_CUBE_H
#define SCATTER_MESH_PER_CUBE_H
#include <fstream>
#include <iostream>
#include <string>
#include <blitz/array.h>
#include <map>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"
#include "dictionary.h"


void scatter_mesh_per_cube(blitz::Array<blitz::Array<int, 1>, 1>& allCubeIntArrays,
                           blitz::Array<blitz::Array<double, 1>, 1>& allCubeDoubleArrays,
                           const string simuDir,
                           const blitz::Array<int, 1>& local_ChunksNumbers,
                           const blitz::Array<int, 1>& local_chunkNumber_N_cubesNumbers,
                           const blitz::Array<int, 1>& local_chunkNumber_to_cubesNumbers,
                           const blitz::Array<int, 1>& local_cubeNumber_to_chunkNumbers);

#endif
