import os.path, sys
from scipy import array, take, argsort, ceil, floor, ones, arange, zeros
from scipy import log2, sqrt
#from scipy import weave
#from scipy.weave import converters

def cube_lower_coord_computation(a, vertexes_coord):
    """This function computes the coordinates of the 1st cube,
    and the number of small cubes in each direction x, y, z.
    a is the length of the side of a cube."""
    x_min = min(vertexes_coord[:, 0])
    x_max = max(vertexes_coord[:, 0])
    Delta_x = x_max - x_min + (sqrt(2.0)-1.0) * a
    y_min = min(vertexes_coord[:, 1])
    y_max = max(vertexes_coord[:, 1])
    Delta_y = y_max - y_min + (sqrt(2.0)-1.0) * a
    z_min = min(vertexes_coord[:, 2])
    z_max = max(vertexes_coord[:, 2])
    Delta_z = z_max - z_min + (sqrt(2.0)-1.0) * a
    big_cube_center_coord = array([x_min + x_max, y_min + y_max, z_min + z_max])/2.0
    N_levels = int(ceil(log2(max(Delta_x, Delta_y, Delta_z)/a)))
    max_N_cubes_1D = 2**N_levels
    length_side_of_big_cube = max(Delta_x, Delta_y, Delta_z)
    big_cube_lower_coord = big_cube_center_coord - length_side_of_big_cube/2.0
    return max_N_cubes_1D, N_levels, big_cube_lower_coord, big_cube_center_coord


def RWGNumber_cubeNumber_computation(a, max_N_cubes_1D, cube_lower_coord, RWGNumber_edgeCentroidCoord):
    """This function finds for each edge the cube to which it belongs.
    a is the length of the side of a cube"""
    RWGNumber_cube = floor((RWGNumber_edgeCentroidCoord - cube_lower_coord)/a).astype('i')
    RWGNumber_cubeNumber = RWGNumber_cube[:, 0] * max_N_cubes_1D**2 + RWGNumber_cube[:, 1] * max_N_cubes_1D + RWGNumber_cube[:, 2]
    RWGNumber_cubeCentroidCoord = cube_lower_coord + a * RWGNumber_cube + ones(3, 'd') * a/2.0
    return RWGNumber_cubeNumber.astype('i'), RWGNumber_cubeCentroidCoord.astype('d')

def cubeIndex_RWGNumbers_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord):
    """each finest-level cube must somehow know which edges it contains.
    This function has the goal of establishing this list for every cube.
    Only the cubes containing edges will be retained.
    We also create a list of the cubes centroids, which will be ordered the same
    way as the cubes_lists_edges_numbers list."""
    E = RWGNumber_cubeNumber.shape[0] # the number of RWGs involved
    ind_sorted_cubes_numbers = argsort(RWGNumber_cubeNumber, kind='mergesort')
    sorted_cubes_numbers = take(RWGNumber_cubeNumber, ind_sorted_cubes_numbers, axis=0)
    sorted_edges_numbers = take(arange(E), ind_sorted_cubes_numbers, axis=0)
    sorted_edges_numbers_cubes_centroids = take(RWGNumber_cubeCentroidCoord, ind_sorted_cubes_numbers, axis=0)
    cubes_lists_edges_numbers = {} # the desired dictionary, renewed for each cube
    cube_list_edges_numbers_tmp = [sorted_edges_numbers[0]] # the temporary list, renewed for each cube
    cubes_centroids = [sorted_edges_numbers_cubes_centroids[0]]
    cubeIndex = 0
    for j in range(E-1): # we cannot go up to (E-1), since (j+1) will then be equal to E (out of bound index)
        if sorted_cubes_numbers[j+1] == sorted_cubes_numbers[j]: # if the next cube number is the same as the current one
            cube_list_edges_numbers_tmp.append(sorted_edges_numbers[j+1]) # add the next element to the temporary list
        else: # if not, we then add the temporary "per-cube" list to the complete list
            cubes_lists_edges_numbers[cubeIndex] = array(cube_list_edges_numbers_tmp)
            cubes_centroids.append(sorted_edges_numbers_cubes_centroids[j+1])
            cube_list_edges_numbers_tmp = [sorted_edges_numbers[j+1]] # init of the temporary list for the next cube
            cubeIndex += 1
    # we must append the last temporary list
    if cubeIndex in cubes_lists_edges_numbers:
        cubes_lists_edges_numbers[cubeIndex+1] = array(cube_list_edges_numbers_tmp)
    else:
        cubes_lists_edges_numbers[cubeIndex] = array(cube_list_edges_numbers_tmp)

    # we transform the "cubes_lists_edges_numbers" in a linear array, useful for the C++ code
    C = len(cubes_lists_edges_numbers)
    cubes_edges_numbers = zeros(E, 'i')
    cube_N_RWGs = zeros(C, 'i')
    startIndex = 0
    for j in range(C):
        length = cubes_lists_edges_numbers[j].shape[0]
        cube_N_RWGs[j] = length
        cubes_edges_numbers[startIndex:startIndex + length] = cubes_lists_edges_numbers[j]
        startIndex += length
    return cubes_edges_numbers, cubes_lists_edges_numbers, cube_N_RWGs.astype('i'), (array(cubes_centroids)).astype('d')

def findCubeNeighbors(max_N_cubes_1D, big_cube_lower_coord, cubes_centroids, a):
    """for each cubes finds its neighbors.
    We use a code similar to Level::searchCubesNeighborsIndexes() from octtree.cpp """
    C = cubes_centroids.shape[0]
    # alternative code
    absoluteCartesianCoord = floor( (cubes_centroids-big_cube_lower_coord)/a )
    CubesNumbers = (absoluteCartesianCoord[:, 0] * max_N_cubes_1D)* max_N_cubes_1D + absoluteCartesianCoord[:, 1] * max_N_cubes_1D + absoluteCartesianCoord[:, 2]
    CubesSortedNumbersToIndexes = zeros((C, 2),'d')
    indSortedCubesNumbers = argsort(CubesNumbers, kind='mergesort')
    CubesSortedNumbersToIndexes[:, 0] = take(CubesNumbers, indSortedCubesNumbers, axis=0)
    CubesSortedNumbersToIndexes[:, 1] = take(arange(C), indSortedCubesNumbers, axis=0)
    cubesNeighborsIndexesTmp2 = zeros((C, 28), 'i') - 1
    wrapping_code2 = """
    int counter;
    for (int i=0 ; i<C ; ++i) {
      blitz::Array<double, 1> absCartCoord(3);
      absCartCoord = absoluteCartesianCoord(i, 0), absoluteCartesianCoord(i, 1), absoluteCartesianCoord(i, 2);
      counter = 1;
      cubesNeighborsIndexesTmp2(i, 0) = i; // we first consider the cube itself
      // we find the neighbors
      for (int x=-1 ; x<2 ; ++x) {
        for (int y=-1 ; y<2 ; ++y) {
          for (int z=-1 ; z<2 ; ++z) {
            int index = -1;
            blitz::Array<double, 1> CandidateAbsCartCoord(3);
            CandidateAbsCartCoord = absCartCoord(0) + x, absCartCoord(1) + y, absCartCoord(2) + z;
            /// no component of (absoluteCartesianCoord(i) + p) -- where i=0,1,2 and p = x,y,z -- can be:
            /// (1) negative or (2) greater than MaxNumberCubes1D.
            int condition = 1;
            for (int j=0 ; j<3 ; ++j) condition *= ( (CandidateAbsCartCoord(j) >= 0) && (CandidateAbsCartCoord(j) < max_N_cubes_1D) );
            // we also do not want to consider the cube itself
            condition *= !((x==0) && (y==0) && (z==0));
            if (condition>0) {
              double candidate_number = (CandidateAbsCartCoord(0) * max_N_cubes_1D)*max_N_cubes_1D + CandidateAbsCartCoord(1) * max_N_cubes_1D + CandidateAbsCartCoord(2);
              { // index search
                if ( (candidate_number < CubesSortedNumbersToIndexes(0, 0)) || (candidate_number > CubesSortedNumbersToIndexes(C-1, 0)) ) index = -1;
                else {
                  int ind_inf = 0, ind_sup = C-1, ind_mid;
                  while(ind_sup-ind_inf > 1) {
                    ind_mid = (ind_sup+ind_inf)/2;
                    if (candidate_number > CubesSortedNumbersToIndexes(ind_mid, 0)) ind_inf = ind_mid;
                    else ind_sup = ind_mid;
                  }
                  if (candidate_number == CubesSortedNumbersToIndexes(ind_inf, 0)) index = CubesSortedNumbersToIndexes(ind_inf, 1);
                  else if (candidate_number == CubesSortedNumbersToIndexes(ind_sup, 0)) index = CubesSortedNumbersToIndexes(ind_sup, 1);
                  else index = -1;
                }
              } // end of index search
            }
            if (index>-1) {cubesNeighborsIndexesTmp2(i, counter) = index; counter++;}
          } // z
        } // y
      } // x
    }
    """
    weave.inline(wrapping_code2,
                 ['C', 'CubesSortedNumbersToIndexes', 'cubesNeighborsIndexesTmp2', 'absoluteCartesianCoord', 'max_N_cubes_1D'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    # construction of "cubesNeighborsIndexes"
    cubes_lists_NeighborsIndexes2 = {}
    N_total_neighbors = 0
    for i in range(C):
        cubes_lists_NeighborsIndexes2[i] = []
    for i in range(C):
        listTmp = []
        j = 0
        while cubesNeighborsIndexesTmp2[i, j] > -1:
            listTmp.append(cubesNeighborsIndexesTmp2[i, j])
            j += 1
        cubes_lists_NeighborsIndexes2[i] = listTmp
        N_total_neighbors += len(listTmp)
    # we also have to save the cubesNeighborsIndexes under a form easily readable by C++ code
    cubes_neighborsIndexes = zeros(N_total_neighbors, 'i')
    cube_N_neighbors = zeros(C, 'i')
    startIndex = 0
    for j in range(C):
        length = len(cubes_lists_NeighborsIndexes2[j])
        cube_N_neighbors[j] = length
        cubes_neighborsIndexes[startIndex:startIndex + length] = cubes_lists_NeighborsIndexes2[j]
        startIndex += length
    return cubes_lists_NeighborsIndexes2, cubes_neighborsIndexes, cube_N_neighbors

def cubes_indexes_to_numbers_computation(a, big_cube_lower_coord, cubes_centroids, N_levels):
    """this function performs the same action as the cube number computation in octtree.cpp
    i.e., it assigns a unique number to each cube, depending upon its position wrt the big cube
    at the coarsest level."""
    cubes_indexes_to_numbers = zeros(cubes_centroids.shape[0], 'i')
    absoluteCartesianCoord = floor( (cubes_centroids-big_cube_lower_coord)/a )
    maxNumberCubes1D = 2.0**(N_levels)
    cubes_indexes_to_numbers = (absoluteCartesianCoord[:,0] * maxNumberCubes1D**2 + absoluteCartesianCoord[:,1] * maxNumberCubes1D + absoluteCartesianCoord[:,2]).astype('i')
    return cubes_indexes_to_numbers

def write_cubes(name, cubes_centroids, a):
    """function that writes the cubes of a given to a file readable by GMSH for viewing."""
    lc = a
    generic_cube = array([[-a/2.0,-a/2.0,-a/2.0],
                          [a/2.0,-a/2.0,-a/2.0],
                          [a/2.0,a/2.0,-a/2.0],
                          [-a/2.0,a/2.0,-a/2.0],
                          [-a/2.0,-a/2.0,a/2.0],
                          [a/2.0,-a/2.0,a/2.0],
                          [a/2.0,a/2.0,a/2.0],
                          [-a/2.0,a/2.0,a/2.0]], 'd')
    lines_generic = array([[1,2], [2,3], [3,4], [4,1], [5,6], [6,7], [7,8], [8,5], [1,5], [2,6], [3,7], [4,8]], 'i')
    f = open(name, 'w')
    f.write('a = ' + str(a) + ';\n')
    f.write('lc = ' + str(lc) + ';\n')
    C = cubes_centroids.shape[0]
    Points = zeros((8,3), 'd')
    for k in range(C):
        for j in range(8):
            Points[j] = cubes_centroids[k] + generic_cube[j]
            PointNumber = str(j+k*8+1)
            Point = str( Points[j, 0] ) + ',' + str( Points[j, 1] ) + ',' + str( Points[j, 2] ) 
            string_to_write = 'Point(' + PointNumber + ') = {' + Point + ', lc};\n'
            f.write(string_to_write)
        for j in range(lines_generic.shape[0]):
            lines = lines_generic + k*8
            LineNumber = str(j+k*12+1)
            string_to_write = 'Line(' + LineNumber + ') = {' + str(lines[j, 0]) + ',' + str(lines[j, 1]) + '};\n'
            f.write(string_to_write)
    f.close()

