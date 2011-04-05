import os.path, sys, time
from scipy import zeros, ones, arange, array, take, reshape, sort, argsort, put, sum, compress, nonzero, prod, floor, mean, sqrt, dot, arccos
from scipy import weave
from scipy.weave import converters
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2
from PyGmsh import executeGmsh, write_geo, findParameter, findParameterValue
from EM_constants import *
import copy
from ReadWriteBlitzArray import *

def edges_computation_C(triangle_vertexes, vertexes_coord):
    """This function builds the edges matrix from the triangles"""

    # For assigning a "number" and a "kind" to an edge (a set of two vertexes),
    # we need to find how many occurrences it has in "edges_vertexes". A costful
    # operation, because for each edge we have to compare the edge to all other
    # that appear in the "edges_vertexes" array.
    #
    # However, if we want efficiency, we can create an array of edges sorted as follows:
    # 1) sort the elements of the first two columns alongside dimension 1; 
    # 2) sort the pairs following their 1st element alongside dimension 0;
    # 3) sort the pairs following their 2nd element alongside dimension 0,
    #    but with keeping the order of their first element.
    # In such an array, all occurrences of an edge would be adjacent to each other. So, for
    # counting the occurrences of an edge, one should only count its similar
    # neighbors, thereby greatly reducing the computational cost of the algorithm

    # Once the elements of the first 2 columns have been sorted alongside dimension 1,
    # we construct a 1-D array of real numbers, with:
    # 1) the entire part formed by the numbers of the first column
    # 2) the decimal part formed by the numbers of the second column

    t10 = time.clock()
    print "    construction of edgeNumber_triangles..."
    sys.stdout.flush()
    
    T = triangle_vertexes.shape[0]
    E = T * 3 # there are 3 edges per triangles
    saveDir = "./geo/" # where we will write the temporary files
    t0 = time.clock()
    print "constructing col_sorted_e_v..."
    sys.stdout.flush() 
    col_sorted_e_v = zeros((E, 2), 'i')
    wrapping_code = """
    const int T = triangle_vertexes.extent(0);
    for (int i=0 ; i<T ; i++) {
      for (int j=0 ; j<3 ; j++) {
        int n_orig = j;
        int n_end = (j<2) ? j+1 : 0;
        int r_orig = triangle_vertexes(i, n_orig);  
        int r_end = triangle_vertexes(i, n_end);
        int index = i*3 + j;
        col_sorted_e_v(index, 0) = std::min(r_orig, r_end);
        col_sorted_e_v(index, 1) = std::max(r_orig, r_end);
      } 
    }
    """
    weave.inline(wrapping_code,
                 ['col_sorted_e_v', 'triangle_vertexes'],
                 type_converters = converters.blitz,
                 include_dirs = [],
                 library_dirs = [],
                 libraries = [],
                 headers = ['<iostream>'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    print "time constructing col_sorted_e_v =", time.clock() - t0
    sys.stdout.flush()

    max_decimal = max(col_sorted_e_v[:, 1]) # the maximum element of "col_sorted_e_v[:, 1]"
    X = 10.0
    while X < max_decimal: # we look for smallest "X" such that "1eX > max_decimal"
        X *= 10.0
    decimal_e_v = col_sorted_e_v[:, 0] + col_sorted_e_v[:, 1]/X
    writeScalarToDisk(len(col_sorted_e_v[:,0]), saveDir + "N_col_sorted_e_v.txt")
    writeBlitzArrayToDisk(col_sorted_e_v, saveDir + "col_sorted_e_v.txt")
    del col_sorted_e_v
    # ind_sorted_e_v
    ind_sorted_e_v = argsort(decimal_e_v, kind='mergesort').astype('i')
    sorted_decimal_e_v = take(decimal_e_v, ind_sorted_e_v, axis=0)
    writeScalarToDisk(len(ind_sorted_e_v), saveDir + "N_ind_sorted_e_v.txt")
    writeBlitzArrayToDisk(ind_sorted_e_v, saveDir + "ind_sorted_e_v.txt")
    del decimal_e_v, ind_sorted_e_v

    diff = ones(sorted_decimal_e_v.shape[0], 'd')
    diff[1:] = abs(sorted_decimal_e_v[1:] - sorted_decimal_e_v[:-1])
    del sorted_decimal_e_v
    indexesEqualPreceding = compress(diff==0.0,arange(len(diff)),axis=0).astype('i')
    writeScalarToDisk(len(indexesEqualPreceding), saveDir + "N_indexesEqualPreceding.txt")
    writeBlitzArrayToDisk(indexesEqualPreceding, saveDir + "indexesEqualPreceding.txt")
    del diff, indexesEqualPreceding
    t0 = time.clock()
    print "Entering the C++ area...",
    sys.stdout.flush()
    triangles_surfaces = (ones(T, 'i') * -1).astype('i')
    wrapping_code = """
    const int T = triangle_vertexes.extent(0);
    std::string SaveDir = saveDir;
    std::cout << std::endl;
    // compute_indexesEqualEdges
    std::cout << "compute_indexesEqualEdges" << std::endl;
    std::flush(std::cout);
    std::vector<std::vector<int> > indexesEqualEdges;
    int N_ind_sorted_e_v, N_indexesEqualPreceding;
    readIntFromASCIIFile(SaveDir + "N_ind_sorted_e_v.txt", N_ind_sorted_e_v);
    readIntFromASCIIFile(SaveDir + "N_indexesEqualPreceding.txt", N_indexesEqualPreceding);
    blitz::Array<int, 1> ind_sorted_e_v(N_ind_sorted_e_v), indexesEqualPreceding(N_indexesEqualPreceding);
    readIntBlitzArray1DFromBinaryFile(SaveDir + "ind_sorted_e_v.txt", ind_sorted_e_v);
    readIntBlitzArray1DFromBinaryFile(SaveDir + "indexesEqualPreceding.txt", indexesEqualPreceding);
    compute_indexesEqualEdges(indexesEqualEdges, indexesEqualPreceding, ind_sorted_e_v);
    ind_sorted_e_v.free();
    indexesEqualPreceding.free();

    std::cout << "edgeNumber_vertexes" << std::endl;
    std::flush(std::cout);
    blitz::Array<int,2> edgeNumber_vertexes;
    int N_col_sorted_e_v;
    readIntFromASCIIFile(SaveDir + "N_col_sorted_e_v.txt", N_col_sorted_e_v);
    blitz::Array<int, 2> col_sorted_e_v(N_col_sorted_e_v, 2);
    readIntBlitzArray2DFromBinaryFile(SaveDir + "col_sorted_e_v.txt", col_sorted_e_v);
    compute_edgeNumber_vertexes(edgeNumber_vertexes, indexesEqualEdges, col_sorted_e_v);
    col_sorted_e_v.free();
    const int N_edges = edgeNumber_vertexes.extent(0);

    std::cout << "edgeNumber_triangles" << std::endl;
    std::flush(std::cout);
    blitz::Array<int, 2> edgeNumber_triangles;
    compute_edgeNumber_triangles(edgeNumber_triangles, indexesEqualEdges);
    indexesEqualEdges.clear(); // not needed anymore

    // compute_triangle_adjacentTriangles
    std::cout << "compute_triangle_adjacentTriangles" << std::endl;
    std::flush(std::cout);
    std::vector<std::vector<int> > triangle_adjacentTriangles;
    compute_triangle_adjacentTriangles(triangle_adjacentTriangles, edgeNumber_triangles, T);

    // reordering vertexes of the triangles
    std::cout << "reorder_triangle_vertexes" << std::endl;
    std::flush(std::cout);
    reorder_triangle_vertexes(triangle_vertexes, triangles_surfaces, vertexes_coord, triangle_adjacentTriangles);
    triangle_adjacentTriangles.clear();

    // finding the open and closed surfaces
    std::cout << "is_surface_closed" << std::endl;
    std::flush(std::cout);
    blitz::Array<int, 1> is_closed_surface;
    blitz::Array<std::vector<int>, 2> connected_surfaces, potential_closed_surfaces;
    is_surface_closed(is_closed_surface, connected_surfaces, potential_closed_surfaces, triangles_surfaces, edgeNumber_triangles);

    // RWGNumber_signedTriangles_computation
    std::cout << "RWGNumber_signedTriangles_computation" << std::endl;
    std::flush(std::cout);
    blitz::Array<int, 2> RWGNumber_signedTriangles, RWGNumber_edgeVertexes;
    RWGNumber_signedTriangles_computation(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, edgeNumber_triangles, edgeNumber_vertexes, triangles_surfaces, is_closed_surface, triangle_vertexes, vertexes_coord);
    edgeNumber_triangles.free();
    edgeNumber_vertexes.free();
    writeIntToASCIIFile(SaveDir + "N_edges.txt", N_edges);
    writeIntToASCIIFile(SaveDir + "N_RWG.txt", RWGNumber_signedTriangles.extent(0));
    writeIntBlitzArray2DToBinaryFile(SaveDir + "RWGNumber_signedTriangles.txt", RWGNumber_signedTriangles);
    writeIntBlitzArray2DToBinaryFile(SaveDir + "RWGNumber_edgeVertexes.txt", RWGNumber_edgeVertexes);
    writeIntBlitzArray1DToASCIIFile(SaveDir + "is_closed_surface.txt", is_closed_surface);
    """
    weave.inline(wrapping_code,
                 ['triangle_vertexes', 'triangles_surfaces', 'vertexes_coord', 'saveDir'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<string>', '<complex>','<vector>','<algorithm>','"mesh.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    print "time C++ execution =", time.clock() - t0
    N_RWG = readIntFromDisk(saveDir + "N_RWG.txt")
    RWGNumber_signedTriangles = readBlitzArrayFromDisk(saveDir + "RWGNumber_signedTriangles.txt", N_RWG, 2, 'i')
    RWGNumber_edgeVertexes = readBlitzArrayFromDisk(saveDir + "RWGNumber_edgeVertexes.txt", N_RWG, 2, 'i')
    is_closed_surface = readASCIIBlitzIntArray1DFromDisk(saveDir + "is_closed_surface.txt")
    print "    edgeNumber_triangles construction cumulated time =", time.clock() - t10
    sys.stdout.flush()
    return triangles_surfaces, is_closed_surface, RWGNumber_signedTriangles, RWGNumber_edgeVertexes

def RWGNumber_oppVertexes_computation_C(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes):
    print "    computation of RWG opposite vertexes...", 
    sys.stdout.flush()
    t5 = time.clock()
    N_RWG = RWGNumber_signedTriangles.shape[0]
    RWGNumber_oppVertexes = zeros((N_RWG, 2), 'i')
    wrapping_code = """
    RWGNumber_oppVertexes_computation(RWGNumber_oppVertexes, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes);
    """
    weave.inline(wrapping_code,
                 ['RWGNumber_oppVertexes', 'RWGNumber_signedTriangles', 'RWGNumber_edgeVertexes', 'triangle_vertexes'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','"mesh.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    print "   time =", time.clock() - t5
    return RWGNumber_oppVertexes.astype('i')

def edgeNumber_triangles_indexes_C(N_triangles, list_of_edges_numbers, RWGNumber_signedTriangles):
    """This function returns a 1-D array of the indexes of the triangles corresponding
    to a 1-D array of edges_numbers. This function is important for creating lists of triangles
    that will participate to the MoM, given a particular criterium concerning the edges.

    Same as above but wraps a C++ function... As we can't resize the 'indexes_of_triangles'
    array in C++ when it has been created in Python, we need to pass as an argument the size
    of 'indexes_of_triangles', which is 'N_triangles'.
    """
    indexes_of_triangles = zeros(N_triangles, 'i')
    wrapping_code = """compute_RWGNumber_trianglesNumbers(indexes_of_triangles, list_of_edges_numbers, RWGNumber_signedTriangles);"""
    weave.inline(wrapping_code,
                 ['indexes_of_triangles', 'list_of_edges_numbers', 'RWGNumber_signedTriangles'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','"mesh.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return indexes_of_triangles.astype('i')

if __name__=="__main__":
    path = './geo'
    targetName = 'cubi'
    f = 2.12e9
    write_geo(path, targetName, 'lc', c/f/10.0)
    write_geo(path, targetName, 'lx', 6.05)
    write_geo(path, targetName, 'ly', 6.05)
    write_geo(path, targetName, 'lz', 6.05)
    write_geo(path, targetName, 'w', 0.02)
    executeGmsh(path, targetName, 0)
    targetDimensions_scaling_factor = 1.0
    z_offset = 0.0
    t0 = time.clock()
    #vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    print "reading mesh time =", time.clock() - t0, "seconds"
    sys.stdout.flush()

    triangles_surfaces_C, is_closed_surface_C, RWGNumber_signedTriangles_C, RWGNumber_edgeVertexes_C = edges_computation_C(triangle_vertexes, vertexes_coord)
    RWGNumber_oppVertexes_C = RWGNumber_oppVertexes_computation_C(RWGNumber_signedTriangles_C, RWGNumber_edgeVertexes_C, triangle_vertexes)

    print "    Number of RWG =", RWGNumber_oppVertexes_C.shape[0]



