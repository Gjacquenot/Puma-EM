import os.path, sys, time, commands
from scipy import zeros, ones, arange, array, take, reshape, sort, argsort, put, sum, compress, nonzero, prod, floor, mean, sqrt, dot, arccos
from scipy import weave
from scipy.weave import converters
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2
from PyGmsh import executeGmsh, write_geo, findParameter, findParameterValue
from EM_constants import *
import copy
from ReadWriteBlitzArray import *


def edges_computation_C(triangle_vertexes, vertexes_coord, saveDir):
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
    V = vertexes_coord.shape[0]
    writeScalarToDisk(T, os.path.join(saveDir, "T.txt"))
    writeScalarToDisk(V, os.path.join(saveDir, "V.txt"))
    writeBlitzArrayToDisk(vertexes_coord, os.path.join(saveDir, 'vertexes_coord') + '.txt')
    writeBlitzArrayToDisk(triangle_vertexes, os.path.join(saveDir, 'triangle_vertexes') + '.txt')

    print commands.getoutput("./code/MoM/mesh_functions_seb " + saveDir + "/")
   
    print "time C++ execution =", time.clock() - t10
    N_RWG = readIntFromDisk(os.path.join(saveDir, "N_RWG.txt"))
    RWGNumber_signedTriangles = readBlitzArrayFromDisk(os.path.join(saveDir, "RWGNumber_signedTriangles.txt"), N_RWG, 2, 'i')
    RWGNumber_edgeVertexes = readBlitzArrayFromDisk(os.path.join(saveDir, "RWGNumber_edgeVertexes.txt"), N_RWG, 2, 'i')
    RWGNumber_oppVertexes = readBlitzArrayFromDisk(os.path.join(saveDir, "RWGNumber_oppVertexes.txt"), N_RWG, 2, 'i')
    is_closed_surface = readASCIIBlitzIntArray1DFromDisk(os.path.join(saveDir, "is_closed_surface.txt"))
    triangles_surfaces = readASCIIBlitzIntArray1DFromDisk(os.path.join(saveDir, "triangles_surfaces.txt"))
    print "    edgeNumber_triangles construction cumulated time =", time.clock() - t10
    sys.stdout.flush()
    return triangles_surfaces, is_closed_surface, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes

def edges_computation_C_old(triangle_vertexes, vertexes_coord, saveDir):
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
    saveDir = saveDir + "/"
    wrapping_code = """
    const int T = triangle_vertexes.extent(0);
    const int E = 3 * T; // there are 3 edges per triangles
    blitz::Array<int, 2> col_sorted_e_v(E, 2);
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
    
    const int max_decimal = max(col_sorted_e_v(blitz::Range::all(), 1));
    double X = 10.0;
    while (X<max_decimal) X *= 10.0; // we look for smallest "X" such that "1eX > max_decimal"

    blitz::Array<double, 1> decimal_e_v(E);
    decimal_e_v = col_sorted_e_v(blitz::Range::all(), 0) + col_sorted_e_v(blitz::Range::all(), 1)/X;

    // we now sort the decimal_e_v
    // we need an argsort type function, given by the Dictionary class (see mesh.h)
    std::vector< Dictionary<double, int> > decimal_e_v_ToIndexes;
    decimal_e_v_ToIndexes.reserve(E);
    for (int j=0 ; j<E ; j++) decimal_e_v_ToIndexes.push_back(Dictionary<double, int> (decimal_e_v(j), j));
    stable_sort(decimal_e_v_ToIndexes.begin(), decimal_e_v_ToIndexes.end());
    blitz::Array<double, 1> sorted_decimal_e_v(E), diff(E);
    blitz::Array<int, 1> ind_sorted_e_v(E);
    for (int j=0 ; j<E ; j++) {
      sorted_decimal_e_v(j) = decimal_e_v(decimal_e_v_ToIndexes[j].getVal());
      ind_sorted_e_v(j) = decimal_e_v_ToIndexes[j].getVal();
    }

    diff = 1.0;
    for (int j=1 ; j<E ; j++) diff(j) = abs(sorted_decimal_e_v(j) - sorted_decimal_e_v(j-1));
    
    blitz::Array<int, 1> indexesEqualPreceding;
    std::vector<int> indexesEqualPrecedingTmp;
    for (int j=0 ; j<E ; j++) {
      if (diff(j)==0.0) indexesEqualPrecedingTmp.push_back(j);
    }
    const int N_indexesEqualPreceding = indexesEqualPrecedingTmp.size();
    indexesEqualPreceding.resize(N_indexesEqualPreceding);
    for (int j=0 ; j<N_indexesEqualPreceding ; j++) indexesEqualPreceding(j) = indexesEqualPrecedingTmp[j];
    indexesEqualPrecedingTmp.clear();

     
    std::string SaveDir = saveDir;
    std::cout << std::endl;
    // compute_indexesEqualEdges
    std::cout << "compute_indexesEqualEdges" << std::endl;
    std::flush(std::cout);
    std::vector<std::vector<int> > indexesEqualEdges;
    compute_indexesEqualEdges(indexesEqualEdges, indexesEqualPreceding, ind_sorted_e_v);
    ind_sorted_e_v.free();
    indexesEqualPreceding.free();

    std::cout << "edgeNumber_vertexes" << std::endl;
    std::flush(std::cout);
    blitz::Array<int,2> edgeNumber_vertexes;
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
    blitz::Array<int, 1> triangles_surfaces(T);
    for (int j=0 ; j<T ; j++) triangles_surfaces(j) = -1;
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
    // computation of opposite vertexes of RWGs in triangles
    blitz::Array<int, 2> RWGNumber_oppVertexes;
    RWGNumber_oppVertexes_computation(RWGNumber_oppVertexes, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes);
    // writing to files
    writeIntToASCIIFile(SaveDir + "N_edges.txt", N_edges);
    writeIntToASCIIFile(SaveDir + "N_RWG.txt", RWGNumber_signedTriangles.extent(0));
    writeIntBlitzArray2DToBinaryFile(SaveDir + "RWGNumber_signedTriangles.txt", RWGNumber_signedTriangles);
    writeIntBlitzArray2DToBinaryFile(SaveDir + "RWGNumber_edgeVertexes.txt", RWGNumber_edgeVertexes);
    writeIntBlitzArray2DToBinaryFile(SaveDir + "RWGNumber_oppVertexes.txt", RWGNumber_oppVertexes);
    writeIntBlitzArray1DToASCIIFile(SaveDir + "is_closed_surface.txt", is_closed_surface);
    writeIntBlitzArray1DToASCIIFile(SaveDir + "triangles_surfaces.txt", triangles_surfaces);
    """
    weave.inline(wrapping_code,
                 ['triangle_vertexes', 'vertexes_coord', 'saveDir'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<string>','<vector>','<algorithm>','"mesh.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    print "time C++ execution =", time.clock() - t10
    N_RWG = readIntFromDisk(saveDir + "N_RWG.txt")
    RWGNumber_signedTriangles = readBlitzArrayFromDisk(saveDir + "RWGNumber_signedTriangles.txt", N_RWG, 2, 'i')
    RWGNumber_edgeVertexes = readBlitzArrayFromDisk(saveDir + "RWGNumber_edgeVertexes.txt", N_RWG, 2, 'i')
    RWGNumber_oppVertexes = readBlitzArrayFromDisk(saveDir + "RWGNumber_oppVertexes.txt", N_RWG, 2, 'i')
    is_closed_surface = readASCIIBlitzIntArray1DFromDisk(saveDir + "is_closed_surface.txt")
    triangles_surfaces = readASCIIBlitzIntArray1DFromDisk(saveDir + "triangles_surfaces.txt")
    print "    edgeNumber_triangles construction cumulated time =", time.clock() - t10
    sys.stdout.flush()
    saveDir = saveDir[:-1]
    return triangles_surfaces, is_closed_surface, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes


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
    targetName = 'strip'
    f = 2.12e9
    write_geo(path, targetName, 'lc', c/f/10.0)
    write_geo(path, targetName, 'lx', 0.07)
    write_geo(path, targetName, 'ly', 0.07)
    write_geo(path, targetName, 'lz', 0.02)
    write_geo(path, targetName, 'w', 0.02)
    executeGmsh(path, targetName, 0)
    targetDimensions_scaling_factor = 1.0
    z_offset = 0.0
    t0 = time.clock()
    #vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    print "reading mesh time =", time.clock() - t0, "seconds"
    sys.stdout.flush()

    triangles_surfaces_C, is_closed_surface_C, RWGNumber_signedTriangles_C, RWGNumber_edgeVertexes_C, RWGNumber_oppVertexes_C = edges_computation_C(triangle_vertexes, vertexes_coord)

    print "    Number of RWG =", RWGNumber_oppVertexes_C.shape[0]


