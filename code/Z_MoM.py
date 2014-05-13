import os.path, time, sys
from scipy import arange, array, zeros, ones, reshape, take, put, sqrt, sum
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    pass
from meshClass import MeshClass, CubeClass
from mesh_functions_seb import edgeNumber_triangles_indexes
from PyGmsh import executeGmsh, write_geo
from EM_constants import *
from ReadWriteBlitzArray import readIntFromDisk, writeScalarToDisk, writeBlitzArrayToDisk, readBlitzArrayFromDisk, read1DBlitzArrayFromDisk

def Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes, vertexes_coord, w, eps_r, mu_r, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION):
    """I don't know yet what's gonna go here.
    Anyway, we use prefentially 2-D triangles arrays in the C++ code"""
    # creation of the local MoM matrices 
    E_test, E_src = list_of_test_edges_numbers.shape[0], list_of_src_edges_numbers.shape[0]
    Z_CFIE_J = zeros((E_test, E_src), 'D')
    Z_CFIE_M = zeros((E_test, E_src), 'D')
    RWGNumber_nodes = zeros((RWGNumber_oppVertexes.shape[0], 4), 'i')
    RWGNumber_nodes[:, 0] = RWGNumber_oppVertexes[:, 0]
    RWGNumber_nodes[:, 1:3] = RWGNumber_edgeVertexes
    RWGNumber_nodes[:, 3] = RWGNumber_oppVertexes[:, 1]
    wrapping_code = """
    Z_CFIE_J_computation(Z_CFIE_J, Z_CFIE_M, CFIE, signSurfObs, signSurfSrc, list_of_test_edges_numbers, list_of_src_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK, RWGNumber_signedTriangles, RWGNumber_nodes, vertexes_coord, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION);
    """
    weave.inline(wrapping_code,
                 ['Z_CFIE_J', 'Z_CFIE_M', 'CFIE', 'signSurfObs', 'signSurfSrc', 'list_of_test_edges_numbers', 'list_of_src_edges_numbers', 'RWGNumber_CFIE_OK', 'RWGNumber_M_CURRENT_OK', 'RWGNumber_signedTriangles', 'RWGNumber_nodes', 'vertexes_coord', 'w', 'eps_r', 'mu_r', 'TDS_APPROX', 'Z_s', 'MOM_FULL_PRECISION'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM','pthread'],
                 headers = ['<iostream>','<complex>','"Z_EJ_Z_HJ.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return Z_CFIE_J, Z_CFIE_M

def Z_MoM_triangles_arraysFromCube(cube, CFIE, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION):
    """I don't know yet what's gonna go here.
    Anyway, we use prefentially 2-D triangles arrays in the C++ code"""
    N_RWG_test, N_RWG_src, N_neighbors, N_nodes, S, testSrc_RWGsNumbers, localTestSrcRWGNumber_signedTriangles, localTestSrcRWGNumber_nodes, localTestRWGNumber_CFIE_OK, localSrcRWGNumber_M_CURRENT_OK, cubeNeighborsIndexes = cube.getIntArrays()
    nodesCoord, rCubeCenter = cube.getDoubleArrays()

    # creation of the local MoM matrices
    Z_CFIE_J = zeros((N_RWG_test, N_RWG_src), 'D')
    #Z_CFIE_M = zeros((N_RWG_test, N_RWG_src), 'D')
    Z_CFIE_M = zeros((1, 1), 'D') # no dielectric in MLFMA yet
    localSrcRWGNumber_M_CURRENT_OK *= 0 # no dielectric in MLFMA yet
    signSurfObs, signSurfSrc = 1.0, 1.0 # no dielectric in MLFMA yet
    wrapping_code = """
    blitz::Array<int, 1> test_RWGsNumbers(N_RWG_test), src_RWGsNumbers(N_RWG_src);
    for (int i=0 ; i<test_RWGsNumbers.size(); i++) test_RWGsNumbers(i) = i;
    for (int i=0 ; i<src_RWGsNumbers.size(); i++) src_RWGsNumbers(i) = i;
    Z_CFIE_J_computation(Z_CFIE_J, Z_CFIE_M, CFIE, signSurfObs, signSurfSrc, test_RWGsNumbers, src_RWGsNumbers, localTestRWGNumber_CFIE_OK, localSrcRWGNumber_M_CURRENT_OK, localTestSrcRWGNumber_signedTriangles, localTestSrcRWGNumber_nodes, nodesCoord, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION);
    """
    weave.inline(wrapping_code,
                 ['Z_CFIE_J', 'Z_CFIE_M', 'CFIE', 'signSurfObs', 'signSurfSrc', 'N_RWG_test', 'N_RWG_src', 'localTestRWGNumber_CFIE_OK', 'localSrcRWGNumber_M_CURRENT_OK', 'localTestSrcRWGNumber_signedTriangles', 'localTestSrcRWGNumber_nodes', 'nodesCoord', 'w', 'eps_r', 'mu_r', 'TDS_APPROX', 'Z_s', 'MOM_FULL_PRECISION'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM','pthread'],
                 headers = ['<iostream>','<complex>','"Z_EJ_Z_HJ.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return Z_CFIE_J

if __name__=="__main__":
    path = './geo'
    targetName = 'sphere'
    f = 2.12e9
    write_geo(path, targetName, 'lc', c/f/9.1)
    write_geo(path, targetName, 'lx', 0.1)
    write_geo(path, targetName, 'ly', 0.1)
    write_geo(path, targetName, 'lz', 0.1)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    languageForMeshConstruction = "Python"
    meshFormat = 'GMSH' 
    meshFileTermination = '.msh'
    target_mesh = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination)
    target_mesh.constructFromGmshFile()
    N_RWG = target_mesh.N_RWG

    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    TDS_APPROX = 0
    Z_s = 0.0
    list_of_test_edges_numbers = arange(N_RWG).astype('i')
    list_of_src_edges_numbers = arange(N_RWG).astype('i')
    CFIE = array([0, 0, 0, 1], 'D')
    signSurfObs, signSurfSrc = 1.0, 1.0
    t0 = time.clock()
    MOM_FULL_PRECISION = 1
    Z_CFIE_J, Z_CFIE_M = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)

    SAVE = False
    if SAVE:
        f = open("Z.txt","w")
        for i in range(Z_CFIE.shape[0]):
            for j in range(Z_CFIE.shape[1]):
                f.write(str(Z_CFIE[i, j]))
                f.write(" ")
            f.write("\n")
        f.close()
    print "time =", time.clock() - t0

