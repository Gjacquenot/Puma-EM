import os.path
from scipy import zeros, ones, arange, reshape, take, put, array, arccos, arcsin, sqrt, dot, sum, real, imag, exp
from scipy import weave, sin, cos, pi
from scipy.weave import converters
from mesh_functions_seb import edgeNumber_triangles_indexes
from meshClass import MeshClass
from EM_constants import *
from PyGmsh import findParameterValue, executeGmsh, write_geo


def G_EJ_G_HJ_cpp(r_dip, r_obs, eps_r, mu_r, k):
    G_EJ = zeros((3,3), 'D')
    G_HJ = zeros((3,3), 'D')
    mu = mu_0*mu_r
    eps = eps_0*eps_r
    wrapping_code = """
    std::vector< std::vector < std::complex<double> > > G_EJ_vector, G_HJ_vector;
    G_EJ_vector.resize(3);
    G_HJ_vector.resize(3);
    for (int i=0; i<3; i++) {
      G_EJ_vector[i].resize(3);
      G_HJ_vector[i].resize(3);
    }
    double rDip[3], rObs[3];
    for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(i);
    for (int i=0 ; i<3 ; ++i) rObs[i] = r_obs(i);
    G_EJ_G_HJ (G_EJ_vector, G_HJ_vector, rDip, rObs, eps, mu, k);
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        G_EJ(i, j) = G_EJ_vector[i][j];
        G_HJ(i, j) = G_HJ_vector[i][j];
      }
    }    
    """
    weave.inline(wrapping_code,
                 ['G_EJ', 'G_HJ', 'r_dip', 'r_obs', 'eps', 'mu', 'k'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>', '<vector>','"V_E_V_H.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return G_EJ, G_HJ

def G_EJ_G_HJ(r_dip, r_obs, eps_r, mu_r, k):
    """inspired from its C++ counterpart"""
    G_EJ = zeros((3,3), 'D')
    G_HJ = zeros((3,3), 'D')
    mu = mu_0*mu_r
    eps = eps_0*eps_r
    r_obs_r_dip = r_obs-r_dip
    R = sqrt(dot(r_obs_r_dip,r_obs_r_dip))
    kRsquare = k*k*R*R
    term_1 = 1.0 + 1.0/(1.j*k*R)
    term_2 = 1.0/R * term_1 + 1.j*k/2.0 * (term_1 - 1.0/kRsquare)
    exp_ikR = exp(-1.j*k*R)
    exp_ikR_R = exp_ikR/R
    x_xp = r_obs_r_dip[0]
    y_yp = r_obs_r_dip[1]
    z_zp = r_obs_r_dip[2]
    x_xp_R_square = (x_xp/R) * (x_xp/R)
    y_yp_R_square = (y_yp/R) * (y_yp/R)
    z_zp_R_square = (z_zp/R) * (z_zp/R)
    G_EJ [0, 1] = term_2*y_yp/R*x_xp/R
    G_EJ [0, 2] = term_2*z_zp/R*x_xp/R 
    G_EJ [1, 2] = term_2*z_zp/R*y_yp/R
    G_EJ [1, 0] = G_EJ [0, 1]
    G_EJ [2, 0] = G_EJ [0, 2]
    G_EJ [2, 1] = G_EJ [1, 2]
    G_EJ [0, 0] = x_xp_R_square*term_1/R-(1.0-x_xp_R_square)*1.j*k/2.0*(term_1 - 1.0/kRsquare)
    G_EJ [1, 1] = y_yp_R_square*term_1/R-(1.0-y_yp_R_square)*1.j*k/2.0*(term_1 - 1.0/kRsquare)
    G_EJ [2, 2] = z_zp_R_square*term_1/R-(1.0-z_zp_R_square)*1.j*k/2.0*(term_1 - 1.0/kRsquare)
    G_EJ = G_EJ*sqrt(mu/eps)/(2.0*pi) * exp_ikR_R

    G_i = exp_ikR/(4.0*pi) * (1.0+1.j*k*R)/(R*R*R)
    G_HJ [0, 1] = (z_zp) * G_i
    G_HJ [1, 0] = -G_HJ [0, 1]
    G_HJ [2, 0] = (y_yp) * G_i
    G_HJ [2, 1] = -1.0*(x_xp) * G_i
    G_HJ [0, 2] = -G_HJ [2, 0]
    G_HJ [1, 2] = -G_HJ [2, 1]

    return G_EJ, G_HJ

def V_EH_dipole(J_dip, r_dip, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes, vertexes_coord, w, eps_r, mu_r):
    """I don't know yet what's gonna go here.
    Anyway, we use prefentially 2-D triangles arrays in the C++ code"""
    # creation of the local V arrays 
    E = list_of_edges_numbers.shape[0]
    V_EH = zeros((E, 4), 'D')
    V_FULL_PRECISION = 1
    # RWGNumber_vertexesCoord
    RWGNumber_vertexesCoord = zeros((E, 6), 'd')
    RWGNumber_vertexesCoord[:, 0:3] = take(vertexes_coord, RWGNumber_edgeVertexes[:,0], axis=0).astype('d')
    RWGNumber_vertexesCoord[:, 3:6] = take(vertexes_coord, RWGNumber_edgeVertexes[:,1], axis=0).astype('d')
    # RWGNumber_oppVertexesCoord
    RWGNumber_oppVertexesCoord = zeros((E, 6), 'd')
    RWGNumber_oppVertexesCoord[:, 0:3] = take(vertexes_coord, RWGNumber_oppVertexes[:,0], axis=0).astype('d')
    RWGNumber_oppVertexesCoord[:, 3:6] = take(vertexes_coord, RWGNumber_oppVertexes[:,1], axis=0).astype('d')
    wrapping_code = """
    std::vector<std::complex<double> > V_tE_J, V_tH_J, V_nE_J, V_nH_J;
    V_tE_J.resize(E);
    V_tH_J.resize(E);
    V_nE_J.resize(E);
    V_nH_J.resize(E);
    double rDip[3];
    std::complex<double> JDip[3];
    for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(i);
    for (int i=0 ; i<3 ; ++i) JDip[i] = J_dip(i);
    V_EJ_HJ_dipole (V_tE_J, V_tH_J, V_nE_J, V_nH_J, JDip, rDip, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_vertexesCoord, RWGNumber_oppVertexesCoord, w, eps_r, mu_r, V_FULL_PRECISION);
    for (int i=0; i<E; i++) {
      V_EH(i, 0) = V_tE_J[i];
      V_EH(i, 1) = V_tH_J[i];
      V_EH(i, 2) = V_nE_J[i];
      V_EH(i, 3) = V_nH_J[i];
    }
    """
    weave.inline(wrapping_code,
                 ['V_EH', 'J_dip', 'r_dip', 'list_of_edges_numbers', 'RWGNumber_CFIE_OK', 'RWGNumber_signedTriangles', 'RWGNumber_vertexesCoord', 'RWGNumber_oppVertexesCoord', 'w', 'eps_r', 'mu_r', 'V_FULL_PRECISION', 'E'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"V_E_V_H.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return V_EH

def V_EH_dipole_alternative(J_dip, r_dip, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes, triangle_vertexes, vertexes_coord, w, eps_r, mu_r):
    pass
    #"""I don't know yet what's gonna go here.
    #Anyway, we use prefentially 2-D triangles arrays in the C++ code"""
    ## creation of the local V arrays 
    #E = list_of_edges_numbers.shape[0]
    #V_EH = zeros((E, 4), 'D')
    #V_FULL_PRECISION = 1
    ## RWGNumber_edgeLength
    #r0_r1 = take(vertexes_coord, RWGNumber_edgeVertexes[:, 0],axis=0) - take(vertexes_coord, RWGNumber_edgeVertexes[:, 1],axis=0)
    #RWGNumber_edgeLength = sqrt(sum(r0_r1 * r0_r1, axis=1)).astype('d')
    ## RWGNumber_oppVertexesCoord
    #RWGNumber_oppVertexesCoord = zeros((E, 6), 'd')
    #RWGNumber_oppVertexesCoord[:, 0:3] = take(vertexes_coord, RWGNumber_oppVertexes[:,0], axis=0).astype('d')
    #RWGNumber_oppVertexesCoord[:, 3:6] = take(vertexes_coord, RWGNumber_oppVertexes[:,1], axis=0).astype('d')
    ## testTriangle_vertexesCoord
    #indexes_test_triangles = edgeNumber_triangles_indexes(list_of_edges_numbers, RWGNumber_signedTriangles).astype('i')
    #testTriangle_vertexes = take(triangle_vertexes, indexes_test_triangles, axis=0)
    #testTriangle_vertexesCoord = zeros((len(indexes_test_triangles), 9), 'd')
    #for i in range(3):
        #testTriangle_vertexesCoord[:, arange(3) + i*3] = take(vertexes_coord, testTriangle_vertexes[:, i], axis=0)
    #wrapping_code = """
    #blitz::Range all = blitz::Range::all();
    #V_EJ_HJ_dipole_alternative (V_EH(all, 0), V_EH(all, 1), V_EH(all, 2), V_EH(all, 3), J_dip, r_dip, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_edgeLength, RWGNumber_oppVertexesCoord, testTriangle_vertexesCoord, w, eps_r, mu_r, V_FULL_PRECISION);
    #"""
    #weave.inline(wrapping_code,
                 #['V_EH', 'J_dip', 'r_dip', 'list_of_edges_numbers', 'RWGNumber_CFIE_OK', 'RWGNumber_signedTriangles', 'RWGNumber_edgeLength', 'RWGNumber_oppVertexesCoord', 'testTriangle_vertexesCoord', 'w', 'eps_r', 'mu_r', 'V_FULL_PRECISION'],
                 #type_converters = converters.blitz,
                 #include_dirs = ['./code/MoM/'],
                 #library_dirs = ['./code/MoM/'],
                 #libraries = ['MoM'],
                 #headers = ['<iostream>','<complex>','"V_E_V_H.h"'],
                 #compiler = 'gcc',
                 #extra_compile_args = ['-O3', '-pthread', '-w'])
    #return V_EH

def V_EH_plane(J_dip, r_dip, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes, vertexes_coord, w, eps_r, mu_r):
    # observation point for the incoming field
    r_ref = zeros(3, 'd') #sum(triangles_centroids, axis=0)/T
    R_hat = (r_dip - r_ref)/sqrt(dot(r_dip - r_ref, r_dip - r_ref))
    k_hat = -R_hat # the propagation vector is indeed opposed to R_hat
    k = w * sqrt(eps_0*eps_r * mu_0*mu_r) # the wavenumber
    G_EJ, G_HJ = G_EJ_G_HJ(r_dip, r_ref, eps_r, mu_r, k)
    E_0 = dot(G_EJ, J_dip).astype('D')
    # creation of the local V arrays
    E = list_of_edges_numbers.shape[0]
    V_EH = zeros((E, 4), 'D')
    V_FULL_PRECISION = 1
    # RWGNumber_vertexesCoord
    RWGNumber_vertexesCoord = zeros((E, 6), 'd')
    RWGNumber_vertexesCoord[:, 0:3] = take(vertexes_coord, RWGNumber_edgeVertexes[:,0], axis=0).astype('d')
    RWGNumber_vertexesCoord[:, 3:6] = take(vertexes_coord, RWGNumber_edgeVertexes[:,1], axis=0).astype('d')
    # RWGNumber_oppVertexesCoord
    RWGNumber_oppVertexesCoord = zeros((E, 6), 'd')
    RWGNumber_oppVertexesCoord[:, 0:3] = take(vertexes_coord, RWGNumber_oppVertexes[:,0], axis=0).astype('d')
    RWGNumber_oppVertexesCoord[:, 3:6] = take(vertexes_coord, RWGNumber_oppVertexes[:,1], axis=0).astype('d')
    wrapping_code = """
    blitz::Range all = blitz::Range::all();
    V_EJ_HJ_plane (V_EH(all, 0), V_EH(all, 1), V_EH(all, 2), V_EH(all, 3), E_0, k_hat, r_ref, list_of_edges_numbers, RWGNumber_CFIE_OK, RWGNumber_signedTriangles, RWGNumber_vertexesCoord, RWGNumber_oppVertexesCoord, w, eps_r, mu_r, V_FULL_PRECISION);
    """
    weave.inline(wrapping_code,
                 ['V_EH', 'E_0', 'k_hat', 'r_ref', 'list_of_edges_numbers', 'RWGNumber_CFIE_OK', 'RWGNumber_signedTriangles', 'RWGNumber_vertexesCoord', 'RWGNumber_oppVertexesCoord', 'w', 'eps_r', 'mu_r', 'V_FULL_PRECISION'],
                 type_converters = converters.blitz,
                 include_dirs = ['./code/MoM/'],
                 library_dirs = ['./code/MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"V_E_V_H.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3', '-pthread', '-w'])
    return V_EH

def computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_edges_numbers, EXCITATION, ELEM_TYPE):
    if EXCITATION=='dipole':
        # V_EH is made of 4 vectors: V_TE_J, V_NE_J, V_TH_J, V_NH_J
        V_EH = V_EH_dipole(J_dip, r_dip, list_of_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r).astype(ELEM_TYPE)
        return V_EH
    elif EXCITATION=='plane':
        # V_EH is made of 4 vectors: V_TE_J, V_NE_J, V_TH_J, V_NH_J
        V_EH = V_EH_plane(J_dip, r_dip, list_of_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r).astype(ELEM_TYPE)
        return V_EH
    elif EXCITATION=='delta_gap':
        print "WARNING!! You asked for delta gap excitation. This is not ready yet. Passing on to plane wave excitation."
        V_EH = V_EH_plane(J_dip, r_dip, list_of_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r).astype(ELEM_TYPE)
        return V_EH
        #if target_mesh.DELTA_GAP:
            #V_EH = V_EH_delta_gap(J_dip, r_dip, list_of_edges_numbers, target_mesh.edges_numbers_triangles, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, target_mesh.triangles_normals, target_mesh.triangles_areas, target_mesh.triangles_centroids, w, eps_r, mu_r).astype(ELEM_TYPE)
            #return V_EH
        #else:
            #print "ERROR!! You asked for delta gap excitation, but none is defined in the file", target_mesh.geoName
            #sys.exit(1)
    else:
        print "ERROR: Wrong excitation setting. Exiting"
        sys.exit(1)
        
        
if __name__=="__main__":
    path = './geo'
    targetName = 'sphere'
    f = 2.12e9
    write_geo(path, targetName, 'lc', c/f/10.1)
    write_geo(path, targetName, 'lx', 0.1)
    write_geo(path, targetName, 'ly', 0.01)
    write_geo(path, targetName, 'lz', 0.0)
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
    MOM_FULL_PRECISION = 1
    list_of_test_edges_numbers = arange(N_RWG).astype('i')
    J_dip = array([1, 0, 0], 'D')
    r_dip = array([0.1, 0.1, 20.0], 'd')
    V_EH = V_EH_dipole(J_dip, r_dip, list_of_test_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r)
    V_EH2 = V_EH_plane(J_dip, r_dip, list_of_test_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r)
    coord = 1
    from pylab import rc, plot, xlabel, ylabel, legend, xticks, yticks, grid, show
    #rc('text', usetex=True)
    FontSize=18
    LineWidth=1
    plot(arange(V_EH[:,coord].shape[0]), real(V_EH[:,coord]), 'b', arange(V_EH[:,coord].shape[0]), real(V_EH2[:,coord]), 'r--', linewidth = LineWidth)
    show()

    r_dip = array([0.1, 0.1, 20.0], 'd')
    r_obs = array([0.15, 0.15, 20.1], 'd')
    mu = mu_0*mu_r
    eps = eps_0*eps_r
    k = w * sqrt(eps_0*eps_r * mu_0*mu_r) # the wavenumber
    G_EJ_cpp, G_HJ_cpp = G_EJ_G_HJ_cpp(r_dip, r_obs, eps_r, mu_r, k)
    G_EJ, G_HJ = G_EJ_G_HJ(r_dip, r_obs, eps_r, mu_r, k)

    print G_EJ-G_EJ_cpp
    print G_HJ_cpp-G_HJ


