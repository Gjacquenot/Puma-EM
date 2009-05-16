import os.path
from scipy import *
from scipy import weave
from weave import converters
from EM_constants import *


# def ITtest_ITsrc_FS(test_vertexes_coord, test_triangles_vertexes, src_vertexes_coord, src_triangles_vertexes, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R):
#     """A Python wrapping of the C++ triangle integrations routines.
#     Useful for C++ library testing..."""
# 
#     k = w * sqrt(eps_0*eps_r*mu_0*mu_r) + 0.j
#     # the 4 scalar quantities computed by ITo_ITs_free will be placed in a vector,
#     # which is "scalar_ITtest_ITsrc_G". This is because we cannot change scalar complex
#     # variables defined in Python and passed to C++ without the "return_val" mechanism
#     scalar_ITtest_ITsrc_G = zeros(4, Complex64) 
#     ITtest_r_ITsrc_G = zeros(3, Complex64)
#     ITtest_ITsrc_G_rprime = zeros(3, Complex64)
#     ITtest_n_hat_X_r_ITsrc_G = zeros(3, Complex64)
#     ITtest_ITsrc_grad_G = zeros(3, Complex64)
#     ITtest_r_X_ITsrc_grad_G = zeros(3, Complex64)
#     ITtest_n_hat_X_r_X_ITsrc_grad_G = zeros(3, Complex64)
#     
#     wrapping_code = """
#     using namespace blitz;
#     complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
#     TinyVector<complex<double>, 3> ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_n_hat_X_r_ITs_G;
#     TinyVector<complex<double>, 3> ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G;
#     cout << "blabla" << endl;
#     triangle Ttest, Tsrc;
#     construct_triangle (Ttest, test_vertexes_coord, test_triangles_vertexes, 0);
#     construct_triangle (Tsrc, src_vertexes_coord, src_triangles_vertexes, 0);
# 
#     ITo_ITs_free (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, Ttest, Tsrc, k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);
#     for (int i=0 ; i<3 ; i++) {
#       ITtest_r_ITsrc_G(i) = ITo_r_ITs_G(i);
#       ITtest_ITsrc_G_rprime(i) = ITo_ITs_G_rprime(i);
#       ITtest_n_hat_X_r_ITsrc_G(i) = ITo_n_hat_X_r_ITs_G(i);
#       ITtest_ITsrc_grad_G(i) = ITo_ITs_grad_G(i);
#       ITtest_r_X_ITsrc_grad_G(i) = ITo_r_X_ITs_grad_G(i);
#       ITtest_n_hat_X_r_X_ITsrc_grad_G(i) = ITo_n_hat_X_r_X_ITs_grad_G(i);
#     }
#     scalar_ITtest_ITsrc_G = ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
#     cout << scalar_ITtest_ITsrc_G << endl;
#     """
#     
#     weave.inline(wrapping_code,
#                  ['test_vertexes_coord', 'test_triangles_vertexes', 'src_vertexes_coord', 'src_triangles_vertexes', 'scalar_ITtest_ITsrc_G', 'ITtest_r_ITsrc_G', 'ITtest_ITsrc_G_rprime', 'ITtest_n_hat_X_r_ITsrc_G', 'ITtest_ITsrc_grad_G', 'ITtest_r_X_ITsrc_grad_G', 'ITtest_n_hat_X_r_X_ITsrc_grad_G', 'k', 'N_points_o', 'N_points_s', 'EXTRACT_1_R', 'EXTRACT_R'],
#                  type_converters = converters.blitz,
#                  include_dirs = ['./MoM/'],
#                  library_dirs = ['./MoM/'],
#                  libraries = ['MoM'],
#                  headers = ['<iostream>','"triangle_struct.h"', '"triangle_int.h"','"vector_functions.h"'],
#                  compiler = 'gcc')
#     
#     ITtest_ITsrc_G, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G = scalar_ITtest_ITsrc_G;
#     return ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G


def ITtest_ITsrc_FS(vertexes_coord, triangles_vertexes, index_test_triangle, index_src_triangle, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R):
    """A Python wrapping of the C++ triangle integrations routines.
    Useful for C++ library testing..."""

    k = w * sqrt(eps_0*eps_r*mu_0*mu_r) + 0.j
    # the 4 scalar quantities computed by ITo_ITs_free will be placed in a vector,
    # which is "scalar_ITtest_ITsrc_G". This is because we cannot change scalar complex
    # variables defined in Python and passed to C++ without the "return_val" mechanism
    scalar_ITtest_ITsrc_G = zeros(4, Complex64) 
    ITtest_r_ITsrc_G = zeros(3, Complex64)
    ITtest_ITsrc_G_rprime = zeros(3, Complex64)
    ITtest_n_hat_X_r_ITsrc_G = zeros(3, Complex64)
    ITtest_ITsrc_grad_G = zeros(3, Complex64)
    ITtest_r_X_ITsrc_grad_G = zeros(3, Complex64)
    ITtest_n_hat_X_r_X_ITsrc_grad_G = zeros(3, Complex64)
    
    wrapping_code = """
    using namespace blitz;
    complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
    TinyVector<complex<double>, 3> ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_n_hat_X_r_ITs_G;
    TinyVector<complex<double>, 3> ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G;

    ITo_ITs_free_wrap (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, vertexes_coord, triangles_vertexes, index_test_triangle, index_src_triangle, k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);
    for (int i=0 ; i<3 ; i++) {
      ITtest_r_ITsrc_G(i) = ITo_r_ITs_G(i);
      ITtest_ITsrc_G_rprime(i) = ITo_ITs_G_rprime(i);
      ITtest_n_hat_X_r_ITsrc_G(i) = ITo_n_hat_X_r_ITs_G(i);
      ITtest_ITsrc_grad_G(i) = ITo_ITs_grad_G(i);
      ITtest_r_X_ITsrc_grad_G(i) = ITo_r_X_ITs_grad_G(i);
      ITtest_n_hat_X_r_X_ITsrc_grad_G(i) = ITo_n_hat_X_r_X_ITs_grad_G(i);
    }
    scalar_ITtest_ITsrc_G = ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
    //cout << scalar_ITtest_ITsrc_G << endl;
    """
    
    weave.inline(wrapping_code,
                 ['scalar_ITtest_ITsrc_G', 'ITtest_r_ITsrc_G', 'ITtest_ITsrc_G_rprime', 'ITtest_n_hat_X_r_ITsrc_G', 'ITtest_ITsrc_grad_G', 'ITtest_r_X_ITsrc_grad_G', 'ITtest_n_hat_X_r_X_ITsrc_grad_G', 'vertexes_coord', 'triangles_vertexes', 'index_test_triangle', 'index_src_triangle', 'k', 'N_points_o', 'N_points_s', 'EXTRACT_1_R', 'EXTRACT_R'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>', '"triangle_struct.h"', '"triangle_int.h"'],
                 compiler = 'gcc')
    
    ITtest_ITsrc_G, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G = scalar_ITtest_ITsrc_G;
    return ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G


if __name__=="__main__":
    
    EXTRACT_1_R = EXTRACT_R = 1
    N_points_o = N_points_s = 13
    eps_r = mu_r = 1
    f = 1.e9
    w = 2*pi*f
    
    test_triangles_vertexes = src_triangles_vertexes = reshape(arange(3),(1,-1))
    test_vertexes_coord = array([[-0.02, 0, 0], [0, 0, 0], [0, 0.03, 0]],Float)
    src_vertexes_coord = array([ [0.03, 0, 0], [0.04, 0, 0], [0.04, 0.02, 0]],Float)
    
#     ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G = ITtest_ITsrc_FS(test_vertexes_coord, test_triangles_vertexes, src_vertexes_coord, src_triangles_vertexes, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R)
#     
#     for j in arange(10000):
#         
# 	ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G = ITtest_ITsrc_FS(test_vertexes_coord, test_triangles_vertexes, src_vertexes_coord, src_triangles_vertexes, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R)
#         
# 	print j

    vertexes_coord = zeros((6, 3), Float)
    vertexes_coord[:3, :] = array([[-0.02, 0, 0], [0, 0, 0], [0, 0.03, 0]],Float)
    vertexes_coord[3:, :] = array([ [0.03, 0, 0], [0.04, 0, 0], [0.04, 0.02, 0]],Float)
    triangles_vertexes = reshape(arange(6), (2, -1))
    index_test_triangle = 0
    index_src_triangle = 1

    ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G = ITtest_ITsrc_FS(vertexes_coord, triangles_vertexes, index_test_triangle, index_src_triangle, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R)
    
    print ITtest_ITsrc_G
    print ITtest_r_ITsrc_G
    print ITtest_ITsrc_G_rprime - ITtest_r_ITsrc_G
    print ITtest_r_dot_ITsrc_G_rprime
    print ITtest_n_hat_X_r_dot_ITsrc_G_rprime
    print ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G
    
