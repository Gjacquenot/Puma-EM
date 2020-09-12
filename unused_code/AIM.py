import os.path
from scipy import *
from weave import converters
from read_mesh_seb import *
from mesh_functions_seb import *
from binomial import *
from EM_constants import *
from ITtest_ITsrc import *

def AIM_Q_numeric(vertexes_coord, triangles_vertexes, triangle_index, R, M, N_points):
    """The wrapping function for 'AIM_Q_numeric.cpp'."""

    Phi = zeros((M, M, M), Float)
    Phi_r = zeros((M, M, M, 3), Float)
    Phi_n_hat_X_r = zeros((M, M, M, 3), Float)

    wrapping_code = """
    using namespace blitz;
    TinyVector<double, 3> R_TinyVector (R(0), R(1), R(2));
    //cout << R_TinyVector << endl;
    AIM_Q_numeric(Phi, Phi_r, Phi_n_hat_X_r, R_TinyVector, vertexes_coord, triangles_vertexes, triangle_index, M, N_points);
    //cout << Phi_r(2,1,1,Range::all()) << endl;
    """

    weave.inline(wrapping_code,
                 ['Phi', 'Phi_r', 'Phi_n_hat_X_r', 'R', 'vertexes_coord', 'triangles_vertexes', 'triangle_index', 'M', 'N_points'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>', '"AIM.h"'],
                 compiler = 'gcc')

    return Phi, Phi_r, Phi_n_hat_X_r

def AIM_S(x):
    """this function wraps the C++ function which computes the inverse of the Vandermonde matrix:

    			[x.^0 ; x.^1 ; ... ; x.^M]

    where x is a horizontal vector of size (M+1) of positions along a line
    """

    S_x = zeros((x.shape[0], x.shape[0]), Float)
    wrapping_code = """AIM_S(S_x, x);"""

    weave.inline(wrapping_code,
                 ['S_x', 'x'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>', '"AIM.h"'],
                 compiler = 'gcc')

    return S_x

def AIM_Lambda(S_x, S_y, S_z, Phi, Phi_r):

    M = S_x.shape[0] - 1
    Lambda = zeros((M+1, M+1, M+1), Float)
    Lambda_r = zeros((M+1, M+1, M+1, 3), Float)
    for i1 in range(M+1):
	    for i2 in range(M+1):
		    for i3 in range(M+1):
			    for j1 in range(M+1):
				    for j2 in range(M+1):
					    for j3 in range(M+1):
						    Lambda[i1, i2, i3] += S_x[i1, j1] * S_y[i2, j2] * S_z[i3, j3] * Phi[j1, j2, j3]
						    Lambda_r[i1, i2, i3, :] += S_x[i1, j1] * S_y[i2, j2] * S_z[i3, j3] * Phi_r[j1, j2, j3, :]

    return Lambda, Lambda_r

if __name__=="__main__":

    M = 4
    N_points = 13
    f = 1.e9
    l = c/f
    d = l*1.5
    eps_r = 1.
    mu_r = 1.
    w = 2*pi*f
    k = w * sqrt(eps_0*eps_r*mu_0*mu_r)

    vertexes_coord = zeros((6, 3), Float)
    vertexes_coord[:3, :] = array([[-0.02, 0, 0], [0, 0, 0], [0, 0.03, 0]],Float)
    vertexes_coord[3:, :] = array([[d, 0, 0], [d+0.02, 0, 0], [d+0.02, 0.01, 0.012]],Float)
    R1 = sum(vertexes_coord[:3])/3.0
    R2 = sum(vertexes_coord[3:])/3.0
    triangles_vertexes = reshape(arange(6),(2,-1))
    index_test_triangle = 0
    index_src_triangle = 1
    Phi_1, Phi_r_1, Phi_n_hat_X_r_1 = AIM_Q_numeric(vertexes_coord, triangles_vertexes, index_test_triangle, R1, M+1, N_points)

    Phi_2, Phi_r_2, Phi_n_hat_X_r_2 = AIM_Q_numeric(vertexes_coord, triangles_vertexes, index_src_triangle, R2, M+1, N_points)

#     Xi_1x = R1[0] + arange(M+1)*l/10.-(M)/2.*l/10
#     Xi_1y = R1[1] + arange(M+1)*l/10.-(M)/2.*l/10
#     Xi_1z = R1[2] + arange(M+1)*l/10.-(M)/2.*l/10
#     Xi_2x = R2[0] + arange(M+1)*l/10.-(M)/2.*l/10
#     Xi_2y = R2[1] + arange(M+1)*l/10.-(M)/2.*l/10
#     Xi_2z = R2[2] + arange(M+1)*l/10.-(M)/2.*l/10
    Xi_1x = arange(M+1)*l/10.-(M)/2.*l/10
    Xi_1y = arange(M+1)*l/10.-(M)/2.*l/10
    Xi_1z = arange(M+1)*l/10.-(M)/2.*l/10
    Xi_2x = arange(M+1)*l/10.-(M)/2.*l/10
    Xi_2y = arange(M+1)*l/10.-(M)/2.*l/10
    Xi_2z = arange(M+1)*l/10.-(M)/2.*l/10

    S_1x = AIM_S(Xi_1x)
    S_1y = AIM_S(Xi_1y)
    S_1z = AIM_S(Xi_1z)
    S_2x = AIM_S(Xi_2x)
    S_2y = AIM_S(Xi_2y)
    S_2z = AIM_S(Xi_2z)
#     W_x = zeros((x.shape[0], x.shape[0]), Float)
#     for i in range(M+1):
# 	    W_x[i, :] = x**i
#     S_x = AIM_S(x)
#     print(W_x )
#     print(S_x )
#     print(matrixmultiply(W_x,S_x))
    Lambda_1, Lambda_r_1 = AIM_Lambda(S_1x, S_1y, S_1z, Phi_1, Phi_r_1)
    Lambda_2, Lambda_r_2 = AIM_Lambda(S_2x, S_2y, S_2z, Phi_2, Phi_r_2)
    u = zeros(3, Float)
    v = zeros(3, Float)
    ITtest_ITsrc_G_AIM = 0
    ITtest_r_ITsrc_G_AIM = zeros(3, Complex)
    ITtest_ITsrc_G_rprime_AIM = zeros(3, Complex)
    ITtest_r_dot_ITsrc_G_rprime_AIM = 0
    for i1 in range(M+1):
	for i2 in range(M+1):
            for i3 in range(M+1):
	        u[0], u[1], u[2] = Xi_1x[i1], Xi_1y[i2], Xi_1z[i3]
		u += R1
	        for j1 in range(M+1):
		    for j2 in range(M+1):
			for j3 in range(M+1):
			    v[0], v[1], v[2] = Xi_2x[j1], Xi_2y[j2], Xi_2z[j3]
			    v += R2
			    R = sqrt(sum((u-v)**2))
			    ITtest_ITsrc_G_AIM += Lambda_1[i1, i2, i3] * exp(-1.j*k*R)/(4*pi*R) * Lambda_2[j1, j2, j3]
			    ITtest_r_ITsrc_G_AIM += Lambda_r_1[i1, i2, i3, :] * exp(-1.j*k*R)/(4*pi*R) * Lambda_2[j1, j2, j3]
			    ITtest_ITsrc_G_rprime_AIM += Lambda_1[i1, i2, i3] * exp(-1.j*k*R)/(4*pi*R) * Lambda_r_2[j1, j2, j3, :]
			    ITtest_r_dot_ITsrc_G_rprime_AIM += sum(Lambda_r_1[i1, i2, i3, :] * exp(-1.j*k*R)/(4*pi*R) * Lambda_r_2[j1, j2, j3, :])
    # now classical integration
    EXTRACT_1_R = EXTRACT_R = 1
    N_points_o = N_points_s = 13
    ITtest_ITsrc_G, ITtest_r_ITsrc_G, ITtest_ITsrc_G_rprime, ITtest_r_dot_ITsrc_G_rprime, ITtest_n_hat_X_r_ITsrc_G, ITtest_n_hat_X_r_dot_ITsrc_G_rprime, ITtest_ITsrc_grad_G, ITtest_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_dot_r_X_ITsrc_grad_G, ITtest_n_hat_X_r_X_ITsrc_grad_G = ITtest_ITsrc_FS(vertexes_coord, triangles_vertexes, index_test_triangle, index_src_triangle, w, eps_r, mu_r, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R)

    print(ITtest_ITsrc_G_AIM)
    print(ITtest_ITsrc_G)
    print(ITtest_r_ITsrc_G_AIM)
    print(ITtest_r_ITsrc_G)
    print(ITtest_ITsrc_G_rprime_AIM)
    print(ITtest_ITsrc_G_rprime)
    print(ITtest_r_dot_ITsrc_G_rprime_AIM)
    print(ITtest_r_dot_ITsrc_G_rprime)