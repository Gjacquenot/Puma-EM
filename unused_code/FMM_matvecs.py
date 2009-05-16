import sys
from scipy import zeros
from scipy import weave
from scipy.weave import converters

def S_EH_J_computation (B_tEJ_src, I_PQ, cubes_edges_numbers, cubes_Ei, ELEM_TYPE):
    C = cubes_Ei.shape[0]
    S_EH_J = zeros((C, B_tEJ_src.shape[1], B_tEJ_src.shape[2], B_tEJ_src.shape[3]), ELEM_TYPE)
    wrapping_code = """S_EH_J_computation (S_EH_J, B_tEJ_src, I_PQ, cubes_edges_numbers, cubes_Ei);"""
    weave.inline(wrapping_code,
                 ['S_EH_J', 'B_tEJ_src', 'I_PQ', 'cubes_edges_numbers', 'cubes_Ei'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 libraries = ['MoM', 'AMOS', 'g2c'],
                 headers = ['<iostream>','<complex>','"FMM.h"'],
                 compiler = 'gcc')
    return S_EH_J

def matvec_Z_PQ_near(Z_PQ_near, I_PQ, pq_array, ELEM_TYPE):
    I_PQ = I_PQ.astype(ELEM_TYPE)
    ZI_PQ_near = zeros(I_PQ.shape[0], ELEM_TYPE)
    wrapping_code = """matvec_Z_PQ_near(ZI_PQ_near, Z_PQ_near, I_PQ, pq_array);"""
    weave.inline(wrapping_code,
                ['ZI_PQ_near', 'Z_PQ_near', 'I_PQ', 'pq_array'],
                type_converters = converters.blitz,
                include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                libraries = ['MoM', 'AMOS', 'g2c'],
                headers = ['<iostream>','<complex>','"FMM.h"'],
                compiler = 'gcc')
    return ZI_PQ_near

def matvec_Z_PQ_far_Ci_Cj(ZI_PQ, B_CFIE_test, S_EH_J_cube_j, edges_numbers_test):
    """Python wrapper for C++ matvec_Z_PQ_far_Ci_Cj function """
    wrapping_code = """matvec_Z_PQ_far_Ci_Cj (ZI_PQ, B_CFIE_test, S_EH_J_cube_j, edges_numbers_test);"""
    weave.inline(wrapping_code,
                 ['ZI_PQ', 'B_CFIE_test', 'S_EH_J_cube_j', 'edges_numbers_test'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 libraries = ['MoM', 'AMOS', 'g2c'],
                 headers = ['<iostream>','<complex>','"FMM.h"'],
                 compiler = 'gcc')


def matvec_Z_PQ_far1(I_PQ, B_CFIE_test, B_tEJ_src, Weights2D, alpha, cubes_lists_edges_numbers, cubes_edges_numbers, cubes_Ei, cubes_centroids, a, ELEM_TYPE):
    N_coord = B_tEJ_src.shape[1]
    ZI_PQ_far = zeros(I_PQ.shape[0], ELEM_TYPE)
    N_Cx = (alpha.shape[0] + 1)/2
    N_Cy = (alpha.shape[1] + 1)/2
    N_Cz = (alpha.shape[2] + 1)/2
    C = cubes_centroids.shape[0]
    S_EH_J = S_EH_J_computation (B_tEJ_src, I_PQ, cubes_edges_numbers, cubes_Ei, ELEM_TYPE)
    sum_S_EH_J = zeros(S_EH_J[0].shape, ELEM_TYPE)
    for i in range(C):
        sum_S_EH_J[:,:,:] = 0.0
        for j in range(C):
            r_ij = cubes_centroids[i]-cubes_centroids[j]
            m, n, p = round(r_ij/a).astype('i') + array([N_Cx-1, N_Cy-1, N_Cz-1], 'i')
            norm_r_ij = sqrt(dot(r_ij, r_ij))
            IS_NOT_NEIGHBOUR = (norm_r_ij > sqrt(3)*a + 1.e-5*a)
            if IS_NOT_NEIGHBOUR:
                for l in range(N_coord):
                    sum_S_EH_J[l, :, :] += S_EH_J[j, l, :, :] * alpha[m, n, p]
        for l in range(N_coord):
            sum_S_EH_J[l, :, :] *= Weights2D
        matvec_Z_PQ_far_Ci_Cj(ZI_PQ_far, B_CFIE_test, sum_S_EH_J, cubes_lists_edges_numbers[i])
    return ZI_PQ_far


def matvec_Z_PQ_far2(I_PQ, B_CFIE_test, B_tEJ_src, Weights2D, alpha, cubes_edges_numbers, cubes_Ei, cubes_centroids, edges_numbers_cubes_centroids, a, L, ELEM_TYPE):
    """Python wrapper for C++ matvec_Z_PQ_far function"""
    I_PQ = I_PQ.astype(ELEM_TYPE)
    ZI_PQ_far = zeros(I_PQ.shape[0], ELEM_TYPE)
    wrapping_code = """matvec_Z_PQ_far(ZI_PQ_far, I_PQ, B_CFIE_test, B_tEJ_src, Weights2D, alpha, cubes_edges_numbers, cubes_Ei, cubes_centroids, edges_numbers_cubes_centroids, a, L);"""
    weave.inline(wrapping_code,
                 ['ZI_PQ_far', 'I_PQ', 'B_CFIE_test', 'B_tEJ_src', 'Weights2D', 'alpha', 'cubes_edges_numbers', 'cubes_Ei', 'cubes_centroids', 'edges_numbers_cubes_centroids', 'a', 'L'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 libraries = ['MoM', 'AMOS', 'g2c'],
                 headers = ['<iostream>','<complex>','"FMM.h"'],
                 compiler = 'gcc',
                 extra_compile_args = ['-O3'])
    return ZI_PQ_far

def matvec_far1(self, x):
    """the far multiplication, python version"""
    return matvec_Z_PQ_far1(x, self.B_CFIE_test, self.B_tEJ_src, self.Weights2D, self.alpha, self.target_mesh.cubes_lists_edges_numbers, self.target_mesh.cubes_edges_numbers, self.target_mesh.cubes_Ei, self.target_mesh.cubes_centroids, self.target_mesh.a, self.ELEM_TYPE)

def matvec_near(self, x):
    """the far multiplication,C++ version"""
    return matvec_Z_PQ_near(self.Z_CFIE_near, x, self.pq_array, self.ELEM_TYPE)

def matvec_FMM1(self, x):
    self.numberOfIterations += 1
    y1 = matvec_Z_PQ_near(self.Z_CFIE_near, x, self.pq_array, self.ELEM_TYPE) + matvec_Z_PQ_far1(x, self.B_CFIE_test, self.B_tEJ_src, self.Weights2D, self.alpha, self.target_mesh.cubes_lists_edges_numbers, self.target_mesh.cubes_edges_numbers, self.target_mesh.cubes_Ei, self.target_mesh.cubes_centroids, self.target_mesh.a, self.ELEM_TYPE)
    sys.stdout.write("\r" + str(self.numberOfIterations) + " " + str(sum(x*self.V_EH[:,0])))
    sys.stdout.flush()
    return y1

def matvec_FMM2(self, x):
    self.numberOfIterations += 1
    y1 = matvec_Z_PQ_near(self.Z_CFIE_near, x, self.pq_array, self.ELEM_TYPE) + matvec_Z_PQ_far2(x, self.B_CFIE_test, self.B_tEJ_src, self.Weights2D, self.alpha, self.target_mesh.cubes_edges_numbers, self.target_mesh.cubes_Ei, self.target_mesh.cubes_centroids, self.target_mesh.edges_numbers_cubes_centroids, self.target_mesh.a, self.L, self.ELEM_TYPE)
    sys.stdout.write("\r" + str(self.numberOfIterations) + " " + str(sum(x*self.V_EH[:,0])))
    sys.stdout.flush()
    return y1
