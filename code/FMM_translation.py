from scipy import zeros, array, sqrt, dot
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    pass

def IT_theta_IT_phi_alpha_C(r_mn, k, L, XcosTheta, Xphi, ELEM_TYPE):
    """wrapper for C++ IT_theta_IT_phi_alpha_C function"""
    N_points_theta, N_points_phi = XcosTheta.shape[0], Xphi.shape[0]
    alpha_mn = zeros((N_points_theta, N_points_phi), ELEM_TYPE)
    wrapping_code = """IT_theta_IT_phi_alpha_C (alpha_mn, r_mn, k, L, XcosTheta, Xphi);"""
    weave.inline(wrapping_code,
                ['alpha_mn', 'r_mn', 'k', 'L', 'XcosTheta', 'Xphi'],
                type_converters = converters.blitz,
                include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                libraries = ['MoM', 'AMOS', 'g2c'],
                headers = ['<iostream>','<complex>','"FMM.h"'],
                compiler = 'gcc',
                extra_compile_args = ['-O3', '-pthread', '-w'])
    return alpha_mn

#def IT_theta_IT_phi_alpha(r_mn, k, L, XcosTheta, Xphi, ELEM_TYPE):
    #"""wrapper for C++ IT_theta_IT_phi_alpha function"""
    #N_points_theta, N_points_phi = XcosTheta.shape[0], Xphi.shape[0]
    #norm_r_mn = sqrt(dot(r_mn, r_mn))
    #h2_sph = sqrt(pi/(2.* k*norm_r_mn)) * hankel2(arange(L+1)+0.5, k*norm_r_mn)
    #alpha_mn = zeros((N_points_theta, N_points_phi), ELEM_TYPE)
    #wrapping_code = """IT_theta_IT_phi_alpha (alpha_mn, h2_sph, r_mn, k, XcosTheta, Xphi);"""
    #weave.inline(wrapping_code,
                #['alpha_mn', 'h2_sph', 'r_mn', 'k', 'XcosTheta', 'Xphi'],
                #type_converters = converters.blitz,
                #include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                #library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                #libraries = ['MoM', 'AMOS', 'g2c'],
                #headers = ['<iostream>','<complex>','"FMM.h"'],
                #compiler = 'gcc')
    #return alpha_mn

def alpha_computation(cubes_centroids, a, k, L, XcosTheta, Xphi, ELEM_TYPE):
    print("alpha computation")
    N_points_theta, N_points_phi = XcosTheta.shape[0], Xphi.shape[0]
    N_Cx = int(round( ( (max(cubes_centroids[:,0]) + a/2.0) - (min(cubes_centroids[:,0]) - a/2.0) )/a ) )
    N_Cy = int(round( ( (max(cubes_centroids[:,1]) + a/2.0) - (min(cubes_centroids[:,1]) - a/2.0) )/a ) )
    N_Cz = int(round( ( (max(cubes_centroids[:,2]) + a/2.0) - (min(cubes_centroids[:,2]) - a/2.0) )/a ) )
    print("N_Cx =", N_Cx, ", N_Cy =", N_Cy, ", N_Cz =", N_Cz)
    alpha = zeros((2*N_Cx-1, 2*N_Cy-1, 2*N_Cz-1, N_points_theta, N_points_phi), ELEM_TYPE)
    for m in range(2*N_Cx-1):
        for n in range(2*N_Cy-1):
            for p in range(2*N_Cz-1):
                r_mnp_cartesian = array([m-(N_Cx-1), n-(N_Cy-1), p-(N_Cz-1)], 'i')
                r_mnp = r_mnp_cartesian * a
                norm_r_mnp = sqrt(dot(r_mnp, r_mnp))
                IS_NOT_NEIGHBOUR = (norm_r_mnp > sqrt(3)*a + 1.e-5*a)
                if IS_NOT_NEIGHBOUR:
                    alpha[m, n, p] = IT_theta_IT_phi_alpha_C(r_mnp, k, L, XcosTheta, Xphi, ELEM_TYPE)
    return alpha

def alpha_computation_alternative(cubes_centroids, a, k, L, XcosTheta, Xphi, ELEM_TYPE):
    print("alpha computation")
    N_points_theta, N_points_phi = XcosTheta.shape[0], Xphi.shape[0]
    N_Cx = int(round( ( (max(cubes_centroids[:,0]) + a/2.0) - (min(cubes_centroids[:,0]) - a/2.0) )/a ) )
    N_Cy = int(round( ( (max(cubes_centroids[:,1]) + a/2.0) - (min(cubes_centroids[:,1]) - a/2.0) )/a ) )
    N_Cz = int(round( ( (max(cubes_centroids[:,2]) + a/2.0) - (min(cubes_centroids[:,2]) - a/2.0) )/a ) )
    print("N_Cx =", N_Cx, ", N_Cy =", N_Cy, ", N_Cz =", N_Cz)
    alpha = zeros((2*N_Cx-1, 2*N_Cy-1, 2*N_Cz-1, N_points_theta, N_points_phi), ELEM_TYPE)
    for m in range(2*N_Cx-1):
        for n in range(2*N_Cy-1):
            for p in range(2*N_Cz-1):
                r_mnp_cartesian = array([m-(N_Cx-1), n-(N_Cy-1), p-(N_Cz-1)], 'i')
                r_mnp = r_mnp_cartesian * a
                norm_r_mnp = sqrt(dot(r_mnp, r_mnp))
                IS_NOT_NEIGHBOUR = (norm_r_mnp > sqrt(3)*a + 1.e-5*a)
                if IS_NOT_NEIGHBOUR:
                    if p<=N_Cz:
                        alpha[m, n, p] = IT_theta_IT_phi_alpha_C(r_mnp, k, L, XcosTheta, Xphi, ELEM_TYPE)
                    else:
                        alpha[m, n, p] = alpha[m, n, p-N_Cz-1,::-1,:]
    return alpha
