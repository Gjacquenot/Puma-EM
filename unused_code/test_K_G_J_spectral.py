from scipy import *
from weave import converters

def K_G_J_spectral(k_rho, v, rho, z, zprime, w, eps_0, mu_0, z_i, mu_i, eps_i):
    """small function for wrapping K_G_J_spectral.cpp"""
    N = size(eps_i)
    m = searchsorted(z_i,z)
    n = searchsorted(z_i,zprime)
    k_0 = w * sqrt(eps_0 * mu_0)
    k_i = k_0 * sqrt(eps_i * mu_i)
    d_i = zeros(N,'d')
    d_i[1:-1] = z_i[1:]-z_i[:-1]

    code = """
    using namespace blitz;
    int i;
    complex<double> KA_xx;
    layers_constants LC;
    LC.N = N;
    LC.w = w;
    LC.k_0 = k_0;
    LC.eps_0 = eps_0;
    LC.mu_0 = mu_0;
    LC.z_i.resize(N-1);
    LC.z_i = z_i;
    LC.d_i.resize(N);
    LC.d_i = d_i;
    LC.mu_i.resize(N);
    LC.mu_i = mu_i;
    LC.eps_i.resize(N);
    LC.eps_i = eps_i;
    LC.k_i.resize(N);
    LC.k_i = k_i;
    for (i=0 ; i<100 ; i++) {
        KA_xx = KA_xx_spectral (k_rho, v, rho, z, zprime, m, n, LC);
        cout << KA_xx << endl;
    }
    return_val = KA_xx;
    """

    KA_xx = weave.inline(code,
                         ['k_rho', 'v', 'rho', 'z', 'zprime', 'm', 'n', 'N', 'w', 'k_0', 'eps_0', 'mu_0', 'z_i', 'd_i', 'mu_i', 'eps_i', 'k_i'],
                         type_converters = converters.blitz,
                         include_dirs = ['./MoM/'],
                         library_dirs = ['./MoM/'],
                         libraries = ['MoM'],
                         headers = ['<iostream>','<complex>','"layers_constants.h"','"K_G_J_spectral.h"'],
                         compiler = 'gcc')
    return KA_xx

if __name__=="__main__":
    k_rho = 40. + 10.j
    v = 0.0
    rho = 0.1
    z = -0.07
    zprime = -0.08
    w = 2.0*pi*1e9
    mu_0 = 4.e-7*pi
    c = 2.99792458e8
    eps_0 = 1/(c**2 * mu_0)
    z_i = array([-0.15, -0.1, 0.0],'d')
    eps_i = array([3.-0.3j, 2.-0.3j, 2., 1.],'D')
    mu_i = array([1., 1., 1., 1.],'D')
    KA_xx = K_G_J_spectral(k_rho, v, rho, z, zprime, w, eps_0, mu_0, z_i, mu_i, eps_i)
    print(KA_xx)
