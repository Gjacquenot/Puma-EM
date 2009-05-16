from scipy import *
from weave import converters

def K_G_J(rho, z, zprime, w, eps_0, mu_0, z_i, mu_i, eps_i):
    """small function for wrapping K_G_J.cpp"""
    N = size(eps_i)
    k_0 = w * sqrt(eps_0 * mu_0)
    k_i = k_0 * sqrt(eps_i * mu_i)
    d_i = zeros(N,'d')
    d_i[1:-1] = z_i[1:]-z_i[:-1]

    
    code = """
    using namespace blitz;
    int i;
    complex<double> result;
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
    result = KA_xx (rho, z, zprime, LC);
    cout << result << endl;

    return_val = result;
    """

    result = weave.inline(code,
                         ['rho', 'z', 'zprime', 'N', 'w', 'k_0', 'eps_0', 'mu_0', 'z_i', 'd_i', 'mu_i', 'eps_i', 'k_i'],
                         type_converters = converters.blitz,
                         include_dirs = ['./MoM/'],
                         library_dirs = ['./MoM /'],
                         libraries = ['MoM'],
                         headers = ['<iostream>','<complex>','"layers_constants.h"','"K_G_J.h"'],
                         compiler = 'gcc')
    return result

if __name__=="__main__":
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
    for rho in arange(0,0.2,0.01):
        result = K_G_J(rho, z, zprime, w, eps_0, mu_0, z_i, mu_i, eps_i)
        #print result

