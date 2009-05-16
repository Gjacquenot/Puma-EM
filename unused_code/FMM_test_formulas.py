import os.path
from scipy import weave, special, Numeric
from weave import converters
from scipy.special import hankel2, jv
from scipy import integrate
from scipy.integrate import *
from read_mesh_seb import *
from mesh_functions_seb import *
from EM_constants import *

def  P_Legendre(L, z):
    """The wrapping function for 'P_Legendre' in 'FMM.cpp'."""
    
    P = zeros(L+1, Float)
    wrapping_code = """P_Legendre(P, z);"""
    weave.inline(wrapping_code,
                 ['P', 'z'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>', '"FMM.h"'],
                 compiler = 'gcc')
    return P

def cross(a, b):

    c = array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]], Float)
    return c

def sphere_integral_theta_identity(theta, phi, k, d, D_hat, l):
    """The sphere identity is formula (3.13) in [Chew 2001]"""
    
    result = zeros(theta.shape[0], Complex)
    for i in range(theta.shape[0]):
        k_hat = array([sin(theta[i])*cos(phi), sin(theta[i])*sin(phi), cos(theta[i])], Float)
        result[i] = exp(-1.j * k * dot(k_hat, d)) * P_Legendre(l, dot(k_hat, D_hat))[-1] * sin(theta[i])
    return result


def sphere_integral_phi_identity(phi, k, d, D_hat, l):
    
    result = zeros(phi.shape[0], Complex)
    for i in range(phi.shape[0]):
        result[i] = fixed_quad(sphere_integral_theta_identity, 0., pi, args=(phi[i], k, d, D_hat, l), n=25)[0]
    return result
    

def sphere_integral_identity(k, d, D_hat, l):

    result = fixed_quad(sphere_integral_phi_identity, 0., 2*pi, args=(k, d, D_hat, l), n=2*25)[0]
    return result


def XGK_WGK():

    XGK15 = array([0.99145537112081263921, 0.94910791234275852453, 0.86486442335976907279, 0.74153118559939443986, 0.58608723546769113029, 0.40584515137739716691, 0.20778495500789846760, 0.00000000000000000000], Float)
    WGK15 = array([0.02293532201052922496, 0.06309209262997855329, 0.10479001032225018384, 0.14065325971552591875, 0.16900472663926790283, 0.19035057806478540991, 0.20443294007529889241, 0.20948214108472782801], Float)
    return XGK15, WGK15


def sphere_integral_theta_identity_GK(theta_1, theta_2, phi, k, d, D_hat, l):

    XGK15, WGK15 = XGK_WGK()
    D_theta = 0.5 * (theta_2 - theta_1)
    theta_center = 0.5 * (theta_2 + theta_1)
    k_hat = array([sin(theta_center)*cos(phi), sin(theta_center)*sin(phi), cos(theta_center)], Float)
    int_f_theta =  exp(-1.j * k * dot(k_hat, d)) * P_Legendre(l, dot(k_hat, D_hat))[-1] * sin(theta_center) * WGK15[7]
    for p in range(7):
        index = p
        theta = D_theta * XGK15[index]
        theta_minus = theta_center - theta
        theta_plus = theta_center + theta
        k_hat_minus = array([sin(theta_minus)*cos(phi), sin(theta_minus)*sin(phi), cos(theta_minus)], Float)
        k_hat_plus = array([sin(theta_plus)*cos(phi), sin(theta_plus)*sin(phi), cos(theta_plus)], Float)
        f_minus = exp(-1.j * k * dot(k_hat_minus, d)) * P_Legendre(l, dot(k_hat_minus, D_hat))[-1] * sin(theta_minus)
        f_plus = exp(-1.j * k * dot(k_hat_plus, d)) * P_Legendre(l, dot(k_hat_plus, D_hat))[-1] * sin(theta_plus)
        int_f_theta += (f_minus + f_plus) * WGK15[index]
    int_f_theta *= abs(D_theta)
    return int_f_theta

def sphere_integral_phi_identity_GK(phi_1, phi_2, k, d, D_hat, l):

    XGK15, WGK15 = XGK_WGK()
    D_phi = 0.5 * (phi_2 - phi_1)
    phi_center = 0.5 * (phi_2 + phi_1)
    int_g_phi =  sphere_integral_theta_identity_GK(0, pi, phi_center, k, d, D_hat, l) * WGK15[7]
    for q in range(7):
        jtw = q
        phi = D_phi * XGK15[jtw]
        phi_minus = phi_center - phi
        phi_plus = phi_center + phi
        g_minus = sphere_integral_theta_identity_GK(0, pi, phi_minus, k, d, D_hat, l)
        g_plus = sphere_integral_theta_identity_GK(0, pi, phi_plus, k, d, D_hat, l)
        int_g_phi += (g_minus + g_plus) * WGK15[jtw]
    int_g_phi *= abs(D_phi)
    return int_g_phi

## def IT_theta_IT_phi_G_r_ji(r_jm, r_in, k):

##     N_points = 3
##     N_intervals = 5
##     G_r_ji = zeros((N_points * N_intervals, 16), Complex)
##     wrapping_code = """IT_theta_IT_phi_G_r_ji (G_r_ji, r_jm, r_in, k, N_points);"""
##     weave.inline(wrapping_code,
##                  ['G_r_ji', 'r_jm', 'r_in', 'k', 'N_points'],
##                  type_converters = converters.blitz,
##                  include_dirs = ['./MoM/'],
##                  library_dirs = ['./MoM/'],
##                  libraries = ['MoM'],
##                  headers = ['<iostream>','<complex>','"FMM.h"'],
##                  compiler = 'gcc')
##     return G_r_ji


if __name__=="__main__":

    D = array([-0.3, 0.2, 0.6], Float)
    norm_D = sqrt(dot(D, D))
    D_hat = D/norm_D
    theta_D = arccos(D_hat[2])
    print "theta_D =", theta_D*180./pi
    if (abs(D_hat[2])<1.0):
        arccos_phi_D = arccos(D_hat[0]/sin(theta_D))
        print "arccos_phi_D =", arccos(D_hat[0]/sin(theta_D))
        print "arccos(1.0)=", arccos(1.0)
        phi_D = arccos_phi_D
        if (D_hat[1]<0.): 
            phi_D = 2.*pi - arccos_phi_D
    else:
        phi_D = 0.
    print "phi_D =", phi_D*180./pi
    d = array([0.03, 0.01, 0.012], Float)
    norm_d = sqrt(dot(d, d))
    d_hat = d/norm_d
    f = 2.e9
    w = 2.*pi*f
    eps_r = 3.-1.j
    mu_r = 1.
    k = w * sqrt(eps_0*eps_r*mu_0*mu_r)
    L  = int(ceil(abs(k*norm_d + 1. * (k*norm_d)**(1./3.))))
    print "L =", L
    G = exp(-1.j * k * sqrt(dot(D+d, D+d))) / (4 * pi * sqrt(dot(D+d, D+d)))
    
##     for L in range(55):
##  	j_sph = sqrt(pi/(2.* k*norm_d)) * jv(arange(L+1)+0.5, k*norm_d)
##  	h2_sph = sqrt(pi/(2.* k*norm_D)) * hankel2(arange(L+1)+0.5, k*norm_D)
##         P_Leg = P_Legendre(L, dot(d_hat, D_hat))
##         coeff = (-1)**arange(L+1) * (2*arange(L+1) + 1)
##         G_addition_theorem = -1.j*k/(4*pi) * sum(coeff * j_sph * h2_sph * P_Leg)
##         print G, G_addition_theorem

    # here we verify formula (3.13) of [Chew 2001]
    l = L+2
    print sphere_integral_identity(k, d, D_hat, l)
    print sphere_integral_phi_identity_GK(0., pi, k, d, D_hat, l) + sphere_integral_phi_identity_GK(pi, 2.*pi, k, d, D_hat, l)
    print 4*pi*(-1.j)**l * sqrt(pi/(2.* k*norm_d)) * jv(l+0.5, k*norm_d) * P_Legendre(l, dot(d_hat, D_hat))[-1]
