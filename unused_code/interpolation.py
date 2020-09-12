from scipy import *
from scipy import gplt
from weave import converters
from integration import *
from pylab import *
import time

def findIndex(x, xi, CYCLIC):
    """finds the lower index of x in xi"""
    N = xi.shape[0]
    if ( (x < xi[0]) ):
        if CYCLIC:
            return -1 # index of last element
        else:
            return 0
    elif (x > xi[N-1]):
        return N-1
    else:
        ind_inf = 0
        ind_sup = N-1
        while(ind_sup-ind_inf > 1):
            ind_mid = (ind_sup+ind_inf)/2
            if (x <= xi[ind_mid]):
                ind_sup = ind_mid
            else:
                ind_inf = ind_mid
    return ind_inf

def gaussInterpCoeffs(x, xi):
    Ni = xi.shape[0]
    coeffs = zeros(Ni, Float)
    for i in range(Ni):
        B, C = 1., 1.
        for j in range(Ni):
            if (j!=i):
                B *= sin(0.5*(x - xi[j]))
                C *= sin(0.5*(xi[i] - xi[j]))
        coeffs[i] = B/C
    return coeffs

def gaussInterp(x, a, b, xi, yi, NOrder, CYCLIC, int_method, IS_THETA, yi_OPP_side, INCLUDE_THETA_BOUNDARIES):
    """
    Gauss interpolation for periodic functions
    a and b are the limits of the interpolation interval
    These boundaries may be already included in xi
    yi_OPP_side is used iif (IS_THETA == 1)
    """
    Ni = xi.shape[0]
    k = findIndex(x, xi, CYCLIC)
    startInd = k - NOrder/2
    if not(CYCLIC):
        if ( startInd<0 ):
            startInd = 0
        elif ( startInd+NOrder>Ni-1 ):
            startInd = (Ni-1) - NOrder
    indexesTmp = arange(startInd, startInd+NOrder+1)
    if not(CYCLIC):
        xiTmp = take(xi, indexesTmp)
    if CYCLIC:
        deltaX = xi[1]-xi[0]
        if int_method in ["PONCELET", "TRAP"]:
            xiTmp = indexesTmp * deltaX + xi[0]
        else: # int_method==GAUSSL
            xiTmp = zeros(indexesTmp.shape[0], Float)
            if (INCLUDE_THETA_BOUNDARIES==0):
                for i in range(xiTmp.shape[0]):
                    ind = indexesTmp[i]
                    if ind<0:
                        xiTmp[i] = xi[ind] + a - b
                    elif ind>Ni-1:
                        xiTmp[i] = xi[ind-Ni] + b - a
                    else:
                        xiTmp[i] = xi[ind]
            else: # (INCLUDE_THETA_BOUNDARIES==1)
                for i in range(xiTmp.shape[0]):
                    ind = indexesTmp[i]
                    if ind<0:
                        xiTmp[i] = xi[ind-1] + a - b # we "jump" over xi[-1] (because xi[-1]==b)
                    elif ind>Ni-1:
                        xiTmp[i] = xi[ind-Ni+1] + b - a # we "jump" over xi[0] (because xi[0]==a)
                    else:
                        xiTmp[i] = xi[ind]

    coeffs = gaussInterpCoeffs(x, xiTmp)
    # !! if cyclic interpolation in theta: if theta<0, F(theta, phi) = -F(|theta|, phi + pi)
    yTmp = zeros(indexesTmp.shape[0], Complex)
    if (IS_THETA==0):
        for i in range(yTmp.shape[0]):
            if (indexesTmp[i]<0):
                yTmp[i] = yi[indexesTmp[i]]
            elif (indexesTmp[i]>Ni-1):
                yTmp[i] = yi[indexesTmp[i]-Ni]
            else:
                yTmp[i] = yi[indexesTmp[i]]
    else:
        Ni_opp = yi_OPP_side.shape[0]
        y_complete = zeros(Ni+Ni_opp, Complex)
        y_complete[:Ni_opp] = -yi_OPP_side[-1::-1]
        y_complete[Ni_opp:] = yi
        for i in range(yTmp.shape[0]):
            if (indexesTmp[i]<0):
                yTmp[i] = y_complete[Ni_opp+indexesTmp[i]]
            elif (indexesTmp[i]>Ni-1):
                yTmp[i] = y_complete[indexesTmp[i]-Ni]
            else:
                yTmp[i] = yi[indexesTmp[i]]
    return sum(yTmp * coeffs)

def gaussInterpVector(x, a, b, xi, yi, NOrder, CYCLIC, int_method, IS_THETA, yi_OPP_side, INCLUDE_THETA_BOUNDARIES):
    N = x.shape[0]
    y = zeros(N, Complex)
    for i in range(N):
        y[i] = gaussInterp(x[i], a, b, xi, yi, NOrder, CYCLIC, int_method, IS_THETA, yi_OPP_side, INCLUDE_THETA_BOUNDARIES)
    return y

def gaussInterpMatrix(x, axi, bxi, xi, NOrderX, CYCLIC_X, int_method_x, y, ayi, byi, yi, NOrderY, CYCLIC_Y, int_method_y, Zi):
    Nx, Nxi, Ny, Nyi = x.shape[0], xi.shape[0], y.shape[0], yi.shape[0]
    print("Ny =", Ny)
    Z = zeros((Nx, Ny), Complex)
    Ztmp = zeros((Nxi, Ny), Complex)
    for i in range(Nxi):
        IS_THETA = 0
        Ztmp[i,:] = gaussInterpVector(y, ayi, byi, yi, Zi[i,:], NOrderY, CYCLIC_Y, int_method_y, IS_THETA, Zi[i,:], 0)
    for j in range(Ny):
        IS_THETA = 1
        INCLUDE_THETA_BOUNDARIES = 0
        if (xi[0]==axi)&(xi[-1]==bxi):
            INCLUDE_THETA_BOUNDARIES = 1
        if (j+Ny/2<Ny):
            if (INCLUDE_THETA_BOUNDARIES == 1):
                Z[:,j] = gaussInterpVector(x, axi, bxi, xi, Ztmp[:,j], NOrderX, CYCLIC_X, int_method_x, IS_THETA, Ztmp[1:-1,j+Ny/2], INCLUDE_THETA_BOUNDARIES)
            else:
                Z[:,j] = gaussInterpVector(x, axi, bxi, xi, Ztmp[:,j], NOrderX, CYCLIC_X, int_method_x, IS_THETA, Ztmp[:,j+Ny/2], INCLUDE_THETA_BOUNDARIES)
        else:
            if (INCLUDE_THETA_BOUNDARIES == 1):
                Z[:,j] = gaussInterpVector(x, axi, bxi, xi, Ztmp[:,j], NOrderX, CYCLIC_X, int_method_x, IS_THETA, Ztmp[1:-1,j-Ny/2], INCLUDE_THETA_BOUNDARIES)
            else:
                Z[:,j] = gaussInterpVector(x, axi, bxi, xi, Ztmp[:,j], NOrderX, CYCLIC_X, int_method_x, IS_THETA, Ztmp[:,j-Ny/2], INCLUDE_THETA_BOUNDARIES)
    return Z


def Lagrange_regular_matrix(n):
    M = zeros((n+1, n+1), Float)
    wrapping_code = """Lagrangian_regular_matrix(n, M);"""
    weave.inline(wrapping_code,
                 ['n', 'M'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    return M

def LagrangeVariableStepInterpolation(x, y_i, x_i):
    coeffs = zeros(x_i.shape[0], Float)
    wrapping_code = """LagrangeVariableStepInterpolationCoeffs(coeffs, x, x_i);"""
    weave.inline(wrapping_code,
                 ['coeffs', 'x', 'x_i'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    #print("Lagrange coeffs =", coeffs)
    return sum(y_i * coeffs)


def Lagrange_fixed_step_point_interpolation(x, xi, yi, n):
    N = xi.shape[0]
    d = xi[1] - xi[0]
    i = int(floor((x-xi[0])/d))
    ind = i + arange(n+1) - n/2
    if ind[0]<0:
        ind -= ind[0]
    elif ind[-1]>N-1:
        ind -= ind[-1]-(N-1)
    s = (x-xi[ind[0]])/d
    S = array([1.0, s, s**2, s**3, s**4])
    M = Lagrange_regular_matrix(n)
    result = sum(take(yi, ind) * matrixmultiply(M, S[0:n+1]))
    return result

def Lagrange_fixed_step_vector_interpolation(x, xi, yi, n):
    y = zeros(x.shape[0], Complex)
    for j in range(x.shape[0]):
        y[j] = Lagrange_fixed_step_point_interpolation(x[j], xi, yi, n)
    return y

def Lagrange_fixed_step_vector_interpolation_C(x, xi, yi, n):
    y = zeros(x.shape[0], Complex)
    M = Lagrange_regular_matrix(n)
    wrapping_code = """Lagrange_vector_fixedstep_interpolation(y, x, xi, yi, M);"""
    weave.inline(wrapping_code,
                 ['y', 'x', 'xi', 'yi', 'M'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    return y

def decimate_2D_Lagrange(theta_i, phi_i, Y_i, n):
    Y = zeros((2*theta_i.shape[0]-1, 2*phi_i.shape[0]-1), Complex)
    wrapping_code = """decimate_2D(Y, Y_i, theta_i, phi_i, n);"""
    weave.inline(wrapping_code,
                 ['Y', 'Y_i', 'theta_i', 'phi_i', 'n'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    return Y


def decimate_2D_lfi(theta, theta_i, A_theta, B_theta, NorderTheta, PERIODIC_Theta, CYCLIC_Theta, phi, phi_i, A_phi, B_phi, NorderPhi, PERIODIC_Phi, CYCLIC_Phi, Y_i):
    Ntheta_i, Nphi_i = Y_i.shape[0], Y_i.shape[1]
    Y_i_reshape = zeros(Ntheta_i * Nphi_i, 'F')
    for i in range(Nphi_i):
        Y_i_reshape[i*Ntheta_i:i*Ntheta_i+Ntheta_i] = Y_i[:,i]
    Y_reshape = zeros(theta.shape[0] * phi.shape[0], 'F')
    INCLUDED_BOUNDARIES_THETA = (abs(theta_i[0]-A_theta)<=1.e-6) & (abs(theta_i[-1]-B_theta)<=1.e-6)
    INCLUDED_BOUNDARIES_PHI = (abs(phi_i[0]-A_phi)<=1.e-6) & (abs(phi_i[-1]-B_phi)<=1.e-6)
    wrapping_code = """
    using namespace blitz;
    LagrangeFastInterpolator2D lfi2D(theta, theta_i, static_cast<float>(A_theta), static_cast<float>(B_theta), INCLUDED_BOUNDARIES_THETA, NorderTheta, PERIODIC_Theta, CYCLIC_Theta, phi, phi_i, static_cast<float>(A_phi), static_cast<float>(B_phi), INCLUDED_BOUNDARIES_PHI, NorderPhi, PERIODIC_Phi, CYCLIC_Phi);
    //cout << lfi2D.getCoefficientsForLinesInterp() << endl;
    //cout << lfi2D.getIndexesForLinesInterp() << endl;
    interpolate2Dlfi(Y_reshape, Y_i_reshape, lfi2D);
    """
    weave.inline(wrapping_code,
                 ['Y_reshape', 'Y_i_reshape', 'theta', 'theta_i', 'A_theta', 'B_theta', 'INCLUDED_BOUNDARIES_THETA', 'NorderTheta', 'PERIODIC_Theta', 'CYCLIC_Theta', 'phi', 'phi_i', 'A_phi', 'B_phi', 'INCLUDED_BOUNDARIES_PHI', 'NorderPhi', 'PERIODIC_Phi', 'CYCLIC_Phi'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    Y = zeros((theta.shape[0], phi.shape[0]), 'F')
    for i in range(Y.shape[1]):
        Y[:,i] = Y_reshape[i*Y.shape[0]:i*Y.shape[0]+Y.shape[0]]
    return Y


def decimate_2D_splint(theta_i, phi_i, Y_i):
    Y = zeros((2*theta_i.shape[0]-1, 2*phi_i.shape[0]-1), Complex)
    theta_int = zeros(2*theta_i.shape[0]-1, Float)
    phi_int = zeros(2*phi_i.shape[0]-1, Float)
    wrapping_code = """decimate_2D_splint(theta_int, phi_int, Y, theta_i, phi_i, Y_i);"""
    weave.inline(wrapping_code,
                 ['theta_int', 'phi_int', 'Y', 'theta_i', 'phi_i', 'Y_i'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"splint.h"'],
                 compiler = 'gcc')
    return Y

def anterDecimate2DLfi(theta, theta_i, NOrderTheta, PERIODIC_Theta, CYCLIC_Theta, phi, phi_i, NOrderPhi, PERIODIC_Phi, CYCLIC_Phi, Y_i):
    Y = zeros((theta_i.shape[0], phi_i.shape[0]), Complex)
    wrapping_code = """
    LagrangeFastInterpolator2D lfi2D(theta, theta_i, NorderTheta, PERIODIC_Theta, CYCLIC_Theta, phi, phi_i, NorderPhi, PERIODIC_Phi, CYCLIC_Phi);
    anterDecimate2DLfi(Y, Y_i, lfi2D);"""
    weave.inline(wrapping_code,
                 ['Y', 'Y_i', 'theta', 'theta_i', 'NOrderTheta', 'PERIODIC_Theta', 'CYCLIC_Theta', 'phi', 'phi_i', 'NOrderPhi', 'PERIODIC_Phi', 'CYCLIC_Phi'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"interpolation.h"'],
                 compiler = 'gcc')
    return Y

def testAnterpolation(N):
    Ntheta, Nphi = 2*N-1, 2*N-1
    Y = zeros((Ntheta, Nphi), Complex)
    Z = zeros((Ntheta, Nphi), Complex)
    W = zeros((Ntheta, Nphi), Complex)
    Xtheta, Wtheta = integr_1D_X_W(0.0, 5.0, Ntheta, "TRAP")
    Xphi, Wphi = integr_1D_X_W(0.0, 5.0, Nphi, "TRAP")
    for i in range(Ntheta):
        Y[i] = Xtheta[i]**(2./5.) * Xphi**(3./5.)
        W[i] = Wtheta[i] * Wphi
    Z = exp(-1.j*Y)
    gplt.surf(real(Z*Y))
    I = sum(sum(W*Y*Z))

    W_a = anterDecimate2DLfi(Xtheta, Xtheta[0::2], 4, Xphi, Xphi[0::2], 4, W)
    print(W_a[4,::2])
    print(W[8,::4])
    Y_a, Z_a = zeros(W_a.shape, Complex), zeros(W_a.shape, Complex)
    Y_a, Z_a = Y[0::2,0::2], Z[0::2,0::2]
    I_a = sum(sum(Y_a*Z_a*W_a))
    print(I, I_a, sum(sum(Y_a*Z[0::2,0::2]*W[0::2,0::2])))

if __name__=="__main__":
    Dtheta = pi/10
    Dphi = 2*pi/40
    theta = arange(0, pi+Dtheta, Dtheta)
    phi = arange(0, 2*pi, Dphi) # we don't include 2*pi, as in testMLFMAinterp.py
    Ntheta, Nphi = theta.shape[0], phi.shape[0]
    print("Ntheta, Nphi =", Ntheta, Nphi)
    Y_original = zeros((Ntheta, Nphi), Complex)
    for i in range(Nphi):
        Y_original[:,i] = cos(theta) * sin(phi[i])
    # lagrange variable step interpolation test
    #x_i = phi_i
    #y_i = sin(x_i)
    #y = LagrangeVariableStepInterpolation(x_i[11], y_i[6:10], x_i[6:10])
    #print("interpolated y =", y, ", real y =", y_i[11])

    # 1-D interpolation
##     y_interp = Lagrange_vector_interpolation(X_i, X_i[0::4], Y_i[0::4], 3)
##     y_interp_C = Lagrange_vector_interpolation_C(X_i, X_i[0::4], Y_i[0::4], 4)
##     #gplt.plot(X_i, Y_i, X_i[0::4], Y_i[0::4], X_i, y_interp, X_i, y_interp_C)
##     gplt.plot(X_i, real(Y_i), X_i, real(y_interp_C))


