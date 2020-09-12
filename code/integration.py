from math import pi
from scipy import zeros, sqrt, sort, real, sum, cos


def integr_1D_X_W(a, b, N_points, METHOD, INCLUDE_BOUDARIES):
    """This function yields the abscissas and weights for 3 integration methods:
    PONCELET (mid-point), TRAPEZOIDAL RULE, GAUSS. The boolean INCLUDE_BOUDARIES specifies if one
    must take the extremities of the interval into account."""
    X, W = zeros(N_points, 'd'), zeros(N_points, 'd')
    if (METHOD == "TRAP"): # trapezoidal rule
        h = (b-a)/(N_points-1)
        for j in range(N_points):
            X[j] = a + j*h
            W[j] = h
        W[0] /= 2.0
        W[N_points-1] /= 2.0
    elif (METHOD == "PONCELET"): # mid-point method
        h = (b-a)/N_points
        for j in range(N_points):
            X[j] = a + j*h + h/2.0
            W[j] = h
    else:
        Dx = 0.5 * (b - a)
        center = 0.5 * (b + a)
        #XGL = real(special.legendre(N_points).weights[:,0])
        #WGL = real(special.legendre(N_points).weights[:,1])
        XGL, WGL = largeXGLWGL(N_points)
        X = center + Dx * XGL
        W = abs(Dx) * WGL
    if ( INCLUDE_BOUDARIES & (METHOD=='GAUSSL') ):
        Xtmp = zeros(N_points+2, 'd')
        Wtmp = zeros(N_points+2, 'd')
        Xtmp[1:-1], Wtmp[1:-1] = X, W
        Xtmp[0], Xtmp[-1] = a, b
        Wtmp[0], Wtmp[-1] = 0.0, 0.0
        X, W = Xtmp, Wtmp
    return X, W


def mlegzo(n):
    """I found this function on the Internet (search for MLEGZO)"""
    n0 = (n+1)/2
    x, w = zeros(n), zeros(n)
    for nr in range(1, int(n0+1)):
        z = cos(pi * (nr-0.25)/n)
        z0 = 0
        while (abs(z-z0) > abs(z) * 1.0e-15) and (z!=0):
            z0 = z
            p = 1.0
            for i in range(0, nr-2):
                p *= (z-x[i])
            f0 = 1.0
            if (nr==n0) and (n!=2*n/2):
                z = 0.0
            f1 = z
            for k in range(2, n+1):
                pf = (2.0 - 1.0/k)*z*f1-(1.0-1.0/k)*f0
                pd = k*(f1-z*pf)/(1.0-z*z)
                f0 = f1
                f1 = pf
            if z==0.0:
                break
            fd = pf/p
            q = 0.0
            for i in range(1, nr):
                wp = 1.0
                for j in range(1, nr):
                    if j!=i:
                        wp = wp*(z-x[j-1])
                q += wp
            gd = (pd - q * fd)/p
            z -= fd/gd
        x[nr-1] = z
        x[n-nr] = -z
        w[nr-1] = 2.0/((1.0-z*z)*pd*pd)
        w[n-nr] = w[nr-1]
    return sort(x, kind='mergesort'), w


def largeXGLWGL(N_points):
    """find Gauss-Legendre abscissas and weights for large N. See page 39 of Guillaume Sylvand thesis"""
    X, W = mlegzo(N_points)
    return X, W


if __name__=="__main__":

    n = 268
    X1, W1 = largeXGLWGL(n)
    X2, W2 = mlegzo(n)
    print(max(abs(X1 - X2)))
    print(max(abs(W1 - W2)))
    n = n-1
    X1, W1 = largeXGLWGL(n)
    X2, W2 = mlegzo(n)
    print(max(abs(X1 - X2)))
    print(max(abs(W1 - W2)))
