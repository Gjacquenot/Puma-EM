from scipy import *
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    pass
from pylab import *
import copy

def norm(x):
    return sqrt(abs(sum(x * conjugate(x))))

def rotmat(a, b):
    # compute the Givens rotation matrix parameters for a and b
    if (b==0.0):
        c, s = 1.0, 0.0
    elif (abs(b) > abs(a)):
        ##temp = a/b
        temp = -a/b
        s = 1.0 / sqrt(1.0 + temp**2)
        c = temp * s
    else:
        ##temp = b/a
        temp = -b/a
        c = 1.0 / sqrt(1.0 + temp**2)
        s = temp * c
    return c, s

def myGmresPython( A, x, b, M, restrt, max_it, tol ):

    iter = 0                                         # initialization
    flag = 0

    bnrm2 = norm(b)
    if  ( bnrm2 == 0.0 ):
        bnrm2 = 1.0

    #r = M \ ( b-A*x );
    r = b - dot(A, x)
    error = norm( r ) / bnrm2
    if ( error < tol ):
        return x, error, iter, flag

    n = A.shape[0]                                  # initialize workspace
    m = restrt
    V = zeros((n,m+1), Complex)
    H = zeros((m+1, m), Complex)
    cs, sn = zeros(m, Complex), zeros(m, Complex)
    e1 = zeros(m+1, Complex)
    e1[0] = 1.0

    for iter in range(max_it):
        #r = M \ ( b-A*x )
        r = b - dot(A, x)
        V[:, 0] = r / norm( r )
        #print('sum(V[:, 0]) =', sum(V[:, 0]))
        s = norm( r )*e1
        # construct orthonormal basis using Gram-Schmidt
        for i in range(m):
            #w = M \ matrixmultiply(A, V[:, i])
            w = dot(A, V[:, i])
            for k in range(i+1):
                H[k, i] = sum(w * conjugate(V[:, k]))
                w = w - H[k, i] * V[:, k]
            #print("sum(w) =", sum(w))
            H[i+1, i] = norm( w )
            V[:, i+1] = w / H[i+1, i]
            #print("sum(V[:, i+1]) =", sum(V[:, i+1]))
            # apply Givens rotation
            for k in range(i):
                ##temp = cs[k] * H[k, i] + sn[k] * H[k+1, i]
                ##H[k+1,i] = -sn[k] * H[k, i] + cs[k] * H[k+1, i]
                temp = conjugate(cs[k]) * H[k, i] - conjugate(sn[k]) * H[k+1, i]
                H[k+1,i] = sn[k] * H[k, i] + cs[k] * H[k+1, i]
                H[k,i] = temp
            # form i-th rotation matrix
            cs[i], sn[i] = rotmat( H[i, i], H[i+1, i] )
            ##H[i, i] = cs[i]*H[i, i] + sn[i]*H[i+1, i]
            H[i, i] = conjugate(cs[i])*H[i, i] - conjugate(sn[i])*H[i+1, i]
            H[i+1,i] = 0.0
            #print("sum(sum(H)) =", sum(sum(H)))
            # approximate residual norm
            ##temp = cs[i] * s[i]
            temp = conjugate(cs[i]) * s[i] - conjugate(sn[i]) * s[i+1]
            ##s[i+1] = -sn[i] * s[i]
            s[i+1] = sn[i] * s[i] + cs[i] * s[i+1]
            s[i] = temp;
            error = abs(s[i+1]) / bnrm2
            #print("error =", error)
            # update approximation
            if ( error <= tol ):
                y = linalg.solve( H[:i+1, :i+1], s[:i+1] )
                x = x + dot(V[:, :i+1], y)
                return x, error, iter, flag

        if ( error <= tol ):
            return x, error, iter, flag
        y = linalg.solve( H[:m, :m], s[:m] )
        #print("sum(y) =", sum(y))
        x = x + dot(V[:, :m], y)          # update approximation
        #print("sum(x) =", sum(x))
        #r = M \ ( b-A*x )                           # compute residual
        r = b - dot(A, x)
        error = abs(s[i+1]) / bnrm2                  # check convergence
        #print("error =", error)
        if ( error <= tol ):
            return x, error, iter, flag

    if ( error > tol ):
        flag = 1
    return x, error, iter, flag


def myGmresC( A, x, b, M, restrt, max_it, tol ):
    error = zeros(1, Float)
    flag = zeros(1, Int)
    iter = zeros(1, Int)
    wrapping_code = """
    // complex<double>
    classicalMatrixMatvec< std::complex<float> > Aptr(A);
    //MatvecFunctor< std::complex<float> > matvec(&Aptr, &classicalMatrixMatvec< std::complex<float> >::classicalMatvec);
    MatvecFunctor< std::complex<float>, classicalMatrixMatvec< std::complex<float> > > matvec(&Aptr, &classicalMatrixMatvec< std::complex<float> >::classicalMatvec);
    gmres(x, error(0), iter(0), flag(0), matvec, b, tol, restrt, max_it);
    """
    weave.inline(wrapping_code,
                 ['x', 'error', 'iter', 'flag', 'A', 'b', 'tol', 'restrt', 'max_it'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/iterative/'],
                 library_dirs = ['./MoM/iterative/'],
                 libraries = ['ITERATIVE'],
                 headers = ['<iostream>','<complex>','"iterative.h"'],
                 compiler = 'gcc')
    return x, error[0], iter[0], flag[0]

def myBiCGSTAB(A, x, b, M, max_it, tol):
    iter = 0                                         # initialization
    flag = 0

    bnrm2 = norm(b)
    if  ( bnrm2 == 0.0 ):
        bnrm2 = 1.0
    r = b - dot(A, x)
    error = norm( r ) / bnrm2
    if ( error < tol ):
        return x, error, iter, flag

    omega  = 1.0
    r_tld = copy.deepcopy(r);
    for iter in range(1, max_it):
        rho = dot(conjugate(r_tld), r)
        print(rho)
        if (rho==0.0):
            return x, error, iter, flag
        if iter==1:
            p = copy.deepcopy(r)
        else:
            beta  = ( rho/rho_1 )*( alpha/omega );
            p = r + beta*( p - omega*v );
        # p_hat = M \ p
        p_hat = copy.deepcopy(p)
        v = dot(A, p_hat)
        alpha = rho / dot( conjugate(r_tld), v );
        s = r - alpha*v
        if ( norm(s) < tol ):
            x = x + alpha*p_hat
            error = norm( s ) / bnrm2
            return x, error, iter, flag
        # stabilizer
        #s_hat = M \ s
        s_hat = copy.deepcopy(s)
        t = dot(A, s_hat)
        omega = dot( conjugate(t), s) / dot( conjugate(t), t )
        # update approximation
        x = x + alpha*p_hat + omega*s_hat

        r = s - omega*t
        error = norm( r ) / bnrm2
        if ( error <= tol ):
            return x, error, iter, flag
        if ( omega == 0.0 ):
            return x, error, iter, flag
        rho_1 = rho

    if ( error <= tol ) or ( norm( s ) <= tol ): # converged
        if ( norm( s ) <= tol ):
            error = norm(s) / bnrm2s
        flag =  0
    elif ( omega == 0.0 ): # breakdown
        flag = -2
    elif ( rho == 0.0 ):
        flag = -1
    else: # no convergence
        flag = 1;
    return x, error, iter, flag


def myBiCGSTABC( A, x, b, M, max_it, tol ):
    error = zeros(1, Float)
    flag = zeros(1, Int)
    iter = zeros(1, Int)
    wrapping_code = """
    // complex<double>
    classicalMatrixMatvec< std::complex<float> > Aptr(A);
    //MatvecFunctor< std::complex<float> > matvec(&Aptr, &classicalMatrixMatvec< std::complex<float> >::classicalMatvec);
    MatvecFunctor< std::complex<float>, classicalMatrixMatvec< std::complex<float> > > matvec(&Aptr, &classicalMatrixMatvec< std::complex<float> >::classicalMatvec);
    bicgstab(x, error(0), iter(0), flag(0), matvec, psolve, b, tol, max_it);
    """
    weave.inline(wrapping_code,
                 ['x', 'error', 'iter', 'flag', 'A', 'b', 'tol', 'max_it'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/iterative/'],
                 library_dirs = ['./MoM/iterative/'],
                 libraries = ['ITERATIVE'],
                 headers = ['<iostream>','<complex>','"iterative.h"'],
                 compiler = 'gcc')
    return x, error[0], iter[0], flag[0]


if __name__=="__main__":
    testBidon()
    N = 20
    Z = rand(N,N)+ 1.j*rand(N,N)
    V = rand(N) + 1.j*rand(N)
    ELEM_TYPE = 'F'
    Z = Z.astype(ELEM_TYPE)
    V = V.astype(ELEM_TYPE)

    Y = zeros(N, ELEM_TYPE)
    X0 = zeros(N, ELEM_TYPE)
    M = ones((N,N), ELEM_TYPE)
    restart = N
    tol = 1.e-3
    maxiter = 50
    # scipy linalg
    I = matrixmultiply(linalg.inv(Z), V)
    # myGmresPython
    I2, error2, iteration2, info2 = myGmresPython( Z, X0, V, M, restart, maxiter, tol )
    print("error2, iteration2, info2 =", error2, iteration2, info2)
    # myGmresC++
    X0 = zeros(N, ELEM_TYPE)
    I3, error3, iteration3, info3 = myGmresC( Z, X0, V, M, restart, maxiter, tol )
    print("error3, iteration3, info3 =", error3, iteration3, info3)

    print(sum(I), sum(I2), sum(I3))

