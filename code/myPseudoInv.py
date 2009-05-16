import sys, os
from scipy import conj, array, zeros, rand, eye, dot, floor, transpose, real, imag, triu
from scipy import weave
from scipy.weave import converters
#from zgels import zgels, dlamch
import copy

def macheps():
    """computes the machine epsilon"""
    machEps, tmp = 1.0, 1.0
    while tmp + machEps != 1.0:
        machEps /= 2.0
    return 2.0*machEps

def computeLwork(M, N, nrhs):
    """computes the space necessary for WORK to.....work!!
    returns the dimension of the array WORK.
    LWORK >= max( 1, MN + max( MN, NRHS ) ).
    For optimal performance,
    LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
    where MN = min(M,N) and NB is the optimum block size.
    """
    NB = 20
    mn = min(M, N)
    return max( 1, mn + max( mn, nrhs )*NB )

def computeMyPinv(A):
    """this routine computes the pseudo inverse of A
    and is a wrapper to the fortran function zgels.f
    By the way, to understand this wrapping structure, you
    better check the comments at the beginning of zgels.f"""
    m, n = A.shape
    lda = A.shape[0]
    ldb = max(n, m)
    nrhs = m
    B = zeros((ldb, nrhs), complex)
    for i in range(min(ldb, nrhs)):
        B[i, i] += 1.
    lwork = computeLwork(m, n, nrhs)
    work = zeros(lwork, complex)
    trans = 'N' # 'C'
    AA = copy.deepcopy(A)
    out = zgels(trans, nrhs, AA, B, lwork)
    info = out[-1]
    if info==0:
        return out[0][:n] # see comments on B in "zgels.f"
    else:
        print "argument", abs(i), "has illegal value. exiting"
        sys.exit(1)

def computeMyPinvCC(A, LIB_G2C):
    """this routine computes the pseudo inverse of A
    and is a wrapper to the fortran function zgels.f
    By the way, to understand this wrapping structure, you
    better check the comments at the beginning of zgels.f"""
    m, n = A.shape
    lda = A.shape[0]
    ldb = max(n, m)
    nrhs = m
    B = zeros((ldb, nrhs), 'D', 2)
    #for i in range(min(ldb, nrhs)):
    #    B[i, i] += 1.
    lwork = computeLwork(m, n, nrhs)
    work = zeros(lwork, 'D')
    trans = 'N' # 'C'
    info = 0
    AA = zeros(A.shape, 'D', 2)
    AA[:] = A
    wrapping_code = """
    char trans = 'N';
    int N = min(ldb, nrhs);
    for (int i=0 ; i<N ; ++i) B(i, i) = 1.;
    zgels(trans, m, n, nrhs, AA, lda, B, ldb, work, lwork, info);
    """
    weave.inline(wrapping_code,
                ['m', 'n', 'nrhs', 'AA', 'lda', 'B', 'ldb', 'work', 'lwork', 'info'],
                type_converters = converters.blitz,
                include_dirs = ['./code/MoM/lapack/', '.'],
                library_dirs = ['./code/MoM/lapack/', '.'],
                libraries = [LIB_G2C, 'm', 'ZGELS'],
                headers = ['<iostream>','<complex>','<blitz/array.h>', '"zgels_interface.h"'],
                compiler = 'gcc',
                extra_compile_args = ['-O3', '-pthread', '-w'])
    if info==0:
        if (m<=n):
            return B
        else:
            return B[:n]
    else:
        print "argument", abs(i), "has illegal value. exiting"
        sys.exit(1)

def computeTriangleUpSolve(A, LIB_G2C):
    """this routine computes the pseudo inverse of A
    and is a wrapper to the fortran function ztrsm.f
    By the way, to understand this wrapping structure, you
    better check the comments at the beginning of zgels.f"""
    m, n = A.shape
    lda = A.shape[0]
    ldb = max(n, m)
    nrhs = m
    N = min(ldb, nrhs)
    B = zeros((ldb, nrhs), 'D', 2)
    #for i in range(min(ldb, nrhs)):
    #    B[i, i] += 1.
    lwork = computeLwork(m, n, nrhs)
    work = zeros(lwork, 'D')
    alpha = 1.0 + 0.j
    info = 0
    AA = zeros(A.shape, 'D', 2)
    AA[:] = A
    wrapping_code = """
    char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
    for (int i=0 ; i<N ; ++i) B(i, i) = 1.;
    ztrsm(side, uplo, transa, diag, m, n, alpha, AA, lda, B, ldb);
    """
    weave.inline(wrapping_code,
                ['m', 'n', 'alpha', 'AA', 'lda', 'B', 'ldb', 'N'],
                type_converters = converters.blitz,
                include_dirs = ['./code/MoM/lapack/', '.'],
                library_dirs = ['./code/MoM/lapack/', '.'],
                libraries = [LIB_G2C, 'm', 'ZGELS'],
                headers = ['<iostream>','<complex>','<blitz/array.h>', '"ztrsm_interface.h"'],
                compiler = 'gcc',
                extra_compile_args = ['-O3', '-pthread', '-w'])
    if info==0:
        if (m<=n):
            return B
        else:
            return B[:n]
    else:
        print "argument", abs(i), "has illegal value. exiting"
        sys.exit(1)

if __name__=="__main__":
    machEps = macheps()
    print "Python eps =", machEps
    #print "FORTRAN eps =", dlamch('E')
    #print "FORTRAN safe min =", dlamch('S')
    #print "FORTRAN base of the machine =", dlamch('B')
    #print "FORTRAN eps*base =", dlamch('P')
    #print "FORTRAN number of (base) digits in the mantissa =", dlamch('N')
    #print "FORTRAN rmax = ", dlamch('O')

    N = 3
    M = 3
    A0 = triu(rand(M, N) + 1.j * rand(M, N))
    A1 = copy.deepcopy(A0)
    A2 = zeros((M, N), 'D', 2)
    A2[:] = A0
    print A0
    #for j in range(A1.shape[0]):
        #A2[j] = A0[j]
    LIB_G2C = 'g2c' # for gcc >= 4.3.0, LIB_G2C = 'gfortran'
    from scipy import linalg
    X1 = linalg.pinv(A1) # the pseudo inverse
    X2 = computeMyPinvCC(A2, LIB_G2C)
    X3 = computeTriangleUpSolve(A1, LIB_G2C)
    print
    print "X1 - X2 = "
    print X1 -X2
    print "X1 - X3 = "
    print X1 - X3
