from Numeric import *
from scipy import *
from scipy import weave
from scipy import linalg
from weave import converters
import time

def multiplyAndSum(A, B):
    code = """ 
           double result;
           result = sum(A*B) + sum(A+B);
           """
    weave.inline(code,
                ['A', 'B'],
                type_converters = converters.blitz,
                include_dirs = ['./MoM/'],
                library_dirs = ['./MoM/'],
                libraries = ['MoM'],
                headers = ['<iostream>','<complex>'],
                compiler = 'gcc')

if __name__=="__main__":
    N = 500
    A = rand(2, N, 2*N)
    B = rand(2, N, 2*N)
    t0 = time.time()
    multiplyAndSum(A, B)
    print time.time() - t0, "seconds..."

    A = rand(2, 2*N*N)
    B = rand(2, 2*N*N)
    t0 = time.time()
    multiplyAndSum(A, B)
    print time.time() - t0, "seconds..."

    A = rand(2*N*N, 2)
    B = rand(2*N*N, 2)
    t0 = time.time()
    multiplyAndSum(A, B)
    print time.time() - t0, "seconds..."
