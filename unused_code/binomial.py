from scipy import *
from weave import converters

def factorial(N):
    """Simple factorial function"""

    wrapping_code = """
	 return_val = factorial (N);
    """

    F = weave.inline(wrapping_code,
                 ['N'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"binomial.h"'],
                 compiler = 'gcc')

    return F

def factorial_N_k(N, k):
    """Simple N-k factorial function"""

    wrapping_code = """
	 return_val = factorial_N_k (N, k);
    """

    F = weave.inline(wrapping_code,
                 ['N', 'k'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"binomial.h"'],
                 compiler = 'gcc')

    return F

def binomial(N, k):
    """Simple N-k binomial function"""

    wrapping_code = """
	 return_val = binomial (N, k);
    """

    F = weave.inline(wrapping_code,
                 ['N', 'k'],
                 type_converters = converters.blitz,
                 include_dirs = ['./MoM/'],
                 library_dirs = ['./MoM/'],
                 libraries = ['MoM'],
                 headers = ['<iostream>','<complex>','"binomial.h"'],
                 compiler = 'gcc')

    return F


if __name__=="__main__":

    for N in range(14):
	for k in range(N):
	    print N, k, factorial(N)/factorial(N-k)/factorial(k), factorial_N_k(N, k)/factorial(k), binomial(N, k)


