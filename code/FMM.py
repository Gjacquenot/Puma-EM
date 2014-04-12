from mpi4py import *
from scipy import log, log10, ceil, floor, arccos, product, where, real, sqrt
from EM_constants import *

def L_computation(k, a, NB_DIGITS):
    """this function computes the number of expansion poles given the wavenumber k and the sidelength a"""
    L_tmp1 = a*real(k)*sqrt(3) + NB_DIGITS * (a*real(k)*sqrt(3))**(1./3.)
    L_tmp2 = a*real(k)*sqrt(3) + NB_DIGITS * log10(a*real(k)*sqrt(3) + pi)
    L_tmp3 = a*real(k)*sqrt(3) + 1 * log(a*real(k)*sqrt(3) + pi)
    L_tmp4 = a*real(k)*sqrt(3) + 1.8 * NB_DIGITS**(2.0/3.0) * (a*real(k)*sqrt(3))**(1.0/3.0)
    #my_id = MPI.COMM_WORLD.Get_rank()
    #if (my_id==0):
        #print "L_tmp1 =", L_tmp1, ", L_tmp2 =", L_tmp2, ", L_tmp3 =", L_tmp3, ", L_tmp4 =", L_tmp4
    L2 = where(L_tmp2<5., 5, int(ceil( L_tmp2 )))
    #return L2
    return int(floor( L_tmp2 ))
    #return int(ceil( L_tmp2 ))

