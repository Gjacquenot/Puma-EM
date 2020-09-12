from mpi4py import MPI
import sys, os, argparse
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from scipy import arccos, log, log10, ceil, floor, where, real, sqrt
from integration import *
from EM_constants import *
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk

def L_computation(k, a, NB_DIGITS):
    """this function computes the number of expansion poles given the wavenumber k and the sidelength a"""
    L_tmp1 = a*real(k)*sqrt(3) + NB_DIGITS * (a*real(k)*sqrt(3))**(1./3.)
    L_tmp2 = a*real(k)*sqrt(3) + NB_DIGITS * log10(a*real(k)*sqrt(3) + pi) # Guillaume SYLVAND: p. 129
    L_tmp3 = a*real(k)*sqrt(3) + 1 * log(a*real(k)*sqrt(3) + pi)
    L_tmp4 = a*real(k)*sqrt(3) + 1.8 * NB_DIGITS**(2.0/3.0) * (a*real(k)*sqrt(3))**(1.0/3.0) # Song, Chew, "Fast and Efficient..."
    #my_id = MPI.COMM_WORLD.Get_rank()
    #if (my_id==0):
        #print("L_tmp1 =", L_tmp1, ", L_tmp2 =", L_tmp2, ", L_tmp3 =", L_tmp3, ", L_tmp4 =", L_tmp4)
    L2 = where(L_tmp2<5., 5, int(ceil( L_tmp2 )))
    #return L2
    return int(floor( L_tmp4 ))

def octtreeXWN_computation(boundary_inf, boundary_sup, L, N_levels, int_method, INCLUDE_BOUNDARIES):
    """This function computes the abscissas and weights for all the levels.
    We do this in python because we don't have the appropriate function in C/C++
    for finding the Gauss-Legendre weights and abscissas."""
    my_id = MPI.COMM_WORLD.Get_rank()
    if int_method in ["GAUSSL"]:
        N_points_max = L[-1] + 1
    else:
        N_points_max = 2*L[-1]
    INCLUDE_BOUNDARIES_2 = ( (int_method=="GAUSSL") & INCLUDE_BOUNDARIES )
    if (INCLUDE_BOUNDARIES_2):
        # we also include 0 and pi at the extremities for enhancing the interpolation
        # this is especially true for gauss-legendre abscissas with non-cyclic interpolation
        N_points_max += 2
    #if (my_id==0):
    #    print("N_points_max =", N_points_max)
    octtreeX = zeros((N_levels-1, N_points_max), 'd')
    octtreeW = zeros((N_levels-1, N_points_max), 'd')
    octtreeN = zeros(N_levels-1, 'i')
    for i in range(N_levels-1):
        N = N_points_max
        if int_method in ["GAUSSL"]:
            N = L[i] + 1
        else:
            N = 2*L[i]
        if int_method in ["GAUSSL", "PONCELET"]:
            octtreeXTmp, octtreeWTmp = integr_1D_X_W(boundary_inf, boundary_sup, N, int_method, INCLUDE_BOUNDARIES_2)
            octtreeN[i] = octtreeXTmp.shape[0]
            octtreeX[i,:octtreeN[i]], octtreeW[i,:octtreeN[i]] = octtreeXTmp, octtreeWTmp
        elif int_method=="TRAP":
            octtreeXTmp, octtreeWTmp = integr_1D_X_W(boundary_inf, boundary_sup, N, "PONCELET", INCLUDE_BOUNDARIES_2)
            octtreeN[i] = octtreeXTmp.shape[0]
            octtreeX[i,:octtreeN[i]], octtreeW[i,:octtreeN[i]] = octtreeXTmp, octtreeWTmp
            octtreeX[i,:octtreeN[i]] -= (octtreeX[i,1] - octtreeX[i,0])/2.
        else:
            print("Error: integration method must be GAUSSL, TRAP or PONCELET.")
    return octtreeX, octtreeW, octtreeN

def directions_zones_calculation(P):
    HZ, VZ = 1, P
    if P<8:
        HZ, VZ = 1, P
    elif (P%2==0):
        HALF_P = P/2
        diviseur = int(floor(sqrt(HALF_P)))
        while diviseur>0:
            if (HALF_P%diviseur==0):
                HZ, VZ = int(diviseur), int(P/diviseur)
                break
            else:
                diviseur -= 1
        if diviseur == 1:
            HZ, VZ = 2, P/2
    else:
        diviseur = int(floor(sqrt(P)))
        while diviseur>0:
            if (P%diviseur==0):
                HZ, VZ = int(diviseur), int(P/diviseur)
                break
            else:
                diviseur -= 1
    return HZ, VZ

def computeTreeParameters(my_id, tmpDirName, a, k, N_levels, params_simu):
    # L computation
    NB_DIGITS = params_simu.NB_DIGITS
    L = zeros(N_levels-1, 'i') # array of poles numbers: 1 number per level
    for i in range(L.shape[0]):
        L[i]  = L_computation(k, a*(2**i), NB_DIGITS)
    if (my_id==0) and (params_simu.VERBOSE == 1):
        print("L = " + str(L))
    # integration and interpolation data
    octtreeXcosThetas, octtreeWthetas, octtreeNthetas = octtreeXWN_computation(-1.0, 1.0, L, N_levels, params_simu.int_method_theta, params_simu.INCLUDE_BOUNDARIES)
    octtreeXthetas = zeros(octtreeXcosThetas.shape, 'd')
    for i in range(octtreeNthetas.shape[0]):
        Npoints = octtreeNthetas[i]
        octtreeXthetas[i,:Npoints] = arccos(octtreeXcosThetas[i, Npoints-1::-1])
    octtreeXphis, octtreeWphis, octtreeNphis = octtreeXWN_computation(0.0, 2.0*pi, L, N_levels, params_simu.int_method_phi, params_simu.INCLUDE_BOUNDARIES)
    #if (my_id==0):
    #    print("Nthetas =", octtreeNthetas)
    #    print("Nphis =", octtreeNphis)
    # order of interpolation
    NOrderInterpTheta = L[0]
    NOrderInterpPhi = L[0]
    # number of zones per theta
    #num_proc = MPI.COMM_WORLD.Get_size()
    #Ntheta_zones, Nphi_zones = directions_zones_calculation(num_proc)
    # now we write the info to disk
    writeScalarToDisk(NOrderInterpTheta, os.path.join(tmpDirName, 'octtree_data/NOrderInterpTheta.txt') )
    writeScalarToDisk(NOrderInterpPhi, os.path.join(tmpDirName, 'octtree_data/NOrderInterpPhi.txt') )
    #writeScalarToDisk(Ntheta_zones, os.path.join(tmpDirName, 'octtree_data/Ntheta_zones.txt') )
    #writeScalarToDisk(Nphi_zones, os.path.join(tmpDirName, 'octtree_data/Nphi_zones.txt') )
    writeASCIIBlitzArrayToDisk(L, os.path.join(tmpDirName, 'octtree_data/LExpansion.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_smoothing_factor, os.path.join(tmpDirName, 'octtree_data/alphaTranslation_smoothing_factor.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_thresholdRelValueMax, os.path.join(tmpDirName, 'octtree_data/alphaTranslation_thresholdRelValueMax.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_RelativeCountAboveThreshold, os.path.join(tmpDirName, 'octtree_data/alphaTranslation_RelativeCountAboveThreshold.txt') )
    writeASCIIBlitzArrayToDisk(octtreeNthetas, os.path.join(tmpDirName, 'octtree_data/octtreeNthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeNphis, os.path.join(tmpDirName, 'octtree_data/octtreeNphis.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXthetas, os.path.join(tmpDirName, 'octtree_data/octtreeXthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXphis, os.path.join(tmpDirName, 'octtree_data/octtreeXphis.txt') )
    writeASCIIBlitzArrayToDisk(octtreeWthetas, os.path.join(tmpDirName, 'octtree_data/octtreeWthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeWphis, os.path.join(tmpDirName, 'octtree_data/octtreeWphis.txt') )
    A_theta, B_theta, A_phi, B_phi = 0., pi, 0., 2.*pi
    N_theta, N_phi = octtreeNthetas[0], octtreeNphis[0]
    INCLUDED_THETA_BOUNDARIES, INCLUDED_PHI_BOUNDARIES = 0, 0
    if (abs(octtreeXthetas[0,0]-A_theta)<=1.e-8) and (abs(octtreeXthetas[0,N_theta-1]-B_theta)<=1.e-8):
        INCLUDED_THETA_BOUNDARIES = 1
    if (abs(octtreeXphis[0,0]-A_phi)<=1.e-8) and (abs(octtreeXphis[0,N_phi-1]-B_phi)<=1.e-8):
        INCLUDED_PHI_BOUNDARIES = 1
    writeScalarToDisk(INCLUDED_THETA_BOUNDARIES, os.path.join(tmpDirName, 'octtree_data/INCLUDED_THETA_BOUNDARIES.txt') )
    writeScalarToDisk(INCLUDED_PHI_BOUNDARIES, os.path.join(tmpDirName, 'octtree_data/INCLUDED_PHI_BOUNDARIES.txt') )

    # we now have to calculate the theta/phi abscissas for the coarsest level
    # These are needed for far-field computation
    L_coarsest = L_computation(k, a*(2**N_levels), NB_DIGITS)
    # theta abscissas
    NpointsTheta = L_coarsest + 1
    DTheta = 0
    if not params_simu.AUTOMATIC_THETAS and (params_simu.USER_DEFINED_NB_THETA > 0):
        NpointsTheta = params_simu.USER_DEFINED_NB_THETA
    else:
        NpointsThetaTmp = NpointsTheta * (params_simu.STOP_THETA - params_simu.START_THETA)/pi
        NpointsTheta = int(ceil(NpointsThetaTmp))+1
    octtreeXthetas_coarsest = zeros(NpointsTheta, 'd')
    if NpointsTheta>1:
        DTheta = (params_simu.STOP_THETA - params_simu.START_THETA)/(NpointsTheta - 1)
        for i in range(NpointsTheta):
            octtreeXthetas_coarsest[i] = params_simu.START_THETA + i*DTheta
        # make sure the last element is params_simu.STOP_THETA
        octtreeXthetas_coarsest[-1] = params_simu.STOP_THETA
    else:
        octtreeXthetas_coarsest[0] = params_simu.START_THETA
    # phis abscissas
    NpointsPhi = 2 * L_coarsest
    DPhi = 0
    if not params_simu.AUTOMATIC_PHIS and (params_simu.USER_DEFINED_NB_PHI > 0):
        NpointsPhi = params_simu.USER_DEFINED_NB_PHI
    else:
        NpointsPhiTmp = NpointsPhi * (params_simu.STOP_PHI - params_simu.START_PHI)/(2.0*pi)
        NpointsPhi = int(ceil(NpointsPhiTmp))+1
    octtreeXphis_coarsest = zeros(NpointsPhi, 'd')
    if NpointsPhi>1:
        DPhi = (params_simu.STOP_PHI - params_simu.START_PHI)/(NpointsPhi-1)
        for i in range(NpointsPhi):
            octtreeXphis_coarsest[i] = params_simu.START_PHI + i*DPhi
        # make sure the last element is params_simu.STOP_PHI
        octtreeXphis_coarsest[-1] = params_simu.STOP_PHI
    else:
        octtreeXphis_coarsest[0] = params_simu.START_PHI
    if (my_id==0):
        print("Summary of sampling points at the coarsest level (used for far-field sampling).")
        print("L_coarsest =", L_coarsest)
        print("For 0 < theta < 180, NpointsTheta = L_coarsest + 1 =", L_coarsest + 1)
        print("For", params_simu.START_THETA/pi*180, "< theta <", params_simu.STOP_THETA/pi*180, ", NpointsTheta =", NpointsTheta, ", DTheta =", DTheta/pi*180, "degrees")
        print("For 0 < phi < 360, NpointsPhi = 2 * L_coarsest =", 2 * L_coarsest)
        print("For", params_simu.START_PHI/pi*180, "< phi <", params_simu.STOP_PHI/pi*180, ", NpointsPhi =", NpointsPhi, ", DPhi =", DPhi/pi*180, "degrees")
    writeASCIIBlitzArrayToDisk(octtreeXthetas_coarsest, os.path.join(tmpDirName, 'octtree_data/octtreeXthetas_coarsest.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXphis_coarsest, os.path.join(tmpDirName, 'octtree_data/octtreeXphis_coarsest.txt') )
    MPI.COMM_WORLD.Barrier()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--inputdir')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    inputDirName = cmdline.inputdir
    simuParams = 'simulation_parameters'

    # the simulation itself
    my_id = MPI.COMM_WORLD.Get_rank()
    if (my_id==0):
        sys.path.append(os.path.abspath(inputDirName))
        exec('from ' + simuParams + ' import *')
    else:
        params_simu = ['blabla']
    params_simu = MPI.COMM_WORLD.bcast(params_simu)
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        my_id = MPI.COMM_WORLD.Get_rank()
        tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
        file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'rb')
        variables = cPickle.load(file)
        file.close()
        computeTreeParameters(my_id, tmpDirName, variables['a'], variables['k'], variables['N_levels'], params_simu)
    else:
        print("you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings.")
        sys.exit(1)


