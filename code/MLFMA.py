from mpi4py import *
import sys, time, copy, commands, pickle, cPickle
from numpy import round, ceil
from scipy import rand, prod, mean, arccos, dot, product, log, log10, ceil, floor, where, real, sqrt
from meshClass import MeshClass, CubeClass
from integration import *
from EM_constants import *
from MoMPostProcessing import *
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk, readIntFromDisk, readFloatFromDisk, read1DBlitzArrayFromDisk, readASCIIBlitzComplexFloatArray2DFromDisk

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
    #    print "N_points_max =", N_points_max
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
            print "Error: integration method must be GAUSSL, TRAP or PONCELET."
    return octtreeX, octtreeW, octtreeN

def computeTreeParameters(my_id, tmpDirName, a, k, N_levels, params_simu):
    # L computation
    NB_DIGITS = params_simu.NB_DIGITS
    L = zeros(N_levels-1, 'i') # array of poles numbers: 1 number per level
    for i in range(L.shape[0]):
        L[i]  = L_computation(k, a*(2**i), NB_DIGITS)
    if (my_id==0) and (params_simu.VERBOSE == 1):
        print "L =", L
    # integration and interpolation data
    octtreeXcosThetas, octtreeWthetas, octtreeNthetas = octtreeXWN_computation(-1.0, 1.0, L, N_levels, params_simu.int_method_theta, params_simu.INCLUDE_BOUNDARIES)
    octtreeXthetas = zeros(octtreeXcosThetas.shape, 'd')
    for i in range(octtreeNthetas.shape[0]):
        Npoints = octtreeNthetas[i]
        octtreeXthetas[i,:Npoints] = arccos(octtreeXcosThetas[i, Npoints-1::-1])
    octtreeXphis, octtreeWphis, octtreeNphis = octtreeXWN_computation(0.0, 2.0*pi, L, N_levels, params_simu.int_method_phi, params_simu.INCLUDE_BOUNDARIES)
    #if (my_id==0):
    #    print "Nthetas =", octtreeNthetas
    #    print "Nphis =", octtreeNphis
    NOrderInterpTheta = L[0]
    NOrderInterpPhi = L[0]
    # now we write the info to disk
    writeScalarToDisk(NOrderInterpTheta, os.path.join(tmpDirName, 'octtree_data/NOrderInterpTheta.txt') )
    writeScalarToDisk(NOrderInterpPhi, os.path.join(tmpDirName, 'octtree_data/NOrderInterpPhi.txt') )
    writeASCIIBlitzArrayToDisk(L, os.path.join(tmpDirName, 'octtree_data/LExpansion.txt') )
    # we need the following line as a workaround for a weave compilation issue
    # that arises sometimes in FMM_precond.py, when not all the weave modules have been compiled
    writeASCIIBlitzArrayToDisk(zeros((2,2), 'i'), os.path.join(tmpDirName, 'octtree_data/nothing.txt') )
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
    writeScalarToDisk(A_theta, os.path.join(tmpDirName, 'octtree_data/A_theta.txt') )
    writeScalarToDisk(B_theta, os.path.join(tmpDirName, 'octtree_data/B_theta.txt') )
    writeScalarToDisk(A_phi, os.path.join(tmpDirName, 'octtree_data/A_phi.txt') )
    writeScalarToDisk(B_phi, os.path.join(tmpDirName, 'octtree_data/B_phi.txt') )
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
    if not params_simu.AUTOMATIC_THETAS and (params_simu.USER_DEFINED_NB_THETA > 0):
        NpointsTheta = params_simu.USER_DEFINED_NB_THETA
    else:
        NpointsTheta *= int(ceil((params_simu.STOP_THETA - params_simu.START_THETA)/pi))
    octtreeXthetas_coarsest = zeros(NpointsTheta, 'd')
    if not NpointsTheta==1:
        DTheta = (params_simu.STOP_THETA - params_simu.START_THETA)/(NpointsTheta - 1)
        for i in range(NpointsTheta):
            octtreeXthetas_coarsest[i] = params_simu.START_THETA + i*DTheta
        # make sure the last element is params_simu.STOP_THETA
        octtreeXthetas_coarsest[-1] = params_simu.STOP_THETA
    else:
        octtreeXthetas_coarsest[0] = params_simu.START_THETA
    # phis abscissas
    NpointsPhi = 2 * L_coarsest
    if not params_simu.AUTOMATIC_PHIS and (params_simu.USER_DEFINED_NB_PHI > 0):
        NpointsPhi = params_simu.USER_DEFINED_NB_PHI
    else:
        NpointsPhi *= int(ceil((params_simu.STOP_PHI - params_simu.START_PHI)/(2.0*pi)))
    octtreeXphis_coarsest = zeros(NpointsPhi, 'd')
    if not NpointsPhi==1:
        DPhi = (params_simu.STOP_PHI - params_simu.START_PHI)/(NpointsPhi-1)
        for i in range(NpointsPhi):
            octtreeXphis_coarsest[i] = params_simu.START_PHI + i*DPhi
        # make sure the last element is params_simu.STOP_PHI
        octtreeXphis_coarsest[-1] = params_simu.STOP_PHI
    else:
        octtreeXphis_coarsest[0] = params_simu.START_PHI
    if (my_id==0) and (params_simu.VERBOSE == 1):
        print "L for coarsest level =", L_coarsest
    writeASCIIBlitzArrayToDisk(octtreeXthetas_coarsest, os.path.join(tmpDirName, 'octtree_data/octtreeXthetas_coarsest.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXphis_coarsest, os.path.join(tmpDirName, 'octtree_data/octtreeXphis_coarsest.txt') )
    MPI.COMM_WORLD.Barrier()

def print_times(params_simu, simuDirName):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    if (my_id == 0):
        # final moves
        file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
        variables = cPickle.load(file)
        file.close()
        numberOfMatvecs = readIntFromDisk(os.path.join(tmpDirName,'iterative_data/numberOfMatvecs.txt'))
        NIterMLFMA = readIntFromDisk(os.path.join(tmpDirName,'iterative_data/iter.txt'))
        average_RWG_length = readFloatFromDisk(os.path.join(tmpDirName,'mesh/average_RWG_length.txt'))
        # CPU_time_GMSH
        file = open(os.path.join(simuDirName,'result/CPU_time_GMSH.txt'), 'r')
        CPU_time_GMSH_tmp = file.readlines()
        file.close()
        CPU_time_GMSH = 0.0
        for line in CPU_time_GMSH_tmp:
            if 'real' in line:
                CPU_time_GMSH = float(line.split()[1])
        # CPU_time_distribute_Z_cubes
        file = open(os.path.join(simuDirName,'result/CPU_time_distribute_Z_cubes.txt'), 'r')
        CPU_time_distribute_Z_cubes_tmp = file.readlines()
        file.close()
        CPU_time_distribute_Z_cubes = 0.0
        for line in CPU_time_distribute_Z_cubes_tmp:
            if 'real' in line:
                CPU_time_distribute_Z_cubes = float(line.split()[1])
        # CPU_time_communicateZnearBlocks
        file = open(os.path.join(simuDirName,'result/CPU_time_communicateZnearBlocks.txt'), 'r')
        CPU_time_communicateZnearBlocks_tmp = file.readlines()
        file.close()
        CPU_time_communicateZnearBlocks = 0.0
        for line in CPU_time_communicateZnearBlocks_tmp:
            if 'real' in line:
                CPU_time_communicateZnearBlocks = float(line.split()[1])
        # CPU_time_MLFMA
        file = open(os.path.join(simuDirName,'result/CPU_time_MLFMA.txt'), 'r')
        CPU_time_MLFMA_tmp = file.readlines()
        file.close()
        CPU_time_MLFMA = 0.0
        for line in CPU_time_MLFMA_tmp:
            if 'real' in line:
                CPU_time_MLFMA = float(line.split()[1])
        if (params_simu.VERBOSE == 1):
            print "N RWG =", variables['N_RWG']
            sys.stdout.write("average RWG length = %.5s" %str(average_RWG_length) + " m = lambda / %.9s" %str((c/params_simu.f)/average_RWG_length) + " \n")
            print CPU_time_GMSH, "CPU time (seconds) for GMSH meshing."
            print CPU_time_distribute_Z_cubes, "CPU time (seconds) for distribute_Z_cubes."
            print variables['CPU_time_Z_near_computation'], "CPU time (seconds) for constructing Z_CFIE_near"
            print variables['Wall_time_Z_near_computation'], "Wall time (seconds) for constructing Z_CFIE_near"
            print CPU_time_communicateZnearBlocks, "CPU time (seconds) for communicateZnearBlocks."
            print variables['CPU_time_Mg_computation'], "CPU time (seconds) for constructing SAI precond"
            print variables['Wall_time_Mg_computation'], "Wall time (seconds) for constructing SAI precond"
            print CPU_time_MLFMA, "CPU time (seconds) for MLFMA iterations and solution. Iterations =", NIterMLFMA
            #print target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for MLFMA iterations and solution. Iterations =", target_MLFMA.NIterMLFMA
            if numberOfMatvecs>0:
                print CPU_time_MLFMA/numberOfMatvecs, "CPU time (seconds) per MLFMA matvec"
                #print target_MLFMA.Wall_time_Target_MLFMA_resolution/target_MLFMA.numberOfMatvecs, "Wall time (seconds) per MLFMA matvec"
            print CPU_time_GMSH + CPU_time_distribute_Z_cubes + variables['CPU_time_Z_near_computation'] + CPU_time_communicateZnearBlocks + variables['CPU_time_Mg_computation'] + CPU_time_MLFMA, "CPU time (seconds) for complete MLFMA solution"
            #print Wall_time_Z_near_computation + Wall_time_Mg_computation + target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for complete MLFMA solution"
        if params_simu.CURRENTS_VISUALIZATION:
            computeCurrentsVisualization(params_simu, variables, simuDirName)
        if params_simu.SAVE_CURRENTS_CENTROIDS:
            saveCurrentsCentroids(params_simu, simuDirName)

def computeCurrentsVisualization(params_simu, variables, simuDirName):
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    if (my_id == 0):
        target_mesh = MeshClass(geoDirName, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join(tmpDirName, "mesh"))
        CURRENTS_VISUALIZATION = params_simu.CURRENTS_VISUALIZATION and (target_mesh.N_RWG<2e6) and (params_simu.BISTATIC == 1)
        if CURRENTS_VISUALIZATION:
            print "............Constructing visualisation of currents"
            if (target_mesh.N_RWG<1e4):
                nbTimeSteps = 48
            else:
                nbTimeSteps = 1
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            norm_J_centroids_triangles = normJMCentroidsTriangles(J_centroids_triangles, variables['w'], nbTimeSteps)
            write_VectorFieldTrianglesCentroidsGMSH(os.path.join(geoDirName, params_simu.targetName) + '.J_centroids_triangles.pos', real(J_centroids_triangles), target_mesh)
            write_ScalarFieldTrianglesCentroidsGMSH(os.path.join(geoDirName, params_simu.targetName) + '.norm_J_centroids_triangles.pos', norm_J_centroids_triangles, target_mesh)
            print "............end of currents visualisation construction."
        else:
            "Error in the computeCurrentsVisualization routine. Probably you required currents visualization computation for a mesh too big"
            sys.exit(1)

def saveCurrentsCentroids(params_simu, simuDirName):
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    if (my_id == 0):
        target_mesh = MeshClass(geoDirName, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join(tmpDirName, "mesh"))
        SAVE_CURRENTS_CENTROIDS = (params_simu.SAVE_CURRENTS_CENTROIDS and (params_simu.BISTATIC == 1))
        if SAVE_CURRENTS_CENTROIDS:
            print "............saving currents at centroids of triangles"
            nbTimeSteps = 1
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            write_VectorFieldTrianglesCentroids(os.path.join(simuDirName, 'result', 'J_centroids_triangles.txt'), J_centroids_triangles, target_mesh)
            print "............end of saving of currents."

def getField(filename):
    my_id = MPI.COMM_WORLD.Get_rank()
    E_obs = zeros(3, 'F')
    if (my_id == 0):
        E_obs = readASCIIBlitzComplexFloatArray2DFromDisk(filename)
    return E_obs


