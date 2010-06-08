from mpi4py import *
import sys, time, copy, commands, pickle, cPickle
from numpy import round, ceil
from scipy import rand, prod, mean, arccos, dot, product, log, log10, ceil, floor, where, real, sqrt
from meshClass import MeshClass, CubeClass
from PyGmsh import executeGmsh, write_geo
from FMM_precond import MgPreconditionerComputation, Mg_CSR, Mg_listsOfZnearBlocks_ToTransmitAndReceive
from FMM_Znear import Z_nearCRS_Computation, Z_nearCRS_Assembling, Z_nearChunksDistribution, Z_near_size_computation
from integration import *
from EM_constants import *
from MoMPostProcessing import *
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk, writeBlitzArrayToDisk, readIntFromDisk, readFloatFromDisk, read1DBlitzArrayFromDisk, readASCIIBlitzComplexFloatArray2DFromDisk, readASCIIBlitzFloatArray2DFromDisk
from runMPIsystemCommand import runMPIsystemCommand, createMPIsystemCommand
from read_dipole_excitation import read_dipole_excitation, read_observation_points

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
    writeScalarToDisk(NOrderInterpTheta, os.path.join('.', tmpDirName, 'octtree_data/NOrderInterpTheta.txt') )
    writeScalarToDisk(NOrderInterpPhi, os.path.join('.', tmpDirName, 'octtree_data/NOrderInterpPhi.txt') )
    writeASCIIBlitzArrayToDisk(L, os.path.join('.', tmpDirName, 'octtree_data/LExpansion.txt') )
    # we need the following line as a workaround for a weave compilation issue
    # that arises sometimes in FMM_precond.py, when not all the weave modules have been compiled
    writeASCIIBlitzArrayToDisk(zeros((2,2), 'i'), os.path.join('.', tmpDirName, 'octtree_data/nothing.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_smoothing_factor, os.path.join('.', tmpDirName, 'octtree_data/alphaTranslation_smoothing_factor.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_thresholdRelValueMax, os.path.join('.', tmpDirName, 'octtree_data/alphaTranslation_thresholdRelValueMax.txt') )
    writeScalarToDisk(params_simu.alphaTranslation_RelativeCountAboveThreshold, os.path.join('.', tmpDirName, 'octtree_data/alphaTranslation_RelativeCountAboveThreshold.txt') )
    writeASCIIBlitzArrayToDisk(octtreeNthetas, os.path.join('.', tmpDirName, 'octtree_data/octtreeNthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeNphis, os.path.join('.', tmpDirName, 'octtree_data/octtreeNphis.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXthetas, os.path.join('.', tmpDirName, 'octtree_data/octtreeXthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXphis, os.path.join('.', tmpDirName, 'octtree_data/octtreeXphis.txt') )
    writeASCIIBlitzArrayToDisk(octtreeWthetas, os.path.join('.', tmpDirName, 'octtree_data/octtreeWthetas.txt') )
    writeASCIIBlitzArrayToDisk(octtreeWphis, os.path.join('.', tmpDirName, 'octtree_data/octtreeWphis.txt') )
    A_theta, B_theta, A_phi, B_phi = 0., pi, 0., 2.*pi
    writeScalarToDisk(A_theta, os.path.join('.', tmpDirName, 'octtree_data/A_theta.txt') )
    writeScalarToDisk(B_theta, os.path.join('.', tmpDirName, 'octtree_data/B_theta.txt') )
    writeScalarToDisk(A_phi, os.path.join('.', tmpDirName, 'octtree_data/A_phi.txt') )
    writeScalarToDisk(B_phi, os.path.join('.', tmpDirName, 'octtree_data/B_phi.txt') )
    N_theta, N_phi = octtreeNthetas[0], octtreeNphis[0]
    INCLUDED_THETA_BOUNDARIES, INCLUDED_PHI_BOUNDARIES = 0, 0
    if (abs(octtreeXthetas[0,0]-A_theta)<=1.e-8) and (abs(octtreeXthetas[0,N_theta-1]-B_theta)<=1.e-8):
        INCLUDED_THETA_BOUNDARIES = 1
    if (abs(octtreeXphis[0,0]-A_phi)<=1.e-8) and (abs(octtreeXphis[0,N_phi-1]-B_phi)<=1.e-8):
        INCLUDED_PHI_BOUNDARIES = 1
    writeScalarToDisk(INCLUDED_THETA_BOUNDARIES, os.path.join('.', tmpDirName, 'octtree_data/INCLUDED_THETA_BOUNDARIES.txt') )
    writeScalarToDisk(INCLUDED_PHI_BOUNDARIES, os.path.join('.', tmpDirName, 'octtree_data/INCLUDED_PHI_BOUNDARIES.txt') )

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
    writeASCIIBlitzArrayToDisk(octtreeXthetas_coarsest, os.path.join('.', tmpDirName, 'octtree_data/octtreeXthetas_coarsest.txt') )
    writeASCIIBlitzArrayToDisk(octtreeXphis_coarsest, os.path.join('.', tmpDirName, 'octtree_data/octtreeXphis_coarsest.txt') )
    MPI.COMM_WORLD.Barrier()

def distributeZcubesAttributions(MAX_BLOCK_SIZE, N_nearPerCube, C, tmpDirName, params_simu):
    my_id = MPI.COMM_WORLD.Get_rank()
    num_proc = MPI.COMM_WORLD.Get_size()
    if (my_id == 0):
        print "distributing the chunks and cubes among processes..."
    writeScalarToDisk( num_proc, os.path.join('.', tmpDirName,'octtree_data/num_procs.txt') )
    ZchunkNumber_to_cubesNumbers, ZcubeNumber_to_chunkNumber, ZchunkNumber_to_processNumber, ZprocessNumber_to_ChunksNumbers = Z_nearChunksDistribution(MAX_BLOCK_SIZE, N_nearPerCube, C, tmpDirName, params_simu)
    ZchunkNumber_to_cubesNumbers = MPI.COMM_WORLD.Bcast(ZchunkNumber_to_cubesNumbers)
    ZcubeNumber_to_chunkNumber = MPI.COMM_WORLD.Bcast(ZcubeNumber_to_chunkNumber)
    ZchunkNumber_to_processNumber = MPI.COMM_WORLD.Bcast(ZchunkNumber_to_processNumber)
    ZprocessNumber_to_ChunksNumbers = MPI.COMM_WORLD.Bcast(ZprocessNumber_to_ChunksNumbers)
    if (my_id == 0):
        print "Process", my_id, ": self.ZprocessNumber_to_ChunksNumbers =", ZprocessNumber_to_ChunksNumbers
    CUBES_DISTRIBUTION = 0 # the cubes distribution is finished
    writeScalarToDisk(CUBES_DISTRIBUTION, os.path.join('.', tmpDirName, 'octtree_data/CUBES_DISTRIBUTION.txt') )
    return ZchunkNumber_to_cubesNumbers, ZcubeNumber_to_chunkNumber, ZchunkNumber_to_processNumber, ZprocessNumber_to_ChunksNumbers

def scatterMesh(target_mesh, ZprocessNumber_to_ChunksNumbers, ZchunkNumber_to_cubesNumbers, tmpDirName, my_id, num_proc):
    # we now exchange the local (cubes) meshes data...
    CPU_time, Wall_time = time.clock(), time.time()
    pathToSaveTo = os.path.join('.', tmpDirName, 'Z_tmp')
    if (my_id == 0 ):
        print "Exchanging local meshes-cubes data..."
    # I create the necessary directories for Z_tmp
    for identity in range(num_proc):
        for chunkNumber in ZprocessNumber_to_ChunksNumbers[identity]:
            if (my_id == 0): # master node
                list_cubeIntArrays, list_cubeDoubleArrays = [], []
            else:
                list_cubeIntArrays, list_cubeDoubleArrays = ['blabla'], ['blabla']
            for cubeNumber in ZchunkNumber_to_cubesNumbers[chunkNumber]:
                # for each one we compute the necessary information for individual Zcube computation
                if (my_id == 0): # master node
                    cubeIntArrays, cubeDoubleArrays = target_mesh.computeCubeLocalArrays(cubeNumber)
                    # we append the cube lists to new lists
                    list_cubeIntArrays.append(cubeIntArrays)
                    list_cubeDoubleArrays.append(cubeDoubleArrays)
            # communicating the arrays
            # we exchange the concatenated arrays
            list_cubeIntArrays = MPI.COMM_WORLD.Bcast(list_cubeIntArrays)
            list_cubeDoubleArrays = MPI.COMM_WORLD.Bcast(list_cubeDoubleArrays)
            # writing the local cube data to disk
            if (my_id==identity):
                pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
                os.mkdir(pathToSaveToChunk)
                for j in range(len(ZchunkNumber_to_cubesNumbers[chunkNumber])):
                    cubeNumber = ZchunkNumber_to_cubesNumbers[chunkNumber][j]
                    cube = CubeClass()
                    cube.cubeIntArrays = list_cubeIntArrays[j]
                    cube.cubeDoubleArrays = list_cubeDoubleArrays[j]
                    cube.writeIntDoubleArraysToFile(pathToSaveToChunk, cubeNumber)

    CPU_time, Wall_time = time.clock() - CPU_time, time.time() - Wall_time
    print "mesh scattering: CPU time =", CPU_time, "sec"
    print "mesh scattering: Wall time =", Wall_time, "sec"

def computeZ_near(tmpDirName, a, C, k, w, N_levels, chunkNumber_to_cubesNumbers, processNumber_to_ChunksNumbers, CFIE, params_simu):
    my_id = MPI.COMM_WORLD.Get_rank()
    num_proc = MPI.COMM_WORLD.Get_size()
    # computation of near interactions matrix
    ELEM_TYPE = 'F'
    Z_TMP_ELEM_TYPE = 'F'
    Wall_t0 = time.time()
    CPU_t0 = time.clock()
    pathToSaveTo = os.path.join('.', tmpDirName, 'Z_tmp')
    # then I can compute the chunks
    if (my_id == 0):
        print "Z_CFIE_near computation..."
    MPI.COMM_WORLD.Barrier()
    Z_nearCRS_Computation(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, CFIE, params_simu.MAX_BLOCK_SIZE, w, params_simu.eps_r, params_simu.mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, params_simu.TDS_APPROX, params_simu.Z_s, params_simu.MOM_FULL_PRECISION, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()
    CPU_time_Z_near_computation = time.clock() - CPU_t0
    Wall_time_Z_near_computation = time.time() - Wall_t0
    return Wall_time_Z_near_computation, CPU_time_Z_near_computation

def compute_SAIpreconditioner(tmpDirName, a, C, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers, MAX_BLOCK_SIZE):
    my_id = MPI.COMM_WORLD.Get_rank()
    num_proc = MPI.COMM_WORLD.Get_size()
    # computation of near interactions matrix
    ELEM_TYPE = 'F'
    Z_TMP_ELEM_TYPE = 'F'
    # computation of the Frobenius preconditioner
    Wall_t0 = time.time()
    CPU_t0 = time.clock()
    pathToReadFrom = os.path.join('.', tmpDirName, 'Z_tmp')
    pathToSaveTo = os.path.join('.', tmpDirName, 'Mg_LeftFrob')
    # we look for the LIB_G2C type
    file = open('makefile.inc', 'r')
    content = file.readlines()
    file.close()
    for elem in content:
        if 'G2C' in elem:
            LIB_G2C = elem.split('=')[1].split()[0]
    if (my_id == 0):
        print "SAI preconditioner computation..."
    MPI.COMM_WORLD.Barrier()
    R_NORM_TYPE_1 = a
    Mg_CSR(my_id, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, a, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, LIB_G2C, pathToReadFrom, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()
    CPU_time_Mg_computation = time.clock() - CPU_t0 
    Wall_time_Mg_computation = time.time() - Wall_t0
    # assembling of near interactions matrix
    pathToReadFrom = os.path.join('.', tmpDirName, 'Z_tmp')
    pathToSaveTo = os.path.join('.', tmpDirName, 'Z_near')
    Z_nearCRS_Assembling(processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, MAX_BLOCK_SIZE, C, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom, pathToSaveTo)
    MPI.COMM_WORLD.Barrier()
    return Wall_time_Mg_computation, CPU_time_Mg_computation

def setup_MLFMA(params_simu):
    """This function provides an easy interface for running an MLFMA simulation.
       params_simu is a class instance that contains the parameters for the simulation.
       COMPUTE_Z_NEAR is a bool that tells if we have to recompute the near-field 
       matrix and the Sparse Approximate Inverse preconditioner.
    """
    status = MPI.Status()
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    targetFileName = os.path.join(params_simu.pathToTarget, params_simu.targetName + params_simu.meshFileTermination)
    tmpDirName = 'tmp' + str(my_id)

    if params_simu.COMPUTE_Z_NEAR==1:
        if params_simu.meshToMake:
            os.system("rm -f " + targetFileName)
            if my_id==0:
                write_geo(params_simu.pathToTarget, params_simu.targetName, 'lc', c/params_simu.f * params_simu.lc_factor)
                write_geo(params_simu.pathToTarget, params_simu.targetName, 'lx', params_simu.lx)
                write_geo(params_simu.pathToTarget, params_simu.targetName, 'ly', params_simu.ly)
                write_geo(params_simu.pathToTarget, params_simu.targetName, 'lz', params_simu.lz)
                executeGmsh(params_simu.pathToTarget, params_simu.targetName, 0)

        commands.getoutput("rm -rf ./tmp*")
        MPI.COMM_WORLD.Barrier()
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName)) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'Z_tmp')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'Z_near')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'Mg_LeftFrob')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'mesh')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'octtree_data')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'V_CFIE')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'ZI')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'iterative_data')) )
        os.mkdir( os.path.abspath(os.path.join('.',tmpDirName,'pickle')) )

    # target_mesh construction
    target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
    # size of cube at finest level
    a = c/params_simu.f * params_simu.a_factor
    MPI.COMM_WORLD.Barrier()
    if (my_id==0):
        if params_simu.COMPUTE_Z_NEAR==1:
            target_mesh.constructFromGmshFile()
        else:
            target_mesh.constructFromSavedArrays(os.path.join("./tmp" + str(my_id), "mesh"))
        average_RWG_length = target_mesh.average_RWG_length
        writeScalarToDisk(average_RWG_length, os.path.join('.',tmpDirName,'mesh/average_RWG_length.txt'))
        if params_simu.VERBOSE==1:
            print "average RWG length =", average_RWG_length, "m = lambda /", (c/params_simu.f)/average_RWG_length
        # target_mesh cubes computation
        target_mesh.cubes_data_computation(a)
        N_nearBlockDiag, target_mesh.N_near, target_mesh.N_nearPerCube  = Z_near_size_computation(target_mesh.cubes_lists_RWGsNumbers, target_mesh.cubesNeighborsIndexes)
        target_mesh.saveToDisk(os.path.join('.', tmpDirName,'mesh'))
    else:
        target_mesh.N_nearPerCube = ['blabla']
    target_mesh.N_nearPerCube = MPI.COMM_WORLD.Bcast(target_mesh.N_nearPerCube)
    if (my_id==0) and (params_simu.COMPUTE_Z_NEAR==1):
        runMPIsystemCommand("./code/MoM", "communicateMeshArrays", num_proc)
    MPI.COMM_WORLD.Barrier()

    N_levels = readIntFromDisk(os.path.join('.', tmpDirName,'mesh', 'N_levels.txt'))
    N_RWG = readIntFromDisk(os.path.join('.', tmpDirName,'mesh', 'E.txt'))
    C = readIntFromDisk(os.path.join('.', tmpDirName,'mesh', 'C.txt'))
    T = readIntFromDisk(os.path.join('.', tmpDirName,'mesh', 'T.txt'))
    if (my_id == 0):
        IS_CLOSED_SURFACE = target_mesh.IS_CLOSED_SURFACE
        big_cube_lower_coord = target_mesh.big_cube_lower_coord
        big_cube_center_coord = target_mesh.big_cube_center_coord
    else:
        IS_CLOSED_SURFACE = ['blabla']
        big_cube_lower_coord = ['blabla']
        big_cube_center_coord = ['blabla']
    big_cube_lower_coord = MPI.COMM_WORLD.Bcast(big_cube_lower_coord)
    big_cube_center_coord = MPI.COMM_WORLD.Bcast(big_cube_center_coord)
    IS_CLOSED_SURFACE = MPI.COMM_WORLD.Bcast(IS_CLOSED_SURFACE)

    CFIE = array([params_simu.nu, 0, 0, -(1.0-params_simu.nu) * 377]).astype('D')
    w = 2. * pi * params_simu.f
    k = w * sqrt(eps_0*params_simu.eps_r*mu_0*params_simu.mu_r) + 1.j * 0.
    # now the observation points
    if (params_simu.r_obs_FROM_FILE == 0):
        if (len(params_simu.r_obs_x)==len(params_simu.r_obs_y)) and (len(params_simu.r_obs_x)==len(params_simu.r_obs_z)):
            N_obs_points = len(params_simu.r_obs_x)
            r_obs = zeros((N_obs_points, 3), 'd')
            for index in range(N_obs_points):
                r_obs[index,:] = array([params_simu.r_obs_x[index], params_simu.r_obs_y[index], params_simu.r_obs_z[index]], 'd')
        else:
            if (my_id==0):
                print "Error in the r_obs parameter. Check your observation points!"
            sys.exit(1)
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join('.',tmpDirName,'V_CFIE/r_obs.txt'))
    elif (params_simu.r_obs_FROM_FILE == 1) and (params_simu.r_obs_FILENAME != ""):
        if (my_id==0): # this file is only on processor 0
            r_obs = read_observation_points(params_simu.r_obs_FILENAME)
        else:
            r_obs = zeros((1, 3), 'd')
        r_obs = MPI.COMM_WORLD.Bcast(r_obs)
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join('.',tmpDirName,'V_CFIE/r_obs.txt'))

    # now the excitations
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_DIPOLES, os.path.join('.',tmpDirName,'V_CFIE/DIPOLES_EXCITATION.txt'))
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_PLANE_WAVE, os.path.join('.',tmpDirName,'V_CFIE/PLANE_WAVE_EXCITATION.txt'))
    # if we have dipoles excitation AND definition of the excitation in simulation_parameters.py
    if (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_DIPOLES_FROM_FILE == 0):
        # first for electric dipoles
        N_src_points = 0
        if (len(params_simu.r_J_src_x)==len(params_simu.r_J_src_y)) and (len(params_simu.r_J_src_x)==len(params_simu.r_J_src_z)):
            N_src_points = len(params_simu.r_J_src_x)
            if len(params_simu.J_src_x)==N_src_points and len(params_simu.J_src_x)==len(params_simu.J_src_y) and len(params_simu.J_src_x)==len(params_simu.J_src_z):
                J_src = zeros((N_src_points, 3), 'D')
                r_J_src = zeros((N_src_points, 3), 'd')
                for index in range(N_src_points):
                    r_J_src[index,:] = array([params_simu.r_J_src_x[index], params_simu.r_J_src_y[index], params_simu.r_J_src_z[index]], 'd')
                    J_src[index,:] = array([params_simu.J_src_x[index], params_simu.J_src_y[index], params_simu.J_src_z[index]], 'D')
            else:
                if (my_id==0):
                    print "Error in the r_J_src parameter. Check your source points!"
                sys.exit(1)
        else:
            if (my_id==0):
                print "Error in the r_J_src parameter. Check your source points!"
            sys.exit(1)
        if N_src_points > 0:
            writeScalarToDisk(1, os.path.join('.',tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join('.',tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join('.',tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join('.',tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
        # then for magnetic dipoles
        N_src_points = 0
        if (len(params_simu.r_M_src_x)==len(params_simu.r_M_src_y)) and (len(params_simu.r_M_src_x)==len(params_simu.r_M_src_z)):
            N_src_points = len(params_simu.r_M_src_x)
            if len(params_simu.M_src_x)==N_src_points and len(params_simu.M_src_x)==len(params_simu.M_src_y) and len(params_simu.M_src_x)==len(params_simu.M_src_z):
                M_src = zeros((N_src_points, 3), 'D')
                r_M_src = zeros((N_src_points, 3), 'd')
                for index in range(N_src_points):
                    r_M_src[index,:] = array([params_simu.r_M_src_x[index], params_simu.r_M_src_y[index], params_simu.r_M_src_z[index]], 'd')
                    M_src[index,:] = array([params_simu.M_src_x[index], params_simu.M_src_y[index], params_simu.M_src_z[index]], 'D')
            else:
                if (my_id==0):
                    print "Error in the r_M_src parameter. Check your source points!"
                sys.exit(1)
        else:
            if (my_id==0):
                print "Error in the r_M_src parameter. Check your source points!"
            sys.exit(1)
        if N_src_points > 0:
            writeScalarToDisk(1, os.path.join('.',tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join('.',tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join('.',tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join('.',tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # if we have dipoles excitation AND definition of the excitation in a user-supplied file
    elif (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_DIPOLES_FROM_FILE == 1):
        if params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                J_src, r_J_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME)
            else:
                J_src, r_J_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            J_src = MPI.COMM_WORLD.Bcast(J_src)
            r_J_src = MPI.COMM_WORLD.Bcast(r_J_src)
            writeScalarToDisk(1, os.path.join('.',tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join('.',tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join('.',tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join('.',tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
        if params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                M_src, r_M_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME)
            else:
                M_src, r_M_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            M_src = MPI.COMM_WORLD.Bcast(M_src)
            r_M_src = MPI.COMM_WORLD.Bcast(r_M_src)
            writeScalarToDisk(1, os.path.join('.',tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join('.',tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join('.',tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join('.',tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # now the plane wave excitation
    if params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1:
        writeScalarToDisk(params_simu.theta_inc, os.path.join('.',tmpDirName,'V_CFIE/theta_inc.txt'))
        writeScalarToDisk(params_simu.phi_inc, os.path.join('.',tmpDirName,'V_CFIE/phi_inc.txt'))
        E_inc = array([params_simu.E_inc_theta, params_simu.E_inc_phi], 'D')
        writeASCIIBlitzArrayToDisk(E_inc, os.path.join('.',tmpDirName,'V_CFIE/E_inc.txt'))
    if (params_simu.BISTATIC_EXCITATION_DIPOLES != 1) and (params_simu.BISTATIC_EXCITATION_PLANE_WAVE != 1):
        if (my_id==0):
            print "incorrect excitation choice. You have to choose dipole and/or plane wave excitation."
        sys.exit(1)
    if params_simu.MONOSTATIC_SAR==1:
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_x_hat, 'd'), os.path.join('.',tmpDirName,'V_CFIE/SAR_local_x_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_y_hat, 'd'), os.path.join('.',tmpDirName,'V_CFIE/SAR_local_y_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_plane_origin, 'd'), os.path.join('.',tmpDirName,'V_CFIE/SAR_plane_origin.txt'))
        writeScalarToDisk(params_simu.SAR_x_span, os.path.join('.',tmpDirName,'V_CFIE/SAR_x_span.txt'))
        writeScalarToDisk(params_simu.SAR_y_span, os.path.join('.',tmpDirName,'V_CFIE/SAR_y_span.txt'))
        writeScalarToDisk(params_simu.SAR_x_span_offset, os.path.join('.',tmpDirName,'V_CFIE/SAR_x_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_y_span_offset, os.path.join('.',tmpDirName,'V_CFIE/SAR_y_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_N_x_points, os.path.join('.',tmpDirName,'V_CFIE/SAR_N_x_points.txt'))
        writeScalarToDisk(params_simu.SAR_N_y_points, os.path.join('.',tmpDirName,'V_CFIE/SAR_N_y_points.txt'))

    writeScalarToDisk( a, os.path.join('.',tmpDirName,'octtree_data/leaf_side_length.txt') )
    writeScalarToDisk(2.0*pi*params_simu.f, os.path.join('.',tmpDirName,'octtree_data/w.txt') )
    writeScalarToDisk(params_simu.eps_r, os.path.join('.',tmpDirName,'octtree_data/eps_r.txt') )
    writeScalarToDisk(params_simu.mu_r, os.path.join('.',tmpDirName,'octtree_data/mu_r.txt') )
    writeScalarToDisk(k, os.path.join('.',tmpDirName,'octtree_data/k.txt') )
    writeASCIIBlitzArrayToDisk(CFIE, os.path.join('.',tmpDirName,'octtree_data/CFIEcoeffs.txt') )
    writeScalarToDisk(N_RWG, os.path.join('.',tmpDirName,'octtree_data/N_RWG.txt') )
    writeScalarToDisk(N_levels-1, os.path.join('.',tmpDirName,'octtree_data/N_active_levels.txt') )
    writeASCIIBlitzArrayToDisk(big_cube_lower_coord, os.path.join('.',tmpDirName,'octtree_data/big_cube_lower_coord.txt') )
    writeASCIIBlitzArrayToDisk(big_cube_center_coord, os.path.join('.',tmpDirName,'octtree_data/big_cube_center_coord.txt') )
    writeScalarToDisk(params_simu.PERIODIC_Theta*1, os.path.join('.',tmpDirName,'octtree_data/PERIODIC_Theta.txt') )
    writeScalarToDisk(params_simu.CYCLIC_Theta*1, os.path.join('.',tmpDirName,'octtree_data/CYCLIC_Theta.txt') )
    writeScalarToDisk(params_simu.PERIODIC_Phi*1, os.path.join('.',tmpDirName,'octtree_data/PERIODIC_Phi.txt') )
    writeScalarToDisk(params_simu.CYCLIC_Phi*1, os.path.join('.',tmpDirName,'octtree_data/CYCLIC_Phi.txt') )
    writeScalarToDisk(params_simu.ALLOW_CEILING_LEVEL*1, os.path.join('.', tmpDirName, 'octtree_data/ALLOW_CEILING_LEVEL.txt') )
    writeScalarToDisk(params_simu.DIRECTIONS_PARALLELIZATION*1, os.path.join('.', tmpDirName, 'octtree_data/DIRECTIONS_PARALLELIZATION.txt') )
    writeScalarToDisk(params_simu.BE_BH_N_Gauss_points, os.path.join('.', tmpDirName, 'octtree_data/N_GaussOnTriangle.txt') )
    writeScalarToDisk(params_simu.MOM_FULL_PRECISION*1, os.path.join('.', tmpDirName, 'octtree_data/MOM_FULL_PRECISION.txt') )
    writeScalarToDisk(params_simu.VERBOSE*1, os.path.join('.', tmpDirName, 'octtree_data/VERBOSE.txt') )
    # what type of simulation are we running?
    writeScalarToDisk(params_simu.BISTATIC*1, os.path.join('.', tmpDirName, 'BISTATIC.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_RCS*1, os.path.join('.', tmpDirName, 'MONOSTATIC_RCS.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_SAR*1, os.path.join('.', tmpDirName, 'MONOSTATIC_SAR.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HH*1, os.path.join('.', tmpDirName, 'COMPUTE_RCS_HH.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_VV*1, os.path.join('.', tmpDirName, 'COMPUTE_RCS_VV.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HV*1, os.path.join('.', tmpDirName, 'COMPUTE_RCS_HV.txt') )
    writeScalarToDisk(params_simu.USE_PREVIOUS_SOLUTION*1, os.path.join('.', tmpDirName, 'USE_PREVIOUS_SOLUTION.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_BY_BISTATIC_APPROX*1, os.path.join('.', tmpDirName, 'MONOSTATIC_BY_BISTATIC_APPROX.txt') )
    writeScalarToDisk(params_simu.MAXIMUM_DELTA_PHASE, os.path.join('.', tmpDirName, 'MAXIMUM_DELTA_PHASE.txt') )
    # writing the iterative solver setup
    restrt = min(params_simu.RESTART, N_RWG)
    writeScalarToDisk(params_simu.MAXITER, os.path.join('.', tmpDirName, 'iterative_data/MAXITER.txt') )
    writeScalarToDisk(restrt, os.path.join('.', tmpDirName, 'iterative_data/RESTART.txt') )
    writeScalarToDisk(params_simu.SOLVER, os.path.join('.', tmpDirName, 'iterative_data/SOLVER.txt') )
    writeScalarToDisk(params_simu.INNER_SOLVER, os.path.join('.', tmpDirName, 'iterative_data/INNER_SOLVER.txt') )
    writeScalarToDisk(params_simu.TOL, os.path.join('.', tmpDirName, 'iterative_data/TOL.txt') )
    writeScalarToDisk(params_simu.INNER_TOL, os.path.join('.', tmpDirName, 'iterative_data/INNER_TOL.txt') )
    writeScalarToDisk(params_simu.INNER_MAXITER, os.path.join('.', tmpDirName, 'iterative_data/INNER_MAXITER.txt') )
    writeScalarToDisk(params_simu.INNER_RESTART, os.path.join('.', tmpDirName, 'iterative_data/INNER_RESTART.txt') )
    writeScalarToDisk(N_RWG, os.path.join('.', tmpDirName, 'ZI/ZI_size.txt') )
    MPI.COMM_WORLD.Barrier()
    computeTreeParameters(my_id, tmpDirName, a, k, N_levels, params_simu)
    if params_simu.COMPUTE_Z_NEAR==1:
        chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = distributeZcubesAttributions(params_simu.MAX_BLOCK_SIZE, target_mesh.N_nearPerCube, C, tmpDirName, params_simu)
        scatterMesh(target_mesh, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, tmpDirName, my_id, num_proc)
    del target_mesh # we gain A LOT of memory here, especially if N_RWG~1e6
    # we now dump-pickle the necessary variables
    variables = {}
    variables['a'] = a
    variables['k'] = k
    variables['w'] = w
    variables['C'] = C
    variables['N_RWG'] = N_RWG
    variables['N_levels'] = N_levels
    variables['CFIE'] = CFIE
    variables['chunkNumber_to_cubesNumbers'] = chunkNumber_to_cubesNumbers
    variables['cubeNumber_to_chunkNumber'] = cubeNumber_to_chunkNumber
    variables['chunkNumber_to_processNumber'] = chunkNumber_to_processNumber
    variables['processNumber_to_ChunksNumbers'] = processNumber_to_ChunksNumbers
    file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

def compute_Z_near(params_simu):
    # computation of near field elements
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    if params_simu.COMPUTE_Z_NEAR==1:
        file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'r')
        variables = cPickle.load(file)
        file.close()
        Wall_time_Z_near_computation, CPU_time_Z_near_computation = computeZ_near(tmpDirName, variables['a'], variables['C'], variables['k'], variables['w'], variables['N_levels'], variables['chunkNumber_to_cubesNumbers'], variables['processNumber_to_ChunksNumbers'], variables['CFIE'], params_simu)
        pathToReadFrom = os.path.join('.', tmpDirName, 'Z_tmp')
        # we exchange the missing Z_near parts for each process
        Mg_listsOfZnearBlocks_ToTransmitAndReceive(variables['chunkNumber_to_cubesNumbers'], variables['cubeNumber_to_chunkNumber'], variables['chunkNumber_to_processNumber'], variables['processNumber_to_ChunksNumbers'], pathToReadFrom, 'F')
        if (my_id == 0):
            createMPIsystemCommand("./code/MoM", "communicateZnearBlocks", num_proc)
    else:
        Wall_time_Z_near_computation, CPU_time_Z_near_computation = 0.0, 0.0
    # we now dump-pickle the necessary variables
    variables['Wall_time_Z_near_computation'] = Wall_time_Z_near_computation
    variables['CPU_time_Z_near_computation'] = CPU_time_Z_near_computation
    file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

def compute_SAI(params_simu):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    if params_simu.COMPUTE_Z_NEAR==1:
        file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'r')
        variables = cPickle.load(file)
        file.close()
        Wall_time_Mg_computation, CPU_time_Mg_computation = compute_SAIpreconditioner(tmpDirName, variables['a'], variables['C'], variables['chunkNumber_to_cubesNumbers'], variables['cubeNumber_to_chunkNumber'], variables['chunkNumber_to_processNumber'], variables['processNumber_to_ChunksNumbers'], params_simu.MAX_BLOCK_SIZE)
    else:
        Wall_time_Mg_computation, CPU_time_Mg_computation = 0.0, 0.0
    variables['Wall_time_Mg_computation'] = Wall_time_Mg_computation
    variables['CPU_time_Mg_computation'] = CPU_time_Mg_computation
    if (my_id == 0) and (params_simu.VERBOSE == 1):
        print variables['CPU_time_Z_near_computation'], "CPU time (seconds) for constructing Z_CFIE_near"
        print variables['Wall_time_Z_near_computation'], "Wall time (seconds) for constructing Z_CFIE_near"
        print variables['CPU_time_Mg_computation'], "CPU time (seconds) for constructing SAI precond"
        print variables['Wall_time_Mg_computation'], "Wall time (seconds) for constructing SAI precond"
    file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()
    if (my_id == 0):
        createMPIsystemCommand("./code/MoM", "mpi_mlfma", num_proc)

def print_times(params_simu):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    if (my_id == 0):
        # final moves
        file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'r')
        variables = cPickle.load(file)
        file.close()
        numberOfMatvecs = readIntFromDisk(os.path.join('.',tmpDirName,'iterative_data/numberOfMatvecs.txt'))
        NIterMLFMA = readIntFromDisk(os.path.join('.',tmpDirName,'iterative_data/iter.txt'))
        average_RWG_length = readFloatFromDisk(os.path.join('.',tmpDirName,'mesh/average_RWG_length.txt'))
        file = open('./result/CPU_time_MLFMA.txt', 'r')
        CPU_time_MLFMA_tmp = file.readlines()
        file.close()
        CPU_time_MLFMA = 0.0
        for line in CPU_time_MLFMA_tmp:
            if 'real' in line:
                CPU_time_MLFMA = float(line.split()[1])
        if (params_simu.VERBOSE == 1):
            print "N RWG =", variables['N_RWG']
            sys.stdout.write("average RWG length = %.5s" %str(average_RWG_length) + " m = lambda / %.9s" %str((c/params_simu.f)/average_RWG_length) + " \n")
            print variables['CPU_time_Z_near_computation'], "CPU time (seconds) for constructing Z_CFIE_near"
            print variables['Wall_time_Z_near_computation'], "Wall time (seconds) for constructing Z_CFIE_near"
            print variables['CPU_time_Mg_computation'], "CPU time (seconds) for constructing SAI precond"
            print variables['Wall_time_Mg_computation'], "Wall time (seconds) for constructing SAI precond"
            print CPU_time_MLFMA, "CPU time (seconds) for MLFMA iterations and solution. Iterations =", NIterMLFMA
            #print target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for MLFMA iterations and solution. Iterations =", target_MLFMA.NIterMLFMA
            if numberOfMatvecs>0:
                print CPU_time_MLFMA/numberOfMatvecs, "CPU time (seconds) per MLFMA matvec"
                #print target_MLFMA.Wall_time_Target_MLFMA_resolution/target_MLFMA.numberOfMatvecs, "Wall time (seconds) per MLFMA matvec"
            print variables['CPU_time_Z_near_computation'] + variables['CPU_time_Mg_computation'] + CPU_time_MLFMA, "CPU time (seconds) for complete MLFMA solution"
            #print Wall_time_Z_near_computation + Wall_time_Mg_computation + target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for complete MLFMA solution"
        if params_simu.CURRENTS_VISUALIZATION:
            computeCurrentsVisualization(params_simu, variables)
        if params_simu.SAVE_CURRENTS_CENTROIDS:
            saveCurrentsCentroids(params_simu)

def computeCurrentsVisualization(params_simu, variables):
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    if (my_id == 0):
        target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join('.', tmpDirName, "mesh"))
        CURRENTS_VISUALIZATION = params_simu.CURRENTS_VISUALIZATION and (target_mesh.N_RWG<2e6) and (params_simu.BISTATIC == 1)
        if CURRENTS_VISUALIZATION:
            print "............Constructing visualisation of currents"
            if (target_mesh.N_RWG<1e4):
                nbTimeSteps = 48
            else:
                nbTimeSteps = 1
            tmpDirName = 'tmp' + str(my_id)
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            norm_J_centroids_triangles = normJMCentroidsTriangles(J_centroids_triangles, variables['w'], nbTimeSteps)
            write_VectorFieldTrianglesCentroidsGMSH(os.path.join(params_simu.pathToTarget, params_simu.targetName) + '.J_centroids_triangles.pos', real(J_centroids_triangles), target_mesh)
            write_ScalarFieldTrianglesCentroidsGMSH(os.path.join(params_simu.pathToTarget, params_simu.targetName) + '.norm_J_centroids_triangles.pos', norm_J_centroids_triangles, target_mesh)
            print "............end of currents visualisation construction."
        else:
            "Error in the computeCurrentsVisualization routine. Probably you required currents visualization computation for a mesh too big"
            sys.exit(1)

def saveCurrentsCentroids(params_simu):
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    if (my_id == 0):
        target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join('.', tmpDirName, "mesh"))
        SAVE_CURRENTS_CENTROIDS = (params_simu.SAVE_CURRENTS_CENTROIDS and (params_simu.BISTATIC == 1))
        if SAVE_CURRENTS_CENTROIDS:
            print "............saving currents at centroids of triangles"
            nbTimeSteps = 1
            tmpDirName = 'tmp' + str(my_id)
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            write_VectorFieldTrianglesCentroids(os.path.join('./result', 'J_centroids_triangles.txt'), J_centroids_triangles, target_mesh)
            print "............end of saving of currents."

def getField(filename):
    my_id = MPI.COMM_WORLD.Get_rank()
    E_obs = zeros(3, 'F')
    if (my_id == 0):
        E_obs = readASCIIBlitzComplexFloatArray2DFromDisk(filename)
    return E_obs

def communicateParamsFile(filename):
    my_id = MPI.COMM_WORLD.Get_rank()
    num_proc = MPI.COMM_WORLD.Get_size()
    # we communicate the simus params files   
    if my_id==0:
        f = open(filename,'r')
        params = f.readlines()
        f.close()
    else:
        params = ['blabla']
    params = MPI.COMM_WORLD.Bcast(params)
    for i in range(num_proc):
        if (my_id!=0) and (my_id==i):
            f = open(filename, 'w')
            for line in params:
                f.write(line)
            f.close()
    MPI.COMM_WORLD.Barrier()

