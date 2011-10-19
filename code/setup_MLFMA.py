import sys, os, argparse
import time, copy, commands, pickle, cPickle
from mpi4py import *
from scipy import zeros, array, sqrt
from meshClass import MeshClass, CubeClass
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk
from MLFMA import computeTreeParameters

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
    params = MPI.COMM_WORLD.bcast(params)
    for i in range(num_proc):
        if (my_id!=0) and (my_id==i):
            f = open(filename, 'w')
            for line in params:
                f.write(line)
            f.close()
    MPI.COMM_WORLD.Barrier()


def setup_excitation(params_simu, tmpDirName):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    # the observation points
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
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join(tmpDirName,'V_CFIE/r_obs.txt'))
    elif (params_simu.r_obs_FROM_FILE == 1) and (params_simu.r_obs_FILENAME != ""):
        if (my_id==0): # this file is only on processor 0
            r_obs = read_observation_points(params_simu.r_obs_FILENAME)
        else:
            r_obs = zeros((1, 3), 'd')
        r_obs = MPI.COMM_WORLD.bcast(r_obs)
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join(tmpDirName,'V_CFIE/r_obs.txt'))

    # now the excitations
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_DIPOLES, os.path.join(tmpDirName,'V_CFIE/DIPOLES_EXCITATION.txt'))
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_PLANE_WAVE, os.path.join(tmpDirName,'V_CFIE/PLANE_WAVE_EXCITATION.txt'))
    writeScalarToDisk(params_simu.V_FULL_PRECISION*1, os.path.join(tmpDirName, 'V_CFIE/V_FULL_PRECISION.txt') )
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
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join(tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join(tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
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
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join(tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join(tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # if we have dipoles excitation AND definition of the excitation in a user-supplied file
    elif (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_DIPOLES_FROM_FILE == 1):
        if params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                J_src, r_J_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME)
            else:
                J_src, r_J_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            J_src = MPI.COMM_WORLD.bcast(J_src)
            r_J_src = MPI.COMM_WORLD.bcast(r_J_src)
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join(tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join(tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
        if params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                M_src, r_M_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME)
            else:
                M_src, r_M_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            M_src = MPI.COMM_WORLD.bcast(M_src)
            r_M_src = MPI.COMM_WORLD.bcast(r_M_src)
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join(tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join(tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # now the plane wave excitation
    if params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1:
        writeScalarToDisk(params_simu.theta_inc, os.path.join(tmpDirName,'V_CFIE/theta_inc.txt'))
        writeScalarToDisk(params_simu.phi_inc, os.path.join(tmpDirName,'V_CFIE/phi_inc.txt'))
        E_inc = array([params_simu.E_inc_theta, params_simu.E_inc_phi], 'D')
        writeASCIIBlitzArrayToDisk(E_inc, os.path.join(tmpDirName,'V_CFIE/E_inc.txt'))
    if (params_simu.BISTATIC_EXCITATION_DIPOLES != 1) and (params_simu.BISTATIC_EXCITATION_PLANE_WAVE != 1):
        if (my_id==0):
            print "incorrect excitation choice. You have to choose dipole and/or plane wave excitation."
        sys.exit(1)
    if params_simu.MONOSTATIC_SAR==1:
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_x_hat, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_local_x_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_y_hat, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_local_y_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_plane_origin, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_plane_origin.txt'))
        writeScalarToDisk(params_simu.SAR_x_span, os.path.join(tmpDirName,'V_CFIE/SAR_x_span.txt'))
        writeScalarToDisk(params_simu.SAR_y_span, os.path.join(tmpDirName,'V_CFIE/SAR_y_span.txt'))
        writeScalarToDisk(params_simu.SAR_x_span_offset, os.path.join(tmpDirName,'V_CFIE/SAR_x_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_y_span_offset, os.path.join(tmpDirName,'V_CFIE/SAR_y_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_N_x_points, os.path.join(tmpDirName,'V_CFIE/SAR_N_x_points.txt'))
        writeScalarToDisk(params_simu.SAR_N_y_points, os.path.join(tmpDirName,'V_CFIE/SAR_N_y_points.txt'))


def setup_MLFMA(params_simu, tmpDirName):
    """Sets up the MLFMA parameters.
       params_simu is a class instance that contains the parameters for the simulation.
    """
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    # target_mesh construction
    target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
    # size of cube at finest level
    a = c/params_simu.f * params_simu.a_factor
    if (my_id==0):
        target_mesh.constructFromGmshFile()
        writeScalarToDisk(target_mesh.average_RWG_length, os.path.join(tmpDirName,'mesh/average_RWG_length.txt'))
        if params_simu.VERBOSE==1:
            print "average RWG length =", target_mesh.average_RWG_length, "m = lambda /", (c/params_simu.f)/target_mesh.average_RWG_length
        # target_mesh cubes computation
        target_mesh.cubes_data_computation(a)
        target_mesh.saveToDisk(os.path.join(tmpDirName,'mesh'))
        IS_CLOSED_SURFACE = target_mesh.IS_CLOSED_SURFACE
        big_cube_lower_coord = target_mesh.big_cube_lower_coord
        big_cube_center_coord = target_mesh.big_cube_center_coord
        N_levels = target_mesh.N_levels
        N_RWG = target_mesh.N_RWG
        T = target_mesh.T
        C = target_mesh.C
    else:
        IS_CLOSED_SURFACE = ['blabla']
        big_cube_lower_coord = ['blabla']
        big_cube_center_coord = ['blabla']
        N_levels = ['blabla']
        N_RWG = ['blabla']
        T = ['blabla']
        C = ['blabla']
    del target_mesh
    big_cube_lower_coord = MPI.COMM_WORLD.bcast(big_cube_lower_coord)
    big_cube_center_coord = MPI.COMM_WORLD.bcast(big_cube_center_coord)
    IS_CLOSED_SURFACE = MPI.COMM_WORLD.bcast(IS_CLOSED_SURFACE)
    N_levels = MPI.COMM_WORLD.bcast(N_levels)
    N_RWG = MPI.COMM_WORLD.bcast(N_RWG)
    C = MPI.COMM_WORLD.bcast(C)
    T = MPI.COMM_WORLD.bcast(T)

    w = 2. * pi * params_simu.f
    k = w * sqrt(eps_0*params_simu.eps_r*mu_0*params_simu.mu_r) + 1.j * 0.
    CFIE = array([params_simu.nu, 0, 0, -(1.0-params_simu.nu) * 377]).astype('D')

    writeScalarToDisk( num_procs, os.path.join(tmpDirName,'octtree_data/num_procs.txt') )
    writeScalarToDisk( a, os.path.join(tmpDirName,'octtree_data/leaf_side_length.txt') )
    writeScalarToDisk(2.0*pi*params_simu.f, os.path.join(tmpDirName,'octtree_data/w.txt') )
    writeScalarToDisk(params_simu.eps_r, os.path.join(tmpDirName,'octtree_data/eps_r.txt') )
    writeScalarToDisk(params_simu.mu_r, os.path.join(tmpDirName,'octtree_data/mu_r.txt') )
    writeScalarToDisk(k, os.path.join(tmpDirName,'octtree_data/k.txt') )
    writeASCIIBlitzArrayToDisk(CFIE, os.path.join(tmpDirName,'octtree_data/CFIEcoeffs.txt') )
    writeScalarToDisk(N_RWG, os.path.join(tmpDirName,'octtree_data/N_RWG.txt') )
    writeScalarToDisk(N_levels-1, os.path.join(tmpDirName,'octtree_data/N_active_levels.txt') )
    writeASCIIBlitzArrayToDisk(big_cube_lower_coord, os.path.join(tmpDirName,'octtree_data/big_cube_lower_coord.txt') )
    writeASCIIBlitzArrayToDisk(big_cube_center_coord, os.path.join(tmpDirName,'octtree_data/big_cube_center_coord.txt') )
    writeScalarToDisk(params_simu.PERIODIC_Theta*1, os.path.join(tmpDirName,'octtree_data/PERIODIC_Theta.txt') )
    writeScalarToDisk(params_simu.CYCLIC_Theta*1, os.path.join(tmpDirName,'octtree_data/CYCLIC_Theta.txt') )
    writeScalarToDisk(params_simu.PERIODIC_Phi*1, os.path.join(tmpDirName,'octtree_data/PERIODIC_Phi.txt') )
    writeScalarToDisk(params_simu.CYCLIC_Phi*1, os.path.join(tmpDirName,'octtree_data/CYCLIC_Phi.txt') )
    writeScalarToDisk(params_simu.ALLOW_CEILING_LEVEL*1, os.path.join(tmpDirName, 'octtree_data/ALLOW_CEILING_LEVEL.txt') )
    writeScalarToDisk(params_simu.DIRECTIONS_PARALLELIZATION*1, os.path.join(tmpDirName, 'octtree_data/DIRECTIONS_PARALLELIZATION.txt') )
    writeScalarToDisk(params_simu.BE_BH_N_Gauss_points, os.path.join(tmpDirName, 'octtree_data/N_GaussOnTriangle.txt') )
    writeScalarToDisk(params_simu.MOM_FULL_PRECISION*1, os.path.join(tmpDirName, 'octtree_data/MOM_FULL_PRECISION.txt') )
    writeScalarToDisk(params_simu.VERBOSE*1, os.path.join(tmpDirName, 'octtree_data/VERBOSE.txt') )
    # what type of simulation are we running?
    writeScalarToDisk(params_simu.BISTATIC*1, os.path.join(tmpDirName, 'BISTATIC.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_RCS*1, os.path.join(tmpDirName, 'MONOSTATIC_RCS.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_SAR*1, os.path.join(tmpDirName, 'MONOSTATIC_SAR.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HH*1, os.path.join(tmpDirName, 'COMPUTE_RCS_HH.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_VV*1, os.path.join(tmpDirName, 'COMPUTE_RCS_VV.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HV*1, os.path.join(tmpDirName, 'COMPUTE_RCS_HV.txt') )
    writeScalarToDisk(params_simu.USE_PREVIOUS_SOLUTION*1, os.path.join(tmpDirName, 'USE_PREVIOUS_SOLUTION.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_BY_BISTATIC_APPROX*1, os.path.join(tmpDirName, 'MONOSTATIC_BY_BISTATIC_APPROX.txt') )
    writeScalarToDisk(params_simu.MAXIMUM_DELTA_PHASE, os.path.join(tmpDirName, 'MAXIMUM_DELTA_PHASE.txt') )
    # writing the iterative solver setup
    restrt = min(params_simu.RESTART, N_RWG)
    writeScalarToDisk(params_simu.MAXITER, os.path.join(tmpDirName, 'iterative_data/MAXITER.txt') )
    writeScalarToDisk(restrt, os.path.join(tmpDirName, 'iterative_data/RESTART.txt') )
    writeScalarToDisk(params_simu.SOLVER, os.path.join(tmpDirName, 'iterative_data/SOLVER.txt') )
    writeScalarToDisk(params_simu.INNER_SOLVER, os.path.join(tmpDirName, 'iterative_data/INNER_SOLVER.txt') )
    writeScalarToDisk(params_simu.TOL, os.path.join(tmpDirName, 'iterative_data/TOL.txt') )
    writeScalarToDisk(params_simu.INNER_TOL, os.path.join(tmpDirName, 'iterative_data/INNER_TOL.txt') )
    writeScalarToDisk(params_simu.INNER_MAXITER, os.path.join(tmpDirName, 'iterative_data/INNER_MAXITER.txt') )
    writeScalarToDisk(params_simu.INNER_RESTART, os.path.join(tmpDirName, 'iterative_data/INNER_RESTART.txt') )
    writeScalarToDisk(N_RWG, os.path.join(tmpDirName, 'ZI/ZI_size.txt') )
    computeTreeParameters(my_id, tmpDirName, a, k, N_levels, params_simu)

    variables = {}
    variables['a'] = a
    variables['k'] = k
    variables['w'] = w
    variables['C'] = C
    variables['N_RWG'] = N_RWG
    variables['N_levels'] = N_levels
    variables['CFIE'] = CFIE
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()    


if __name__=='__main__':
    my_id = MPI.COMM_WORLD.Get_rank()
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    if simuDirName==None:
        simuDirName = '.'
    if (my_id==0):
        if 'result' not in os.listdir(simuDirName):
            os.mkdir(os.path.join(simuDirName, 'result'))
    # MLFMA parameters
    # we communicate to the others processors the params files
    communicateParamsFile("MLFMA_parameters.py")
    communicateParamsFile("simulation_parameters.py")
    # creation of the directories
    sys.path.append(os.path.abspath('.'))
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    os.mkdir( tmpDirName )
    os.mkdir( os.path.join(tmpDirName,'Z_tmp') )
    os.mkdir( os.path.join(tmpDirName,'Z_near') )
    os.mkdir( os.path.join(tmpDirName,'Mg_LeftFrob') )
    os.mkdir( os.path.join(tmpDirName,'mesh') )
    os.mkdir( os.path.join(tmpDirName,'octtree_data') )
    os.mkdir( os.path.join(tmpDirName,'V_CFIE') )
    os.mkdir( os.path.join(tmpDirName,'ZI') )
    os.mkdir( os.path.join(tmpDirName,'iterative_data') )
    os.mkdir( os.path.join(tmpDirName,'pickle') )
    # the simulation itself
    from simulation_parameters import *
    if (my_id==0):
        params_simu.display()
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_excitation(params_simu, tmpDirName)
        setup_MLFMA(params_simu, tmpDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    MPI.Finalize()

