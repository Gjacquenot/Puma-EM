import sys, os, argparse, time
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from mpi4py import MPI
from scipy import array, sqrt, pi
from ReadWriteBlitzArray import readIntFromDisk, writeScalarToDisk, writeASCIIBlitzArrayToDisk, read1DBlitzArrayFromDisk

def setup_mesh(params_simu, simuDirName):
    """Sets up the mesh.
       params_simu is a class instance that contains the parameters for the simulation.
    """
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    meshPath = os.path.join(tmpDirName, "mesh")

    # size of cube at finest level
    a = params_simu.c/params_simu.f * params_simu.a_factor
    if (my_id==0):
        N_RWG = readIntFromDisk(os.path.join(meshPath, "N_RWG.txt"))
        N_levels = readIntFromDisk(os.path.join(meshPath,'N_levels.txt'))
        max_N_cubes_1D = readIntFromDisk(os.path.join(meshPath,'max_N_cubes_1D.txt'))
        C = readIntFromDisk(os.path.join(meshPath,'C.txt'))
        big_cube_center_coord = read1DBlitzArrayFromDisk(os.path.join(meshPath, "big_cube_center_coord.txt"), 'd')
        big_cube_lower_coord = read1DBlitzArrayFromDisk(os.path.join(meshPath, "big_cube_lower_coord.txt"), 'd')

        # writing some data
        print("N_levels = " + str(N_levels))
        print("max_N_cubes_1D = " + str(max_N_cubes_1D))
        print("big_cube_center_coord = " + str(big_cube_center_coord))
        print("big_cube_lower_coord = " + str(big_cube_lower_coord))

    else:
        big_cube_lower_coord = ['blabla']
        big_cube_center_coord = ['blabla']
        N_levels = ['blabla']
        N_RWG = ['blabla']
        C = ['blabla']
    big_cube_lower_coord = MPI.COMM_WORLD.bcast(big_cube_lower_coord)
    big_cube_center_coord = MPI.COMM_WORLD.bcast(big_cube_center_coord)
    N_levels = MPI.COMM_WORLD.bcast(N_levels)
    N_RWG = MPI.COMM_WORLD.bcast(N_RWG)
    C = MPI.COMM_WORLD.bcast(C)

    w = 2. * pi * params_simu.f
    k = w * sqrt(params_simu.eps_0*params_simu.eps_r*params_simu.mu_0*params_simu.mu_r) + 1.j * 0.
    CFIE = array(params_simu.CFIE).astype('D')

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
    writeScalarToDisk(params_simu.TDS_APPROX*1, os.path.join(tmpDirName, 'octtree_data/TDS_APPROX.txt') )
    writeScalarToDisk(params_simu.Z_s, os.path.join(tmpDirName, 'octtree_data/Z_s.txt') )
    # what type of simulation are we running?
    writeScalarToDisk(params_simu.BISTATIC*1, os.path.join(tmpDirName, 'BISTATIC.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_RCS*1, os.path.join(tmpDirName, 'MONOSTATIC_RCS.txt') )
    writeScalarToDisk(params_simu.MONOSTATIC_SAR*1, os.path.join(tmpDirName, 'MONOSTATIC_SAR.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HH*1, os.path.join(tmpDirName, 'COMPUTE_RCS_HH.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_VV*1, os.path.join(tmpDirName, 'COMPUTE_RCS_VV.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_HV*1, os.path.join(tmpDirName, 'COMPUTE_RCS_HV.txt') )
    writeScalarToDisk(params_simu.COMPUTE_RCS_VH*1, os.path.join(tmpDirName, 'COMPUTE_RCS_VH.txt') )
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

    variables = {}
    variables['a'] = a
    variables['k'] = k
    variables['w'] = w
    variables['C'] = C
    variables['N_RWG'] = N_RWG
    variables['N_levels'] = N_levels
    variables['CFIE'] = CFIE
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'wb')
    cPickle.dump(variables, file)
    file.close()    


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
        params_simu.c = c
        params_simu.eps_0 = eps_0
        params_simu.mu_0 = mu_0
    else:
        params_simu = ['blabla']
    params_simu = MPI.COMM_WORLD.bcast(params_simu)
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_mesh(params_simu, simuDirName)
    else:
        print("you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings.")
        sys.exit(1)
    #MPI.Finalize()

