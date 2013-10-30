import sys, os, argparse, time, pickle, cPickle, commands
from mpi4py import MPI
from scipy import array, sqrt, take, sum, mean, zeros
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk, writeBlitzArrayToDisk, read1DBlitzArrayFromDisk
from ReadWriteBlitzArray import readIntFromDisk, readASCIIBlitzIntArray1DFromDisk, readBlitzArrayFromDisk
from MLFMA import computeTreeParameters
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2, read_mesh_GiD, read_mesh_ANSYS
from mesh_functions_seb import compute_RWG_meanEdgeLength

def readMesh(path, name, params_simu):
    if params_simu.meshFormat == 'GMSH':
        #vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    elif params_simu.meshFormat == 'GiD':
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GiD(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    elif params_simu.meshFormat == 'ANSYS':
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_ANSYS(path, name, params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    else:
        print "setup_MLFMA_mesh.py : error on the mesh format. Enter a correct one please."
    return vertexes_coord, triangle_vertexes, triangles_physicalSurface

def compute_RWG_CFIE_OK(triangles_surfaces, RWGNumber_signedTriangles, IS_CLOSED_SURFACE):
    RWGNumber_CFIE_OK_tmp1 = take(triangles_surfaces, RWGNumber_signedTriangles, axis=0)
    RWGNumber_CFIE_OK_tmp2 = take(IS_CLOSED_SURFACE, RWGNumber_CFIE_OK_tmp1, axis=0)
    # following the Taskinen et al. paper in PIER
    # we can have a CFIE on a junction straddling a dielectric and metallic body
    RWGNumber_CFIE_OK = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1)>=1) * 1).astype('i')
    # We cannot have M on a junction between a dielectric and metallic body!
    # The following expression is lacking the fact that a surface can be metallic
    # or dielectric. If metallic, then there is no M current, even if surface is closed
    RWGNumber_M_CURRENT_OK = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1)>1) * 1).astype('i')
    return RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK

def setup_mesh(params_simu, simuDirName):
    """Sets up the mesh.
       params_simu is a class instance that contains the parameters for the simulation.
    """
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    # size of cube at finest level
    a = c/params_simu.f * params_simu.a_factor
    if (my_id==0):
        # reading the mesh
        path, name = geoDirName, params_simu.targetName + params_simu.meshFileTermination
        meshPath = os.path.join(tmpDirName, "mesh")
        print "  reading", os.path.join(path, name), "...",
        t0 = time.clock()
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = readMesh(path, name, params_simu)
        time_reading = time.clock()-t0
        print "reading mesh time =", time_reading, "seconds"
        T = triangle_vertexes.shape[0]
        V = vertexes_coord.shape[0]
        print "  number of triangles =", T
        print "  edges classification..."
        sys.stdout.flush()

        writeScalarToDisk(T, os.path.join(meshPath, "T.txt"))
        writeScalarToDisk(V, os.path.join(meshPath, "V.txt"))
        writeScalarToDisk(a, os.path.join(meshPath, "a.txt"))
        writeBlitzArrayToDisk(vertexes_coord, os.path.join(meshPath, 'vertexes_coord.txt'))
        writeBlitzArrayToDisk(triangle_vertexes, os.path.join(meshPath, 'triangle_vertexes.txt'))

        t10 = time.clock()
        print commands.getoutput("./code/MoM/mesh_functions_seb " + meshPath + "/")
        print "time C++ execution =", time.clock() - t10
        print "edgeNumber_triangles construction cumulated time =", time.clock() - t10

        is_closed_surface = readASCIIBlitzIntArray1DFromDisk(os.path.join(meshPath, "is_closed_surface.txt"))
        S = len(is_closed_surface)
        print "test of the closed surfaces :", is_closed_surface

        N_RWG = readIntFromDisk(os.path.join(meshPath, "N_RWG.txt"))
        print "Number of RWG =", N_RWG
        sys.stdout.flush()

        triangles_surfaces = readASCIIBlitzIntArray1DFromDisk(os.path.join(meshPath, "triangles_surfaces.txt"))
        RWGNumber_signedTriangles = readBlitzArrayFromDisk(os.path.join(meshPath, "RWGNumber_signedTriangles.txt"), N_RWG, 2, 'i')
        RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK = compute_RWG_CFIE_OK(triangles_surfaces, RWGNumber_signedTriangles, is_closed_surface)
        if N_RWG<1e4:
            stride = 1
        else:
            stride = N_RWG/100
        RWGNumber_edgeVertexes = readBlitzArrayFromDisk(os.path.join(meshPath, "RWGNumber_edgeVertexes.txt"), N_RWG, 2, 'i')
        average_RWG_length = compute_RWG_meanEdgeLength(vertexes_coord, RWGNumber_edgeVertexes, stride)
        writeScalarToDisk(average_RWG_length, os.path.join(meshPath,'average_RWG_length.txt'))
        if params_simu.VERBOSE==1:
            print "average RWG length =", average_RWG_length, "m = lambda /", (c/params_simu.f)/average_RWG_length

        # cubes computation
        WEAVE = 0
        if WEAVE != 0:
            print "Using good old weave!"
            from Cubes import cube_lower_coord_computation, RWGNumber_cubeNumber_computation, cubeIndex_RWGNumbers_computation, findCubeNeighbors
            from mesh_functions_seb import compute_RWGNumber_edgeCentroidCoord
            max_N_cubes_1D, N_levels, big_cube_lower_coord, big_cube_center_coord = cube_lower_coord_computation(a, vertexes_coord)
            N_levels = max(N_levels, 2)
            RWGNumber_edgeCentroidCoord = compute_RWGNumber_edgeCentroidCoord(vertexes_coord, RWGNumber_edgeVertexes)
            RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord = RWGNumber_cubeNumber_computation(a, max_N_cubes_1D, big_cube_lower_coord, RWGNumber_edgeCentroidCoord)
            cubes_RWGsNumbers, cubes_lists_RWGsNumbers, cube_N_RWGs, cubes_centroids = cubeIndex_RWGNumbers_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord)
            print "Average number of RWGs per cube:", mean(cube_N_RWGs)
            C = cubes_centroids.shape[0]
            cubes_lists_NeighborsIndexes, cubes_neighborsIndexes, cube_N_neighbors = findCubeNeighbors(max_N_cubes_1D, big_cube_lower_coord, cubes_centroids, a)
            writeScalarToDisk(C, os.path.join(meshPath, "C.txt"))
            writeBlitzArrayToDisk(cubes_centroids, os.path.join(meshPath, 'cubes_centroids') + '.txt')
            writeBlitzArrayToDisk(cubes_RWGsNumbers, os.path.join(meshPath, 'cubes_RWGsNumbers') + '.txt')
            writeBlitzArrayToDisk(cube_N_RWGs, os.path.join(meshPath, 'cube_N_RWGs') + '.txt')
            writeBlitzArrayToDisk(cubes_neighborsIndexes, os.path.join(meshPath, 'cubes_neighborsIndexes') + '.txt')
            writeBlitzArrayToDisk(cube_N_neighbors, os.path.join(meshPath, 'cube_N_neighbors') + '.txt')
            writeScalarToDisk(N_levels, os.path.join(meshPath, "N_levels.txt"))
        else:
            print "Using new mesh_cubes.cpp code"
            print commands.getoutput("./code/MoM/mesh_cubes " + meshPath + "/")
            N_levels = readIntFromDisk(os.path.join(meshPath,'N_levels.txt'))
            max_N_cubes_1D = readIntFromDisk(os.path.join(meshPath,'max_N_cubes_1D.txt'))
            C = readIntFromDisk(os.path.join(meshPath,'C.txt'))
            big_cube_center_coord = read1DBlitzArrayFromDisk(os.path.join(meshPath, "big_cube_center_coord.txt"), 'd')
            big_cube_lower_coord = read1DBlitzArrayFromDisk(os.path.join(meshPath, "big_cube_lower_coord.txt"), 'd')
            # making of cubes_lists_RWGsNumbers
            cubes_RWGsNumbers = read1DBlitzArrayFromDisk(os.path.join(meshPath, "cubes_RWGsNumbers.txt"), 'i')
            cube_N_RWGs = read1DBlitzArrayFromDisk(os.path.join(meshPath, "cube_N_RWGs.txt"), 'i')
            print "Average number of RWGs per cube:", mean(cube_N_RWGs)
            cubes_lists_RWGsNumbers = {}
            index = 0
            for i in range(C):
                array_tmp = zeros(cube_N_RWGs[i], 'i')
                for j in range(cube_N_RWGs[i]):
                    array_tmp[j] = cubes_RWGsNumbers[index]
                    index += 1
                cubes_lists_RWGsNumbers[i] = array_tmp
            # making of cubes_lists_NeighborsIndexes
            cubesNeighborsIndexesTmp2 = readBlitzArrayFromDisk(os.path.join(meshPath, "cubesNeighborsIndexesTmp2.txt"), C, 28, 'i')
            cubes_lists_NeighborsIndexes = {}
            for i in range(C):
                cubes_lists_NeighborsIndexes[i] = [elem for elem in cubesNeighborsIndexesTmp2[i] if elem>-1]

        # writing some data
        print "N_levels = ", N_levels
        print "max_N_cubes_1D = ", max_N_cubes_1D
        print "big_cube_center_coord = ", big_cube_center_coord
        print "big_cube_lower_coord = ", big_cube_lower_coord
        file = open(os.path.join(meshPath, 'cubes_lists_RWGsNumbers.txt'), 'w')
        cPickle.dump(cubes_lists_RWGsNumbers, file)
        file.close()
        file = open(os.path.join(meshPath, 'cubes_lists_NeighborsIndexes.txt'), 'w')
        cPickle.dump(cubes_lists_NeighborsIndexes, file)
        file.close()
        writeScalarToDisk(S, os.path.join(meshPath, "S.txt"))
        writeBlitzArrayToDisk(RWGNumber_CFIE_OK, os.path.join(meshPath, 'RWGNumber_CFIE_OK') + '.txt')
        writeBlitzArrayToDisk(RWGNumber_M_CURRENT_OK, os.path.join(meshPath, 'RWGNumber_M_CURRENT_OK') + '.txt')

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
    parser.add_argument('--simuparams')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    simuParams = cmdline.simuparams
    if simuDirName==None:
        simuDirName = '.'
    if simuParams==None:
        simuParams = 'simulation_parameters'

    # the simulation itself
    sys.path.append(os.path.abspath('.'))
    exec 'from ' + simuParams + ' import *'
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_mesh(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()

