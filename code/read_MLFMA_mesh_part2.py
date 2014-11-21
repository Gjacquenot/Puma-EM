import sys, os, argparse, time
try:
    import commands
except ImportError:
    import subprocess as commands
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from scipy import take, sum, mean, zeros
from ReadWriteBlitzArray import writeScalarToDisk, writeBlitzArrayToDisk, read1DBlitzArrayFromDisk
from ReadWriteBlitzArray import readIntFromDisk, readFloatFromDisk, readASCIIBlitzIntArray1DFromDisk, readBlitzArrayFromDisk

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
    my_id = 0

    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    meshPath = os.path.join(tmpDirName, "mesh")
    # size of cube at finest level
    a = c/params_simu.f * params_simu.a_factor

    is_closed_surface = readASCIIBlitzIntArray1DFromDisk(os.path.join(meshPath, "is_closed_surface.txt"))
    S = len(is_closed_surface)
    print("test of the closed surfaces : " + str(is_closed_surface))

    N_RWG = readIntFromDisk(os.path.join(meshPath, "N_RWG.txt"))
    print("Number of RWG = " + str(N_RWG))
    sys.stdout.flush()

    triangles_surfaces = readASCIIBlitzIntArray1DFromDisk(os.path.join(meshPath, "triangles_surfaces.txt"))
    RWGNumber_signedTriangles = readBlitzArrayFromDisk(os.path.join(meshPath, "RWGNumber_signedTriangles.txt"), N_RWG, 2, 'i')
    RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK = compute_RWG_CFIE_OK(triangles_surfaces, RWGNumber_signedTriangles, is_closed_surface)
    RWGNumber_edgeVertexes = readBlitzArrayFromDisk(os.path.join(meshPath, "RWGNumber_edgeVertexes.txt"), N_RWG, 2, 'i')
    average_RWG_length = readFloatFromDisk(os.path.join(meshPath,'average_RWG_length.txt'))
    if params_simu.VERBOSE==1:
        print("average RWG length = " + str(average_RWG_length) + "m = lambda /" + str((c/params_simu.f)/average_RWG_length))

    # cubes computation
    WEAVE = 0
    if WEAVE != 0:
        print("Using good old weave!")
        from Cubes import cube_lower_coord_computation, RWGNumber_cubeNumber_computation, cubeIndex_RWGNumbers_computation, findCubeNeighbors
        from mesh_functions_seb import compute_RWGNumber_edgeCentroidCoord
        max_N_cubes_1D, N_levels, big_cube_lower_coord, big_cube_center_coord = cube_lower_coord_computation(a, vertexes_coord)
        N_levels = max(N_levels, 2)
        RWGNumber_edgeCentroidCoord = compute_RWGNumber_edgeCentroidCoord(vertexes_coord, RWGNumber_edgeVertexes)
        RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord = RWGNumber_cubeNumber_computation(a, max_N_cubes_1D, big_cube_lower_coord, RWGNumber_edgeCentroidCoord)
        cubes_RWGsNumbers, cubes_lists_RWGsNumbers, cube_N_RWGs, cubes_centroids = cubeIndex_RWGNumbers_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord)
        print("Average number of RWGs per cube: " + str(mean(cube_N_RWGs)))
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
        print("Using new mesh_cubes.cpp code")
        print(commands.getoutput("./code/MoM/mesh_cubes " + meshPath + "/"))
        C = readIntFromDisk(os.path.join(meshPath,'C.txt'))
        # making of cubes_lists_RWGsNumbers
        cubes_RWGsNumbers = read1DBlitzArrayFromDisk(os.path.join(meshPath, "cubes_RWGsNumbers.txt"), 'i')
        cube_N_RWGs = read1DBlitzArrayFromDisk(os.path.join(meshPath, "cube_N_RWGs.txt"), 'i')
        print("Average number of RWGs per cube: " + str(mean(cube_N_RWGs)))
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
    file = open(os.path.join(meshPath, 'cubes_lists_RWGsNumbers.txt'), 'wb')
    cPickle.dump(cubes_lists_RWGsNumbers, file)
    file.close()
    file = open(os.path.join(meshPath, 'cubes_lists_NeighborsIndexes.txt'), 'wb')
    cPickle.dump(cubes_lists_NeighborsIndexes, file)
    file.close()
    writeScalarToDisk(S, os.path.join(meshPath, "S.txt"))
    writeBlitzArrayToDisk(RWGNumber_CFIE_OK, os.path.join(meshPath, 'RWGNumber_CFIE_OK') + '.txt')
    writeBlitzArrayToDisk(RWGNumber_M_CURRENT_OK, os.path.join(meshPath, 'RWGNumber_M_CURRENT_OK') + '.txt')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--inputdir')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    inputDirName = cmdline.inputdir
    simuParams = 'simulation_parameters'

    # the simulation itself
    sys.path.append(os.path.abspath(inputDirName))
    exec('from ' + simuParams + ' import *')
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_mesh(params_simu, simuDirName)
    else:
        print("you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings.")
        sys.exit(1)

