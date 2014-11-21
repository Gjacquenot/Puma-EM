import sys, os, argparse, time
from ReadWriteBlitzArray import writeScalarToDisk, writeBlitzArrayToDisk
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2, read_mesh_GiD, read_mesh_ANSYS

def readMesh(path, name, params_simu):
    if params_simu.meshFormat == 'GMSH':
        #vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    elif params_simu.meshFormat == 'GiD':
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GiD(os.path.join(path, name), params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    elif params_simu.meshFormat == 'ANSYS':
        vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_ANSYS(path, name, params_simu.targetDimensions_scaling_factor, params_simu.z_offset)
    else:
        print("setup_MLFMA_mesh.py : error on the mesh format. Enter a correct one please.")
    return vertexes_coord, triangle_vertexes, triangles_physicalSurface

def setup_mesh(params_simu, simuDirName):
    """Sets up the mesh.
       params_simu is a class instance that contains the parameters for the simulation.
    """
    my_id = 0
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    # size of cube at finest level
    a = c/params_simu.f * params_simu.a_factor

    # reading the mesh
    path, name = geoDirName, params_simu.targetName + params_simu.meshFileTermination
    meshPath = os.path.join(tmpDirName, "mesh")
    print("read_MLFMA_mesh_part1.py: reading" + os.path.join(path, name) + "...")
    t0 = time.clock()
    vertexes_coord, triangle_vertexes, triangles_physicalSurface = readMesh(path, name, params_simu)
    time_reading = time.clock()-t0
    print("reading mesh time = " + str(time_reading) + " seconds")
    T = triangle_vertexes.shape[0]
    V = vertexes_coord.shape[0]
    print("number of triangles = " + str(T) + "\n")
    sys.stdout.flush()

    writeScalarToDisk(T, os.path.join(meshPath, "T.txt"))
    writeScalarToDisk(V, os.path.join(meshPath, "V.txt"))
    writeScalarToDisk(a, os.path.join(meshPath, "a.txt"))
    writeBlitzArrayToDisk(vertexes_coord, os.path.join(meshPath, 'vertexes_coord.txt'))
    writeBlitzArrayToDisk(triangle_vertexes, os.path.join(meshPath, 'triangle_vertexes.txt'))

	
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
