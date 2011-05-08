import sys, os
import time, copy, commands, pickle, cPickle
from mpi4py import *
from meshClass import MeshClass, CubeClass
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk
from FMM_Znear import Z_near_size_computation, Z_nearChunksDistribution

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

def distribute_Chunks_and_Mesh(params_simu):
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = 'tmp' + str(my_id)
    
    # creating the mesh
    target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
    if my_id==0:
        target_mesh.constructFromSavedArrays(os.path.join('.', tmpDirName, "mesh"))
        N_nearBlockDiag, N_near, N_nearPerCube = Z_near_size_computation(target_mesh.cubes_lists_RWGsNumbers, target_mesh.cubesNeighborsIndexes)
    else:
        N_nearPerCube = ['blabla']
    N_nearPerCube = MPI.COMM_WORLD.Bcast(N_nearPerCube)
    file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = Z_nearChunksDistribution(params_simu.MAX_BLOCK_SIZE, N_nearPerCube, variables['C'], tmpDirName)
    chunkNumber_to_cubesNumbers = MPI.COMM_WORLD.Bcast(chunkNumber_to_cubesNumbers)
    cubeNumber_to_chunkNumber = MPI.COMM_WORLD.Bcast(cubeNumber_to_chunkNumber)
    chunkNumber_to_processNumber = MPI.COMM_WORLD.Bcast(chunkNumber_to_processNumber)
    processNumber_to_ChunksNumbers = MPI.COMM_WORLD.Bcast(processNumber_to_ChunksNumbers)
    # distributing chunks of the mesh
    scatterMesh(target_mesh, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, tmpDirName, my_id, num_procs)
    del target_mesh
    # we now dump-pickle the necessary variables
    variables['chunkNumber_to_cubesNumbers'] = chunkNumber_to_cubesNumbers
    variables['cubeNumber_to_chunkNumber'] = cubeNumber_to_chunkNumber
    variables['chunkNumber_to_processNumber'] = chunkNumber_to_processNumber
    variables['processNumber_to_ChunksNumbers'] = processNumber_to_ChunksNumbers
    file = open(os.path.join('.', tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

if __name__=='__main__':
    sys.path.append(os.path.abspath('.'))
    from simulation_parameters import *
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        distribute_Chunks_and_Mesh(params_simu)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    MPI.Finalize()
