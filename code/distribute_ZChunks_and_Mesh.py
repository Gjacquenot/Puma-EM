import sys, os, argparse
import time, copy, commands, pickle, cPickle
from mpi4py import MPI
from meshClass import MeshClass, CubeClass
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk
from FMM_Znear import Z_near_size_computation, Z_nearChunksDistribution

def scatterMesh(target_mesh, ZprocessNumber_to_ChunksNumbers, ZchunkNumber_to_cubesNumbers, tmpDirName, my_id, num_proc):
    # we now exchange the local (cubes) meshes data...
    CPU_time, Wall_time = time.clock(), time.time()
    pathToSaveTo = os.path.join(tmpDirName, 'Z_tmp')
    if (my_id == 0 ):
        print "Exchanging local meshes-cubes data..."
    # I create the necessary directories for Z_tmp
    if (my_id == 0): # master node
        for recv_id in range(num_proc-1, -1, -1):
            for chunkNumber in ZprocessNumber_to_ChunksNumbers[recv_id]:
                list_cubeIntArrays, list_cubeDoubleArrays = [], []
                for cubeNumber in ZchunkNumber_to_cubesNumbers[chunkNumber]:
                    # for each one we compute the necessary information for individual Zcube computation
                    #cubeIntArrays, cubeDoubleArrays = target_mesh.computeCubeLocalArrays(cubeNumber)
                    cubeIntArrays, cubeDoubleArrays = target_mesh.computeCubeLocalArrays_C(cubeNumber)
                    # we append the cube lists to new lists
                    list_cubeIntArrays.append(cubeIntArrays)
                    list_cubeDoubleArrays.append(cubeDoubleArrays)
                # communicating the arrays
                # we exchange the concatenated arrays
                if (recv_id!=0): # master node
                    MPI.COMM_WORLD.send(list_cubeIntArrays, dest=recv_id, tag=chunkNumber)
                    MPI.COMM_WORLD.send(list_cubeDoubleArrays, dest=recv_id, tag=chunkNumber+1)
                else:
                    pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
                    os.mkdir(pathToSaveToChunk)
                    for j in range(len(ZchunkNumber_to_cubesNumbers[chunkNumber])):
                        cubeNumber = ZchunkNumber_to_cubesNumbers[chunkNumber][j]
                        cube = CubeClass()
                        cube.cubeIntArrays = list_cubeIntArrays[j]
                        cube.cubeDoubleArrays = list_cubeDoubleArrays[j]
                        cube.writeIntDoubleArraysToFile(pathToSaveToChunk, cubeNumber)
    else:
        for chunkNumber in ZprocessNumber_to_ChunksNumbers[my_id]:
            list_cubeIntArrays, list_cubeDoubleArrays = ['blabla'], ['blabla']
            list_cubeIntArrays = MPI.COMM_WORLD.recv(list_cubeIntArrays, source=0, tag=chunkNumber)
            list_cubeDoubleArrays = MPI.COMM_WORLD.recv(list_cubeDoubleArrays, source=0, tag=chunkNumber+1)
            # writing the local cube data to disk
            pathToSaveToChunk = os.path.join(pathToSaveTo, "chunk" + str(chunkNumber))
            os.mkdir(pathToSaveToChunk)
            for j in range(len(ZchunkNumber_to_cubesNumbers[chunkNumber])):
                cubeNumber = ZchunkNumber_to_cubesNumbers[chunkNumber][j]
                cube = CubeClass()
                cube.cubeIntArrays = list_cubeIntArrays[j]
                cube.cubeDoubleArrays = list_cubeDoubleArrays[j]
                cube.writeIntDoubleArraysToFile(pathToSaveToChunk, cubeNumber)

    CPU_time, Wall_time = time.clock() - CPU_time, time.time() - Wall_time
    print "Process", my_id, "mesh scattering: CPU time =", CPU_time, "sec"
    print "Process", my_id, "mesh scattering: Wall time =", Wall_time, "sec"

def distribute_Chunks_and_Mesh(params_simu, simuDirName):
    num_procs = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    # creating the mesh
    target_mesh = MeshClass(geoDirName, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
    if my_id==0:
        target_mesh.constructFromSavedArrays(os.path.join(tmpDirName, "mesh"))
        N_nearBlockDiag, N_near, N_nearPerCube = Z_near_size_computation(target_mesh.cubes_lists_RWGsNumbers, target_mesh.cubesNeighborsIndexes)
    else:
        N_nearPerCube = ['blabla']
    N_nearPerCube = MPI.COMM_WORLD.bcast(N_nearPerCube)
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'r')
    variables = cPickle.load(file)
    file.close()
    chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers = Z_nearChunksDistribution(params_simu.MAX_BLOCK_SIZE, N_nearPerCube, variables['C'], tmpDirName)
    chunkNumber_to_cubesNumbers = MPI.COMM_WORLD.bcast(chunkNumber_to_cubesNumbers)
    cubeNumber_to_chunkNumber = MPI.COMM_WORLD.bcast(cubeNumber_to_chunkNumber)
    chunkNumber_to_processNumber = MPI.COMM_WORLD.bcast(chunkNumber_to_processNumber)
    processNumber_to_ChunksNumbers = MPI.COMM_WORLD.bcast(processNumber_to_ChunksNumbers)
    # distributing chunks of the mesh
    scatterMesh(target_mesh, processNumber_to_ChunksNumbers, chunkNumber_to_cubesNumbers, tmpDirName, my_id, num_procs)
    del target_mesh
    # we now dump-pickle the necessary variables
    variables['chunkNumber_to_cubesNumbers'] = chunkNumber_to_cubesNumbers
    variables['cubeNumber_to_chunkNumber'] = cubeNumber_to_chunkNumber
    variables['chunkNumber_to_processNumber'] = chunkNumber_to_processNumber
    variables['processNumber_to_ChunksNumbers'] = processNumber_to_ChunksNumbers
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'w')
    cPickle.dump(variables, file)
    file.close()

if __name__=='__main__':
    sys.path.append(os.path.abspath('.'))
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir
    if simuDirName==None:
        simuDirName = '.'
    from simulation_parameters import *
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        distribute_Chunks_and_Mesh(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()
