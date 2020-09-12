import sys, time, os
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from EM_constants import *
from meshClass import MeshClass
from MoMPostProcessing import *
from ReadWriteBlitzArray import readIntFromDisk, readFloatFromDisk, read1DBlitzArrayFromDisk


def read_CPU_time(filename):
    file = open(filename, 'r')
    CPU_time_tmp = file.readlines()
    file.close()
    CPU_time = 0.0
    for line in CPU_time_tmp:
        if 'real' in line:
            CPU_time = float(line.split()[1])
    return CPU_time

def print_times(params_simu, simuDirName):
    my_id = 0
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    # final moves
    file = open(os.path.join(tmpDirName, 'pickle', 'variables.txt'), 'rb')
    variables = cPickle.load(file)
    file.close()
    numberOfMatvecs = readIntFromDisk(os.path.join(tmpDirName,'iterative_data/numberOfMatvecs.txt'))
    NIterMLFMA = readIntFromDisk(os.path.join(tmpDirName,'iterative_data/iter.txt'))
    average_RWG_length = readFloatFromDisk(os.path.join(tmpDirName,'mesh/average_RWG_length.txt'))
    # CPU_time_GMSH
    CPU_time_GMSH = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_GMSH.txt'))
    # CPU_time_read_MLFMA_mesh_part1
    CPU_time_read_MLFMA_mesh_part1 = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_read_MLFMA_mesh_part1.txt'))
    # CPU_time_mesh_functions_seb
    CPU_time_mesh_functions_seb = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_mesh_functions_seb.txt'))
    # CPU_time_read_MLFMA_mesh_part2
    CPU_time_read_MLFMA_mesh_part2 = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_read_MLFMA_mesh_part2.txt'))
    # CPU_time_distribute_Z_cubes
    CPU_time_distribute_Z_cubes = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_distribute_Z_cubes.txt'))
    # CPU_time_compute_Z_near
    CPU_time_compute_Z_near = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_compute_Z_near.txt'))
    # CPU_time_compute_SAI_precond
    CPU_time_compute_SAI_precond = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_compute_SAI_precond.txt'))
    # CPU_time_communicateZnearBlocks
    CPU_time_communicateZnearBlocks = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_communicateZnearBlocks.txt'))
    # CPU_time_RWGs_renumbering
    CPU_time_RWGs_renumbering = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_RWGs_renumbering.txt'))
    # CPU_time_MLFMA
    CPU_time_MLFMA = read_CPU_time(os.path.join(simuDirName,'result/CPU_time_MLFMA.txt'))

    print("N RWG = " + str(variables['N_RWG']))
    sys.stdout.write("average RWG length = %.5s" %str(average_RWG_length) + " m = lambda / %.9s" %str((c/params_simu.f)/average_RWG_length) + " \n")
    print(str(CPU_time_GMSH) + " CPU time (seconds) for GMSH meshing")
    print(str(CPU_time_read_MLFMA_mesh_part1) + " CPU time (seconds) for read_MLFMA_mesh_part1")
    print(str(CPU_time_mesh_functions_seb) + " CPU time (seconds) for mesh_functions_seb")
    print(str(CPU_time_read_MLFMA_mesh_part2) + " CPU time (seconds) for read_MLFMA_mesh_part2")
    print(str(CPU_time_distribute_Z_cubes) + " CPU time (seconds) for distribute_Z_cubes")
    print(str(variables['CPU_time_distribute_ZChunks_and_cubes']) + " CPU time (seconds) for distributing Z chunks and cubes")
    #print(str(variables['Wall_time_distribute_ZChunks_and_cubes']) + " Wall time (seconds) for distributing Z chunks and cubes")
    print(str(CPU_time_compute_Z_near) + " CPU time (seconds) for compute_Z_near (C++)")
    print(str(variables['CPU_time_Z_near_computation']) + " CPU time (seconds) for constructing Z_CFIE_near")
    #print(str(variables['Wall_time_Z_near_computation']) + " Wall time (seconds) for constructing Z_CFIE_near")
    print(str(CPU_time_communicateZnearBlocks) + " CPU time (seconds) for communicateZnearBlocks")
    print(str(CPU_time_compute_SAI_precond) + " CPU time (seconds) for computing SAI precond (C++)")
    #print(str(variables['CPU_time_Mg_computation']) + " CPU time (seconds) for constructing SAI precond")
    #print(str(variables['Wall_time_Mg_computation']) + " Wall time (seconds) for constructing SAI precond")
    print(str(CPU_time_RWGs_renumbering) + " CPU time (seconds) for RWGs_renumbering")
    print(str(CPU_time_MLFMA) + " CPU time (seconds) for MLFMA iterations and solution. Iterations = " + str(NIterMLFMA))
    #print(target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for MLFMA iterations and solution. Iterations =", target_MLFMA.NIterMLFMA)
    if numberOfMatvecs>0:
        print(str(CPU_time_MLFMA/numberOfMatvecs) + " CPU time (seconds) per MLFMA matvec")
        #print(target_MLFMA.Wall_time_Target_MLFMA_resolution/target_MLFMA.numberOfMatvecs, "Wall time (seconds) per MLFMA matvec")
    print(str(CPU_time_read_MLFMA_mesh_part1 + CPU_time_mesh_functions_seb + CPU_time_read_MLFMA_mesh_part2 + CPU_time_distribute_Z_cubes + variables['CPU_time_distribute_ZChunks_and_cubes'] + CPU_time_compute_Z_near + variables['CPU_time_Z_near_computation'] + CPU_time_communicateZnearBlocks + variables['CPU_time_Mg_computation'] + CPU_time_compute_SAI_precond + CPU_time_RWGs_renumbering + CPU_time_MLFMA) + " CPU time (seconds) for complete MLFMA solution (GMSH excluded).")
    #print(Wall_time_Z_near_computation + Wall_time_Mg_computation + target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for complete MLFMA solution")
    if params_simu.CURRENTS_VISUALIZATION:
        computeCurrentsVisualization(params_simu, variables, simuDirName)
    if params_simu.SAVE_CURRENTS_CENTROIDS:
        saveCurrentsCentroids(params_simu, simuDirName)

def computeCurrentsVisualization(params_simu, variables, simuDirName):
    my_id = 0
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    if (my_id == 0):
        target_mesh = MeshClass(geoDirName, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join(tmpDirName, "mesh"))
        N_RWG = readIntFromDisk(os.path.join(tmpDirName, "mesh", "N_RWG.txt"))
        CURRENTS_VISUALIZATION = params_simu.CURRENTS_VISUALIZATION and (N_RWG<2e6) and (params_simu.BISTATIC == 1)
        if CURRENTS_VISUALIZATION:
            print("............Constructing visualisation of currents")
            if (N_RWG<5e4):
                nbTimeSteps = 48
            else:
                nbTimeSteps = 1
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            norm_J_centroids_triangles = normJMCentroidsTriangles(J_centroids_triangles)
            norm_J_centroids_triangles_timeDomain = normJMCentroidsTriangles_timeDomain(J_centroids_triangles, variables['w'], nbTimeSteps)
            write_VectorFieldTrianglesCentroidsGMSH(os.path.join(simuDirName, 'result', params_simu.targetName) + '.J_centroids_triangles.pos', real(J_centroids_triangles), target_mesh)
            write_ScalarFieldTrianglesCentroidsGMSH(os.path.join(simuDirName, 'result', params_simu.targetName) + '.norm_J_centroids_triangles_timeDomain.pos', norm_J_centroids_triangles_timeDomain, target_mesh)
            write_ScalarFieldTrianglesCentroidsGMSH(os.path.join(simuDirName, 'result', params_simu.targetName) + '.norm_J_centroids_triangles.pos', norm_J_centroids_triangles, target_mesh)
            print("............end of currents visualisation construction. Open with gmsh the *.pos files in result directory.")
        else:
            "Error in the computeCurrentsVisualization routine. Probably you required currents visualization computation for a mesh too big"
            sys.exit(1)

def saveCurrentsCentroids(params_simu, simuDirName):
    my_id = 0
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    geoDirName = os.path.join(simuDirName, 'geo')
    if (my_id == 0):
        target_mesh = MeshClass(geoDirName, params_simu.targetName, params_simu.targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction, params_simu.meshFormat, params_simu.meshFileTermination)
        # target_mesh.constructFromGmshFile()
        target_mesh.constructFromSavedArrays(os.path.join(tmpDirName, "mesh"))
        SAVE_CURRENTS_CENTROIDS = (params_simu.SAVE_CURRENTS_CENTROIDS and (params_simu.BISTATIC == 1))
        if SAVE_CURRENTS_CENTROIDS:
            print("............saving currents at centroids of triangles")
            nbTimeSteps = 1
            I_coeff = read1DBlitzArrayFromDisk(os.path.join(tmpDirName, 'ZI/ZI.txt'), 'F')
            J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)
            write_VectorFieldTrianglesCentroids(os.path.join(simuDirName, 'result', 'J_centroids_triangles.txt'), J_centroids_triangles, target_mesh)
            print("............end of saving of currents. Located in result directory.")


