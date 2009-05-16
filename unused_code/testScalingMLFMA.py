from mpi4py import *
#MPI.Init()
import sys, os
from scipy import zeros, ones, sqrt, array, arange, dot
from MLFMA import run_MLFMA
from PyGmsh import executeGmsh, write_geo
from EM_constants import *
from mesh_functions_seb import MeshClass
from loadParamsSimu import loadParamsSimu

class resClass:
    def __init__(self, parameter_values):
        resClass.parameter_values = parameter_values
        resClass.numberOfRWGs = zeros(len(resClass.parameter_values),'i')
        resClass.CPU_initializationTimes = zeros(len(resClass.parameter_values)).astype('f')
        resClass.CPU_resolutionTimes = zeros(len(resClass.parameter_values)).astype('f')
        resClass.Wall_initializationTimes = zeros(len(resClass.parameter_values)).astype('f')
        resClass.Wall_resolutionTimes = zeros(len(resClass.parameter_values)).astype('f')
        resClass.numbersOfMatvecs = zeros(len(resClass.parameter_values),'i')
        resClass.numbersOfGMRESIterations = zeros(len(resClass.parameter_values),'i')
        resClass.MLFMA_RCS = zeros(len(resClass.parameter_values)).astype('F')

    def save(self, path, fileName):
        saveRes = open(os.path.join(path, fileName), 'w')
        saveRes.write('parameter_values = array(' + str(resClass.parameter_values) + ')\n')
        saveRes.write('numberOfRWGs = array(' + str(resClass.numberOfRWGs.tolist()) + ')\n')
        saveRes.write('CPU_initializationTimes = array(' + str(resClass.CPU_initializationTimes.tolist()) + ')\n')
        saveRes.write('CPU_resolutionTimes = array(' + str(resClass.CPU_resolutionTimes.tolist()) + ')\n')
        saveRes.write('Wall_initializationTimes = array(' + str(resClass.Wall_initializationTimes.tolist()) + ')\n')
        saveRes.write('Wall_resolutionTimes = array(' + str(resClass.Wall_resolutionTimes.tolist()) + ')\n')
        saveRes.write('numbersOfMatvecs = array(' + str(resClass.numbersOfMatvecs.tolist()) + ')\n')
        saveRes.write('numbersOfGMRESIterations = array(' + str(resClass.numbersOfGMRESIterations.tolist()) + ')\n')
        saveRes.write('MLFMA_RCS = array(' + str(resClass.MLFMA_RCS.tolist()) + ')\n')
        saveRes.close()


if __name__=="__main__":
    #######################################################################
    # initialization
    #######################################################################
    os.system("rm -rf tmp*")
    # now we load the parameters of the simulation
    filename = "./parameters/parameters.ini"
    params_simu = loadParamsSimu(filename)

    # and now the scaling tests
    parameter_values = [0.1, 0.11, 0.12]
    #parameter_values = [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.2]
    res = resClass(parameter_values)
    for index in range(len(res.parameter_values)):
        fileName = os.path.join(params_simu.pathToTarget, params_simu.targetName)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lc', c/params_simu.f/10.)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lx', res.parameter_values[index])
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'ly', res.parameter_values[index])
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lz', res.parameter_values[index])
        executeGmsh(params_simu.pathToTarget, params_simu.targetName, 0)
        targetDimensions_scaling_factor = 1.0

        target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction)
        N_RWG = target_mesh.E
        sys.stdout.write("average edge length = %.5s" %str(target_mesh.average_edge_length/(c/params_simu.f)) + " lambdas \n")
        ## cubes computation
        target_mesh.cubes_data_computation(a)
        target_mesh.write_cubes()
        # target_MLFMA construction
        target_MLFMA = Target_MLFMA(num_proc, my_id, N_RWG, T, C, N_levels, a, CFIE, params_simu)
        # computation of near field elements
        Wall_time_Z_near_Mg_computation, CPU_time_Z_near_Mg_computation = computeZ_near(tmpDirName, a, k, w, N_levels, my_id, num_proc, chunkNumber_to_cubesNumbers, cubeNumber_to_chunkNumber, chunkNumber_to_processNumber, processNumber_to_ChunksNumbers, CFIE, params_simu)
        w = 2. * pi * params_simu.f
        k = w * sqrt(eps_0*params_simu.eps_r*mu_0*params_simu.mu_r) + 1.j * 0.
        J_dip = array([params_simu.J_src_x[0], params_simu.J_src_y[0], params_simu.J_src_z[0]], 'D')
        r_dip = array([params_simu.r_src_x[0], params_simu.r_src_y[0], params_simu.r_src_z[0]], 'f')
        print target_MLFMA.CPU_time_Target_MLFMA_making, "CPU time (seconds) for constructing MLFMA data"
        del target_mesh # we gain memory here
        target_MLFMA.solve(params_simu)
        print "N RWG =", N_RWG
        sys.stdout.write("average edge length = %.5s" %str(target_MLFMA.average_edge_length) + " m = %.5s" %str(target_MLFMA.average_edge_length/(c/params_simu.f)) + " lambdas \n")
        print target_MLFMA.CPU_time_Target_MLFMA_making, "CPU time (seconds) for constructing MLFMA data"
        print target_MLFMA.CPU_time_Target_MLFMA_resolution, "CPU time (seconds) for MLFMA iterations and solution. Iterations =", target_MLFMA.NIterMLFMA
        print target_MLFMA.CPU_time_Target_MLFMA_resolution/target_MLFMA.numberOfMatvecs, "CPU time (seconds) per MLFMA matvec"
        print target_MLFMA.CPU_time_Target_MLFMA_making + target_MLFMA.CPU_time_Target_MLFMA_resolution, "CPU time (seconds) for complete MLFMA solution"
        print target_MLFMA.Wall_time_Target_MLFMA_making, "Wall time (seconds) for constructing MLFMA data"
        print target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for MLFMA iterations and solution. Iterations =", target_MLFMA.NIterMLFMA
        print target_MLFMA.Wall_time_Target_MLFMA_resolution/target_MLFMA.numberOfMatvecs, "Wall time (seconds) per MLFMA matvec"
        print target_MLFMA.Wall_time_Target_MLFMA_making + target_MLFMA.Wall_time_Target_MLFMA_resolution, "Wall time (seconds) for complete MLFMA solution"
        print "computing the RCS..."
        target_mesh = MeshClass(params_simu.pathToTarget, params_simu.targetName, targetDimensions_scaling_factor, params_simu.z_offset, params_simu.languageForMeshConstruction)
        J_obs = array([params_simu.J_obs_x[0], params_simu.J_obs_y[0], params_simu.J_obs_z[0]], 'D')
        r_obs = array([params_simu.r_obs_x[0], params_simu.r_obs_y[0], params_simu.r_obs_z[0]], 'f')
        target_MLFMA.computeV_TE(target_mesh, J_obs, r_obs, params_simu.EXCITATION)
        res.MLFMA_RCS[index] = dot(target_MLFMA.ZI_MLFMA, target_MLFMA.V_TE)
        print "MLFMA RCS =", res.MLFMA_RCS[index]

        res.numberOfRWGs[index] = N_RWG
        res.CPU_initializationTimes[index] = target_MLFMA.CPU_time_Target_MLFMA_making
        res.CPU_resolutionTimes[index] = target_MLFMA.CPU_time_Target_MLFMA_resolution
        res.Wall_initializationTimes[index] = target_MLFMA.Wall_time_Target_MLFMA_making
        res.Wall_resolutionTimes[index] = target_MLFMA.Wall_time_Target_MLFMA_resolution
        res.numbersOfMatvecs[index] = target_MLFMA.numberOfMatvecs
        res.numbersOfGMRESIterations[index] = target_MLFMA.NIterMLFMA
        res.CPU_matvecTimes = res.CPU_resolutionTimes/res.numbersOfMatvecs
        res.CPU_totalTimes = res.CPU_initializationTimes + res.CPU_resolutionTimes
        res.save('./result', 'resultScaling' + targetName + '.txt')

    from pylab import plot, rc, subplot, xlabel, ylabel, legend, xticks, yticks, grid, gca, setp, show, title
    #rc('text', usetex=True)
    FontSize=20
    LineWidth=1
    subplot(211)
    plot(res.numberOfRWGs, res.CPU_initializationTimes, 'bo-', res.numberOfRWGs, res.CPU_resolutionTimes,'ro-',res.numberOfRWGs,res.CPU_totalTimes,'ko-')
    xlabel(r'Number of RWG functions',fontsize=FontSize+2)
    ylabel(r'time (s)',fontsize=FontSize+2)
    legend(['init','iter','total'],2)
    xticks(fontsize=FontSize)
    yticks(fontsize=FontSize)
    grid(True)
    subplot(212)
    plot(res.numberOfRWGs, res.CPU_matvecTimes, 'bo-')
    xlabel(r'Number of RWG functions',fontsize=FontSize+2)
    ylabel(r'matvec time (s)',fontsize=FontSize+2)
    xticks(fontsize=FontSize)
    yticks(fontsize=FontSize)
    grid(True)
    show()
    #MPI.Finalize()

