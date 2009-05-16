import sys, os
from mpi4py import *
from MLFMA import getField, print_times
from scipy import cos, sin, conj, log10, real, sum, dot, pi, sqrt, exp
from scipy import array, arange, zeros, ones
from V_EH import G_EJ_G_HJ
from EM_constants import *
from ReadWriteBlitzArray import readASCIIBlitzFloatArray1DFromDisk, readBlitzArrayFromDisk, readASCIIBlitzFloatArray2DFromDisk

def monostatic_SAR(params_simu):
    # data generation of Volakis article, IEEE Antennas and Propagation Magazine, December 1992
    # for this purpose we use the MLFMA routines
    my_id = MPI.COMM_WORLD.Get_rank()
    if my_id==0:
        r_SAR = readASCIIBlitzFloatArray2DFromDisk('./result/r_SAR.txt')
        SAR_RCS_HH = readASCIIBlitzFloatArray1DFromDisk('./result/SAR_RCS_HH_ASCII.txt')
        SAR_RCS_HV = readASCIIBlitzFloatArray1DFromDisk('./result/SAR_RCS_HV_ASCII.txt')
        SAR_RCS_VV = readASCIIBlitzFloatArray1DFromDisk('./result/SAR_RCS_VV_ASCII.txt')
    else:
        r_SAR = zeros((1, 1))
        SAR_RCS_HH = zeros(1)
        SAR_RCS_HV = zeros(1)
        SAR_RCS_VV = zeros(1)
    return SAR_RCS_HH, SAR_RCS_VV, SAR_RCS_HV, r_SAR

def monostatic_RCS(params_simu):
    # data generation of Volakis article, IEEE Antennas and Propagation Magazine, December 1992
    # for this purpose we use the MLFMA routines
    my_id = MPI.COMM_WORLD.Get_rank()
    if my_id==0:
        phis_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk('./result/phis_far_field_ASCII.txt')
        thetas_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk('./result/thetas_far_field_ASCII.txt')
        RCS_HH = readASCIIBlitzFloatArray2DFromDisk('./result/RCS_HH_ASCII.txt')
        RCS_HV = readASCIIBlitzFloatArray2DFromDisk('./result/RCS_HV_ASCII.txt')
        RCS_VV = readASCIIBlitzFloatArray2DFromDisk('./result/RCS_VV_ASCII.txt')
    else:
        phis_far_field = zeros(1)
        thetas_far_field = zeros(1)
        RCS_HH = zeros((1, 1))
        RCS_HV = zeros((1, 1))
        RCS_VV = zeros((1, 1))
    return RCS_HH, RCS_VV, RCS_HV, thetas_far_field, phis_far_field

def bistatic_RCS(params_simu):
    my_id = MPI.COMM_WORLD.Get_rank()
    if my_id==0:
        phis_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk('./result/phis_far_field_ASCII.txt')
        thetas_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk('./result/thetas_far_field_ASCII.txt')
        e_phi = readBlitzArrayFromDisk('./result/e_phi_far_Binary.txt', thetas_far_field.shape[0], phis_far_field.shape[0], 'F')
        e_theta = readBlitzArrayFromDisk('./result/e_theta_far_Binary.txt', thetas_far_field.shape[0], phis_far_field.shape[0], 'F')
        p_scatt_theta = real(e_theta * conj(e_theta))
        p_scatt_phi = real(e_phi*conj(e_phi))
        R_cube_center = readASCIIBlitzFloatArray1DFromDisk('./tmp' + str(my_id) +  '/octtree_data/big_cube_center_coord.txt')
        w = 2. * pi * params_simu.f
        eps, mu = params_simu.eps_r * eps_0, params_simu.mu_r * mu_0
        k = w * sqrt(mu * eps) # the wavenumber
        if params_simu.EXCITATION=='dipole':
            r_dip_src = array([params_simu.r_src_x, params_simu.r_src_y, params_simu.r_src_z], 'd')
            G_EJ_inc, G_HJ_inc = G_EJ_G_HJ(r_dip_src, R_cube_center, eps, mu, k)
            J_dip_src = array([params_simu.J_src_x, params_simu.J_src_y, params_simu.J_src_z], 'D')
            E_inc = dot(G_EJ_inc, J_dip_src)
        elif params_simu.EXCITATION=='plane':
            E_inc = array([params_simu.E_inc_theta, params_simu.E_inc_phi], 'D')
        else:
            pass
        P_inc = real(dot(E_inc, conj(E_inc)))
        sigma_phi = p_scatt_phi/(P_inc * 4.0 * pi)
        sigma_theta = p_scatt_theta/(P_inc * 4.0 * pi)
    else:
        sigma_theta, sigma_phi = ones((1,1)), ones((1,1))
        thetas_far_field, phis_far_field = zeros(1), zeros(1)
    return sigma_theta, sigma_phi, thetas_far_field, phis_far_field

def testVolakisMLFMA(params_simu):

    # Testing of the code for PEC targets
    #
    # For more info about the tests, see:
    # [1] A.C. Woo et al., "Benchmark Plate Radar Targets for the Validation
    #     of Computational Electromagnetics Programs", IEEE Antennas and 
    #     Propagation Magazine, Vol. 34, No. 6, December 1992
    # [2] A.C. Woo et al., "Benchmark Radar Targets for the Validation
    #     of Computational Electromagnetics Programs", IEEE Antennas and 
    #     Propagation Magazine, Vol. 35, No. 1, February 1993
    #
    # The code can also be compared to the RCS computed in:
    # [3] Gang Kan et al., "A Novel Grid-Robust Higher Order Vector Basis Function
    #     for the Method of Moments", IEEE Trans. Ant. Prop., vol. 49, No. 6, 
    #     June 2001, pp 908--915
    #
    # There are tests on 2 types of targets: plane targets and 3D targets.
    # These targets have been defined by the Electromagnetic Code Consortium.
    #
    # The EMCC targets described in [1], [2] have been implemented 
    # and are in the "Puma-EM/geo" directory.
    #
    # certain data will be automatically created/overwritten for you, i.e.:
    #   - params_simu.thetas_obs
    #
    # However, you still have to set your own desired frequency, this for allowing you to easily 
    # generate the double-frequency results of reference [2] (or create new results on your own).
    # 
    # The frequency doesn't matter for plate targets, as their dimensions are calculated in wavelengths.
    #
    # For the 3D targets, the frequencies used in [2] are:
    #   - 1.19 GHz, 7 GHz and 9.92 GHz for the almond (also 5 GHz in [3], Fig. 8)
    #   - 1.18 GHz and 9 GHz for the ogive
    #   - 1.57 GHz and 9 GHz for the double ogive
    #   - 0.869 GHz and 9 GHz for the cone-sphere and the cone-gap-sphere

    params_simu.COMPUTE_RCS_HH = 1
    params_simu.COMPUTE_RCS_VV = 1
    params_simu.MONOSTATIC_RCS = 1 # we have to compute the monostatic RCS
    params_simu.MOM_FULL_PRECISION = 1 # we need precision
    params_simu.lc_factor = 1.0/10.0 # we need precision
    params_simu.dB = 1.0
    if ('EMCC' in params_simu.targetName) and ('plate' in params_simu.targetName):
        params_simu.lx = c/params_simu.f
        params_simu.USER_DEFINED_THETAS_OBS = 1
        params_simu.thetas_obs = array([80.0/180.0]) * pi
        params_simu.USER_DEFINED_PHIS_OBS = 0
        params_simu.dB = 1.0 / (c/params_simu.f)**2
    elif 'EMCC' in params_simu.targetName:
        params_simu.USER_DEFINED_THETAS_OBS = 1
        params_simu.thetas_obs = array([90.0/180.0]) * pi
        params_simu.USER_DEFINED_PHIS_OBS = 0
    RCS_HH, RCS_VV, RCS_HV, thetas_far_field, phis_far_field = monostatic_RCS(params_simu)
    return RCS_HH, RCS_VV, RCS_HV, thetas_far_field, phis_far_field

if __name__=='__main__':
    MPI.Init()
    my_id = MPI.COMM_WORLD.Get_rank()
    sys.path.append(os.path.abspath('.'))
    # the simulation itself
    from simulation_parameters import *
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        print_times(params_simu)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    if params_simu.MONOSTATIC_RCS==1:
        RCS_HH, RCS_VV, RCS_HV, thetas_far_field, phis_far_field = monostatic_RCS(params_simu)
        if (my_id==0):
            nameOfFileToSaveTo = os.path.join('./result', "simulation_parameters.txt") 
            params_simu.saveTo(nameOfFileToSaveTo)
        if (my_id==0) and params_simu.SHOW_FIGURE:
            from pylab import rc, plot, show, xlabel, ylabel, xticks, yticks, grid, legend, title
            #rc('text', usetex=True)
            FontSize=18
            LineWidth=2
            LEGEND = []
            if (params_simu.COMPUTE_RCS_HH==1):
                plot(phis_far_field, 10 * log10(RCS_HH[0]), 'bo-', linewidth = LineWidth)
                LEGEND.append(r'RCS$_{HH}$')
            if (params_simu.COMPUTE_RCS_VV==1):
                plot(phis_far_field, 10 * log10(RCS_VV[0]), 'rs-', linewidth = LineWidth)
                LEGEND.append(r'RCS$_{VV}$')
            if (params_simu.COMPUTE_RCS_HV==1):
                plot(phis_far_field, 10 * log10(RCS_HV[0]), 'gv-', linewidth = LineWidth)
                LEGEND.append(r'RCS$_{HV}$')
            figureTitle = ""
            for elem in params_simu.targetName.split("_"):
                figureTitle += elem + " "
            figureTitle += ", f = " + str(params_simu.f/1.e9) + " GHz"
            title(figureTitle,fontsize=FontSize+2)
            xlabel(r'azimuthal angle $\phi$',fontsize=FontSize+2)
            ylabel(r'$\sigma = 4 \pi R^2  P_s/P_i$ [dB]',fontsize=FontSize+2)
            legend(LEGEND)
            xticks(fontsize=FontSize)
            yticks(fontsize=FontSize)
            grid(True)
            show()
    if params_simu.BISTATIC==1:
        params_simu.VERBOSE = 1
        sigma_theta, sigma_phi, thetas_far_field, phis_far_field = bistatic_RCS(params_simu)
        if (my_id==0):
            E_field = getField("./result/E_obs.txt")
            r_obs = readASCIIBlitzFloatArray2DFromDisk("./result/r_obs.txt")
            print
            print "MLFMA E_field ="
            print E_field
            print
            print "at observation points r_obs ="
            print r_obs
            print
            print "See 'E_obs.txt' and 'r_obs.txt' in './result' directory for the results of the computation, where you will also find the far field values."
            print
            nameOfFileToSaveTo = os.path.join('./result', "simulation_parameters.txt") 
            params_simu.saveTo(nameOfFileToSaveTo)
        if (my_id==0) and params_simu.SHOW_FIGURE:
            from pylab import rc, plot, show, xlabel, ylabel, xticks, yticks, grid, legend, title
            #rc('text', usetex=True)
            FontSize=18
            LineWidth=2
            LEGEND = []
            if (params_simu.COMPUTE_RCS_HH==1):
                plot(phis_far_field, 10 * log10(sigma_phi[0]), 'bo-', linewidth = LineWidth)
                LEGEND.append(r'$\sigma_{\phi}$')
            if (params_simu.COMPUTE_RCS_VV==1):
                plot(phis_far_field, 10 * log10(sigma_theta[0]), 'rs-', linewidth = LineWidth)
                LEGEND.append(r'$\sigma_{\theta}$')
            figureTitle = ""
            for elem in params_simu.targetName.split("_"):
                figureTitle += elem + " "
            figureTitle += ", f = " + str(params_simu.f/1.e9) + " GHz"
            title(figureTitle,fontsize=FontSize+2)
            xlabel(r'azimuthal angle $\phi$',fontsize=FontSize+2)
            ylabel(r'$\sigma = 4 \pi R^2  P_s/P_i$ [dB]',fontsize=FontSize+2)
            legend(LEGEND)
            xticks(fontsize=FontSize)
            yticks(fontsize=FontSize)
            grid(True)
            show()

    if params_simu.MONOSTATIC_SAR==1:
        SAR_RCS_HH, SAR_RCS_VV, SAR_RCS_HV, r_SAR = monostatic_SAR(params_simu)

    MPI.Finalize()

