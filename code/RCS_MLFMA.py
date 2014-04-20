import sys, os, argparse
from MLFMA import getField, print_times
from scipy import cos, sin, conj, log10, real, sum, dot, pi, sqrt, exp
from scipy import array, arange, zeros, ones
from V_EH import G_EJ_G_HJ
from EM_constants import *
from ReadWriteBlitzArray import *
from read_dipole_excitation import read_dipole_excitation, read_input_angles

def monostatic_SAR(params_simu, simuDirName):
    r_SAR = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/r_SAR.txt'))
    SAR_RCS_HH = readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/SAR_RCS_HH_ASCII.txt'))
    SAR_RCS_HV = readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/SAR_RCS_HV_ASCII.txt'))
    SAR_RCS_VH = readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/SAR_RCS_VH_ASCII.txt'))
    SAR_RCS_VV = readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/SAR_RCS_VV_ASCII.txt'))
    return SAR_RCS_HH, SAR_RCS_VV, SAR_RCS_HV, SAR_RCS_VH, r_SAR

def monostatic_RCS(params_simu, simuDirName):
    phis_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/phis_far_field_ASCII.txt'))
    thetas_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/thetas_far_field_ASCII.txt'))
    RCS_HH = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/RCS_HH_ASCII.txt'))
    RCS_HV = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/RCS_HV_ASCII.txt'))
    RCS_VH = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/RCS_VH_ASCII.txt'))
    RCS_VV = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/RCS_VV_ASCII.txt'))
    return RCS_HH, RCS_VV, RCS_HV, RCS_VH, thetas_far_field, phis_far_field

def bistatic_RCS(params_simu, simuDirName):
    phis_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/phis_far_field_ASCII.txt'))
    thetas_far_field = 180./pi * readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'result/thetas_far_field_ASCII.txt'))
    e_phi = readBlitzArrayFromDisk(os.path.join(simuDirName, 'result/e_phi_far_Binary.txt'), thetas_far_field.shape[0], phis_far_field.shape[0], 'F')
    e_theta = readBlitzArrayFromDisk(os.path.join(simuDirName, 'result/e_theta_far_Binary.txt'), thetas_far_field.shape[0], phis_far_field.shape[0], 'F')
    p_scatt_theta = real(e_theta * conj(e_theta))
    p_scatt_phi = real(e_phi*conj(e_phi))
    R_cube_center = readASCIIBlitzFloatArray1DFromDisk(os.path.join(simuDirName, 'tmp' + str(0) +  '/octtree_data/big_cube_center_coord.txt'))
    w = 2. * pi * params_simu.f
    eps, mu = params_simu.eps_r * eps_0, params_simu.mu_r * mu_0
    k = w * sqrt(mu * eps) # the wavenumber
    P_inc = 0.0
    if (params_simu.BISTATIC_EXCITATION_DIPOLES == 1):
        J_src, r_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME)
        r_dip_src = r_src[0,:]
        J_dip_src = J_src[0,:]
        G_EJ_inc, G_HJ_inc = G_EJ_G_HJ(r_dip_src, R_cube_center, eps, mu, k)
        E_inc = dot(G_EJ_inc, J_dip_src)
        P_inc += real(dot(E_inc, conj(E_inc)))
    if (params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1):
        E_inc = array([params_simu.E_inc_theta, params_simu.E_inc_phi], 'D')
        P_inc += real(dot(E_inc, conj(E_inc)))            
    if (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1):
        print("WARNING: you have dipole and plane wave excitation simultaneously. Is it what you intended??")
    sigma_phi = p_scatt_phi/(P_inc * 4.0 * pi)
    sigma_theta = p_scatt_theta/(P_inc * 4.0 * pi)
    return sigma_theta, sigma_phi, thetas_far_field, phis_far_field

if __name__=='__main__':
    sys.path.append(os.path.abspath('.'))
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
    exec('from ' + simuParams + ' import *')
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        print_times(params_simu, simuDirName)
    else:
        print("you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings.")
        sys.exit(1)
    if params_simu.MONOSTATIC_RCS==1:
        if params_simu.ANGLES_FROM_FILE == 1:
            monostatic_angles = 180./pi * readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, 'result/monostatic_angles_ASCII.txt'))
            print("monostatic angles (degrees) = ")
            print(monostatic_angles)
        else:
            RCS_HH, RCS_VV, RCS_HV, RCS_VH, thetas_far_field, phis_far_field = monostatic_RCS(params_simu, simuDirName)
            nameOfFileToSaveTo = os.path.join(simuDirName, 'result', "simulation_parameters.txt") 
            params_simu.saveTo(nameOfFileToSaveTo)
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
            if (params_simu.COMPUTE_RCS_VH==1):
                plot(phis_far_field, 10 * log10(RCS_VH[0]), 'gv-', linewidth = LineWidth)
                LEGEND.append(r'RCS$_{VH}$')
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
        # user supplied R_OBS
        if params_simu.BISTATIC_R_OBS==1:
            E_field = getField(os.path.join(simuDirName, "result", "E_obs.txt"))
            r_obs = readASCIIBlitzFloatArray2DFromDisk(os.path.join(simuDirName, "result", "r_obs.txt"))
            print('\r')
            print("MLFMA E_field = ") 
            print(E_field)
            print('\r')
            print("at observation points r_obs = ")
            print(r_obs)
            print('\r')
            print("See 'E_obs.txt' and 'r_obs.txt' in './result' directory for the results of the computation, where you will also find the far field values.")
            print('\r')
        nameOfFileToSaveTo = os.path.join(simuDirName, 'result', "simulation_parameters.txt") 
        params_simu.saveTo(nameOfFileToSaveTo)

        # user supplied observation angles
        if (params_simu.BISTATIC_ANGLES_OBS == 1) and (params_simu.BISTATIC_ANGLES_OBS_FILENAME != ""):
            bistatic_angles_obs = read_input_angles(params_simu.BISTATIC_ANGLES_OBS_FILENAME)
            thetas_obs = bistatic_angles_obs[:,0]
            phis_obs = bistatic_angles_obs[:,1]
            # the regular far field grid
            sigma_theta, sigma_phi, thetas_far_field, phis_far_field = bistatic_RCS(params_simu, simuDirName)
            import scipy.interpolate
            spline_sigma_theta = scipy.interpolate.RectBivariateSpline(thetas_far_field, phis_far_field, sigma_theta)
            spline_sigma_phi = scipy.interpolate.RectBivariateSpline(thetas_far_field, phis_far_field, sigma_phi)
            sigma_theta_obs = spline_sigma_theta.ev(thetas_obs, phis_obs)
            sigma_phi_obs = spline_sigma_phi.ev(thetas_obs, phis_obs)
            writeASCIIBlitzArrayToDisk(sigma_theta_obs, os.path.join(simuDirName, 'result', "sigma_theta_obs.txt"))
            writeASCIIBlitzArrayToDisk(sigma_phi_obs, os.path.join(simuDirName, 'result', "sigma_phi_obs.txt"))
            writeASCIIBlitzArrayToDisk(thetas_obs, os.path.join(simuDirName, 'result', "thetas_obs.txt"))
            writeASCIIBlitzArrayToDisk(phis_obs, os.path.join(simuDirName, 'result', "phis_obs.txt"))

        # automatic far field computations
        sigma_theta, sigma_phi, thetas_far_field, phis_far_field = bistatic_RCS(params_simu, simuDirName)
        dimensions = 2
        if (len(thetas_far_field)==1) or (len(phis_far_field)==1):
            dimensions = 1
        if params_simu.SHOW_FIGURE:
            if dimensions==1:
                from pylab import rc, plot, show, xlabel, ylabel, xticks, yticks, grid, legend, title
                #rc('text', usetex=True)
                FontSize=18
                LineWidth=2
                LEGEND = []
                angles_far_field = thetas_far_field
                if (len(thetas_far_field)==1):
                    angles_far_field = phis_far_field
                if (params_simu.COMPUTE_RCS_HH==1):
                    plot(angles_far_field, 10 * log10(sigma_phi[0]), 'bo-', linewidth = LineWidth)
                    LEGEND.append(r'$\sigma_{\phi}$')
                if (params_simu.COMPUTE_RCS_VV==1):
                    plot(angles_far_field, 10 * log10(sigma_theta[0]), 'rs-', linewidth = LineWidth)
                    LEGEND.append(r'$\sigma_{\theta}$')
                figureTitle = ""
                for elem in params_simu.targetName.split("_"):
                    figureTitle += elem + " "
                figureTitle += ", f = " + str(params_simu.f/1.e9) + " GHz"
                title(figureTitle,fontsize=FontSize+2)
                if (len(thetas_far_field)==1):
                    xlabel(r'azimuthal angle $\phi$',fontsize=FontSize+2)
                else:
                    xlabel(r'elevation angle $\theta$',fontsize=FontSize+2)
                ylabel(r'$\sigma = 4 \pi R^2  P_s/P_i$ [dB]',fontsize=FontSize+2)
                legend(LEGEND)
                xticks(fontsize=FontSize)
                yticks(fontsize=FontSize)
                grid(True)
                show()
            else:
                from mpl_toolkits.mplot3d import Axes3D
                from matplotlib import cm
                from matplotlib.ticker import LinearLocator, FormatStrFormatter
                import matplotlib.pyplot as plt
                import numpy as np
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                X = phis_far_field
                Y = thetas_far_field
                X, Y = np.meshgrid(X, Y)
                surf = ax.plot_surface(X, Y, sigma_theta, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                ax.zaxis.set_major_locator(LinearLocator(10))
                ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
                fig.colorbar(surf, shrink=0.5, aspect=5)
                plt.show()

    if params_simu.MONOSTATIC_SAR==1:
        SAR_RCS_HH, SAR_RCS_VV, SAR_RCS_HV, SAR_RCS_VH, r_SAR = monostatic_SAR(params_simu, simuDirName)


