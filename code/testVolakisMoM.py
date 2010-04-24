import sys, os
from MOM import *
from mesh_functions_seb import triangles_centroids_computation, triangles_areas_normals_computation
from V_EH import computeV_EH

def testVolakisMoM(path, targetName, f, M0M_FULL_PRECISION):
    # data generation of Volakis article, IEEE Antennas and Propagation Magazine, December 1992
    # for this purpose we use the MoM routines with no FMM or MLFMA
    write_geo(path, targetName, 'lc', c/f/10.0)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    languageForMeshConstruction = "Python"
    target_mesh = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, 'GMSH', '.msh')
    target_mesh.constructFromGmshFile()
    N_RWG = target_mesh.N_RWG
    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    coeff = 0.2
    TDS_APPROX = 0
    Z_s = 0.0
    CFIE = array([coeff, 0.0, 0.0, -(1.0-coeff) * 377], 'D')
    list_of_test_edges_numbers = arange(N_RWG,dtype='i')#.astype('i')#[0::6]
    list_of_src_edges_numbers = arange(N_RWG,dtype='i')#.astype('i')#[2::6]
    print "Target_MoM instanciation"
    target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
    print "Matrix inversion..."
    target_MoM.compute_Y_CFIE()
    T = target_mesh.T
    if 'plate' in targetName:
        theta = 80.0/180.0 * pi # cone angle in radians
    else:
        theta = 90.0/180.0 * pi
    phis = arange(0., 181, 1) / 180.0 * pi
    sigmas_HH = zeros(len(phis), 'f')
    sigmas_VV = zeros(len(phis), 'f')
    for polarization in ['HH', 'VV']:
        print "polarization " + polarization
        index = 0
        r_obs = zeros(3).astype('d') # observation point for the incoming field
        for phi in phis:
            sys.stdout.write("\r" + "%.4s" %str(phis[0]/pi*180.) + " -> phi = %.4s" %str(phi/pi*180.) + " -> %.4s" %str(phis[-1]/pi*180.))
            sys.stdout.flush()
            # unit vectors computation
            r_hat = array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)], 'd');
            theta_hat = array([cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)], 'd');
            phi_hat = array([-sin(phi), cos(phi), 0.0], 'f');
            # excitation parameters computation
            J_dip_factor = -(1.0 + 0.0j)
            if polarization=='HH':
                J_dip = J_dip_factor * phi_hat # HH polarization
            else:
                J_dip = J_dip_factor * theta_hat # VV polarization
            R_dip = 300.0 * c/f # we wanna be in far field
            r_dip = R_dip * r_hat
            # excitation vector computation
            V_EH = computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, 'plane', 'F')
            V_CFIE = zeros(V_EH.shape[0], complex)
            for i in range(4):
                V_CFIE += V_EH[:,i] * CFIE[i]
            # scattered field computation
            I_CFIE = dot(target_MoM.Y_CFIE, V_CFIE)
            # we compute the scattered fields by reciprocity
            V_EH = computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, 'dipole', 'F')
            E_scatt = sum(I_CFIE * V_EH[:,0]/J_dip_factor)
            # incoming field computation
            Delta_r = (r_dip - r_obs)
            k = w/c # the wavenumber
            eps = eps_r * eps_0
            mu = mu_r * mu_0
            G_EJ_inc, G_HJ_inc = G_EJ_G_HJ(r_dip, r_obs, eps, mu, k)
            E_inc = dot(G_EJ_inc, J_dip)
            # sigma computation
            P_scatt = real(E_scatt*conj(E_scatt))
            P_inc = real(dot(E_inc, conj(E_inc)))
            if polarization=='HH':
                sigmas_HH[index] = 4*pi* sum(Delta_r**2) * P_scatt/P_inc
            else:
                sigmas_VV[index] = 4*pi* sum(Delta_r**2) * P_scatt/P_inc
            index = index + 1
        print
    return sigmas_HH, sigmas_VV

if __name__=='__main__':
    path = './geo'
    targetName = 'EMCC_almond'
    if 'plate' in targetName:
        f = 2.12e9
        dB = 1.0 / (c/f)**2
    else:
        dB = 1.0
        if 'ogive' in targetName and not 'double-ogive' in targetName:
            f = 1.18e9
        elif 'double-ogive' in targetName:
            f = 1.57e9
        elif 'almond' in targetName:
            f = 1.19e9
        elif 'cone-sphere' in targetName or 'cone-gap-sphere' in targetName:
            f = 0.869e9
        else:
            print "No target named", targetName, ". Exiting"
    MOM_FULL_PRECISION = 1
    sigmas_HH, sigmas_VV = testVolakisMoM(path, targetName, f, MOM_FULL_PRECISION)

    from pylab import rc, plot, xlabel, ylabel, legend, xticks, yticks, grid, show
    #rc('text', usetex=True)
    FontSize=18
    LineWidth=1
    plot(arange(sigmas_HH.shape[0]), 10 * log10(sigmas_HH * dB), 'b', arange(sigmas_VV.shape[0]), 10 * log10(sigmas_VV * dB), 'r--', linewidth = LineWidth)
    xlabel(r'azimuthal angle $\phi$',fontsize=FontSize+2)
    ylabel(r'$\sigma = 4 \pi R^2  P_s/P_i$',fontsize=FontSize+2)
    legend([r'$\sigma_{HH}$', r'$\sigma_{VV}$'])
    xticks(fontsize=FontSize)
    yticks(fontsize=FontSize)
    grid(True)
    show()
