import sys, os, argparse
from mpi4py import MPI
from scipy import zeros, array
from ReadWriteBlitzArray import writeScalarToDisk, writeASCIIBlitzArrayToDisk
from read_dipole_excitation import read_dipole_excitation, read_observation_points, read_input_angles

def setup_excitation(params_simu, simuDirName):
    num_proc = MPI.COMM_WORLD.Get_size()
    my_id = MPI.COMM_WORLD.Get_rank()

    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    # the observation points
    if (params_simu.r_obs_FROM_FILE == 0):
        if (len(params_simu.r_obs_x)==len(params_simu.r_obs_y)) and (len(params_simu.r_obs_x)==len(params_simu.r_obs_z)):
            N_obs_points = len(params_simu.r_obs_x)
            r_obs = zeros((N_obs_points, 3), 'd')
            for index in range(N_obs_points):
                r_obs[index,:] = array([params_simu.r_obs_x[index], params_simu.r_obs_y[index], params_simu.r_obs_z[index]], 'd')
        else:
            if (my_id==0):
                print "Error in the r_obs parameter. Check your observation points!"
            sys.exit(1)
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join(tmpDirName,'V_CFIE/r_obs.txt'))
    elif (params_simu.r_obs_FROM_FILE == 1) and (params_simu.r_obs_FILENAME != ""):
        if (my_id==0): # this file is only on processor 0
            r_obs = read_observation_points(params_simu.r_obs_FILENAME)
        else:
            r_obs = zeros((1, 3), 'd')
        r_obs = MPI.COMM_WORLD.bcast(r_obs)
        writeASCIIBlitzArrayToDisk(r_obs, os.path.join(tmpDirName,'V_CFIE/r_obs.txt'))

    # now the excitations
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_DIPOLES, os.path.join(tmpDirName,'V_CFIE/DIPOLES_EXCITATION.txt'))
    writeScalarToDisk(params_simu.BISTATIC_EXCITATION_PLANE_WAVE, os.path.join(tmpDirName,'V_CFIE/PLANE_WAVE_EXCITATION.txt'))
    writeScalarToDisk(params_simu.V_FULL_PRECISION*1, os.path.join(tmpDirName, 'V_CFIE/V_FULL_PRECISION.txt') )
    # if we have dipoles excitation AND definition of the excitation in simulation_parameters.py
    if (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_DIPOLES_FROM_FILE == 0):
        # first for electric dipoles
        N_src_points = 0
        if (len(params_simu.r_J_src_x)==len(params_simu.r_J_src_y)) and (len(params_simu.r_J_src_x)==len(params_simu.r_J_src_z)):
            N_src_points = len(params_simu.r_J_src_x)
            if len(params_simu.J_src_x)==N_src_points and len(params_simu.J_src_x)==len(params_simu.J_src_y) and len(params_simu.J_src_x)==len(params_simu.J_src_z):
                J_src = zeros((N_src_points, 3), 'D')
                r_J_src = zeros((N_src_points, 3), 'd')
                for index in range(N_src_points):
                    r_J_src[index,:] = array([params_simu.r_J_src_x[index], params_simu.r_J_src_y[index], params_simu.r_J_src_z[index]], 'd')
                    J_src[index,:] = array([params_simu.J_src_x[index], params_simu.J_src_y[index], params_simu.J_src_z[index]], 'D')
            else:
                if (my_id==0):
                    print "Error in the r_J_src parameter. Check your source points!"
                sys.exit(1)
        else:
            if (my_id==0):
                print "Error in the r_J_src parameter. Check your source points!"
            sys.exit(1)
        if N_src_points > 0:
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join(tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join(tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
        # then for magnetic dipoles
        N_src_points = 0
        if (len(params_simu.r_M_src_x)==len(params_simu.r_M_src_y)) and (len(params_simu.r_M_src_x)==len(params_simu.r_M_src_z)):
            N_src_points = len(params_simu.r_M_src_x)
            if len(params_simu.M_src_x)==N_src_points and len(params_simu.M_src_x)==len(params_simu.M_src_y) and len(params_simu.M_src_x)==len(params_simu.M_src_z):
                M_src = zeros((N_src_points, 3), 'D')
                r_M_src = zeros((N_src_points, 3), 'd')
                for index in range(N_src_points):
                    r_M_src[index,:] = array([params_simu.r_M_src_x[index], params_simu.r_M_src_y[index], params_simu.r_M_src_z[index]], 'd')
                    M_src[index,:] = array([params_simu.M_src_x[index], params_simu.M_src_y[index], params_simu.M_src_z[index]], 'D')
            else:
                if (my_id==0):
                    print "Error in the r_M_src parameter. Check your source points!"
                sys.exit(1)
        else:
            if (my_id==0):
                print "Error in the r_M_src parameter. Check your source points!"
            sys.exit(1)
        if N_src_points > 0:
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join(tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join(tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # if we have dipoles excitation AND definition of the excitation in a user-supplied file
    elif (params_simu.BISTATIC_EXCITATION_DIPOLES == 1) and (params_simu.BISTATIC_EXCITATION_DIPOLES_FROM_FILE == 1):
        if params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                J_src, r_J_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME)
            else:
                J_src, r_J_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            J_src = MPI.COMM_WORLD.bcast(J_src)
            r_J_src = MPI.COMM_WORLD.bcast(r_J_src)
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(J_src, os.path.join(tmpDirName,'V_CFIE/J_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_J_src, os.path.join(tmpDirName,'V_CFIE/r_J_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/J_DIPOLES_EXCITATION.txt'))
        if params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME != "":
            if (my_id==0): # this file is only on processor 0
                M_src, r_M_src = read_dipole_excitation(params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME)
            else:
                M_src, r_M_src = zeros((1, 3), 'D'), zeros((1, 3), 'd')
            M_src = MPI.COMM_WORLD.bcast(M_src)
            r_M_src = MPI.COMM_WORLD.bcast(r_M_src)
            writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
            writeASCIIBlitzArrayToDisk(M_src, os.path.join(tmpDirName,'V_CFIE/M_dip.txt'))
            writeASCIIBlitzArrayToDisk(r_M_src, os.path.join(tmpDirName,'V_CFIE/r_M_dip.txt'))
        else:
            writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/M_DIPOLES_EXCITATION.txt'))
    # now the plane wave excitation
    if params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1:
        writeScalarToDisk(params_simu.theta_inc, os.path.join(tmpDirName,'V_CFIE/theta_inc.txt'))
        writeScalarToDisk(params_simu.phi_inc, os.path.join(tmpDirName,'V_CFIE/phi_inc.txt'))
        E_inc = array([params_simu.E_inc_theta, params_simu.E_inc_phi], 'D')
        writeASCIIBlitzArrayToDisk(E_inc, os.path.join(tmpDirName,'V_CFIE/E_inc.txt'))
    if (params_simu.BISTATIC_EXCITATION_DIPOLES != 1) and (params_simu.BISTATIC_EXCITATION_PLANE_WAVE != 1):
        if (my_id==0):
            print "incorrect excitation choice. You have to choose dipole and/or plane wave excitation."
        sys.exit(1)

    if (params_simu.MONOSTATIC_RCS == 1) and (params_simu.ANGLES_FROM_FILE == 1) and (params_simu.ANGLES_FILENAME != ""):
        if (my_id==0): # this file is only on processor 0
            angles = read_input_angles(params_simu.ANGLES_FILENAME)
        else:
            angles = zeros((1, 2), 'd')
        angles = MPI.COMM_WORLD.bcast(angles)
        writeASCIIBlitzArrayToDisk(angles, os.path.join(tmpDirName,'V_CFIE/monostatic_angles.txt'))
        writeScalarToDisk(1, os.path.join(tmpDirName,'V_CFIE/ANGLES_FROM_FILE.txt'))
    elif (params_simu.MONOSTATIC_RCS == 1) and ((params_simu.ANGLES_FROM_FILE == 0) or (params_simu.ANGLES_FILENAME == "")):
        writeScalarToDisk(0, os.path.join(tmpDirName,'V_CFIE/ANGLES_FROM_FILE.txt'))

    if params_simu.MONOSTATIC_SAR==1:
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_x_hat, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_local_x_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_local_y_hat, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_local_y_hat.txt'))
        writeASCIIBlitzArrayToDisk(array(params_simu.SAR_plane_origin, 'd'), os.path.join(tmpDirName,'V_CFIE/SAR_plane_origin.txt'))
        writeScalarToDisk(params_simu.SAR_x_span, os.path.join(tmpDirName,'V_CFIE/SAR_x_span.txt'))
        writeScalarToDisk(params_simu.SAR_y_span, os.path.join(tmpDirName,'V_CFIE/SAR_y_span.txt'))
        writeScalarToDisk(params_simu.SAR_x_span_offset, os.path.join(tmpDirName,'V_CFIE/SAR_x_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_y_span_offset, os.path.join(tmpDirName,'V_CFIE/SAR_y_span_offset.txt'))
        writeScalarToDisk(params_simu.SAR_N_x_points, os.path.join(tmpDirName,'V_CFIE/SAR_N_x_points.txt'))
        writeScalarToDisk(params_simu.SAR_N_y_points, os.path.join(tmpDirName,'V_CFIE/SAR_N_y_points.txt'))


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
    if (my_id==0):
        params_simu.display()
    if (params_simu.MONOSTATIC_RCS==1) or (params_simu.MONOSTATIC_SAR==1) or (params_simu.BISTATIC==1):
        setup_excitation(params_simu, simuDirName)
    else:
        print "you should select monostatic RCS or monostatic SAR or bistatic computation, or a combination of these computations. Check the simulation settings."
        sys.exit(1)
    #MPI.Finalize()

