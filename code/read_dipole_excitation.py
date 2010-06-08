from scipy import array, zeros

def read_dipole_excitation(filename):
    # the structure of the excitation file must be as follows:
    # 1 line per dipole, as many lines as there are dipoles
    # each line has 9 columns that must be arranged as follows:
    #
    # real(J_x) imag(J_x) real(J_y) imag(J_y) real(J_z) imag(J_z) r_x r_y r_z
    #
    # where J = [J_x J_y J_z] is the dipole and r = [r_x r_y r_z] its origin.
    r_tmp, J_tmp = [], []
    excitation_file = open(filename, 'r')
    for line in excitation_file:
        elems = line.split()
        if (len(elems)==9) and (elems!=[]):
            r_tmp.append(map(float, elems[6:]))
            J_tmp.append(map(float, elems[:6]))
    excitation_file.close()
    if r_tmp == []:
        J_src, r_src = zeros((0,0), 'D'), zeros((0,0), 'd')
    else:
        r_src = array(r_tmp, 'd')
        J_tmp2 = array(J_tmp, 'd')
        J_src = array(J_tmp2[:,0:6:2], 'd') + array(J_tmp2[:,1:6:2], 'd') * 1.j
    return J_src, r_src

def read_observation_points(filename):
    # the structure of the r_obs file MUST BE AS FOLLOWS:
    # 1 line per observation point, as many lines as there are points
    # each line has 3 columns that must be arranged as follows:
    #
    # r_obs_x r_obs_y r_obs_z

    r_tmp = []
    observation_file = open(filename, 'r')
    for line in observation_file:
        elems = line.split()
        if (len(elems)==3) and (elems!=[]):
            r_tmp.append(map(float, elems[0:]))
    observation_file.close()
    if r_tmp == []:
        r_obs = zeros((0, 0), 'd')
    else:
        r_obs = array(r_tmp, 'd')
    return r_obs
    
