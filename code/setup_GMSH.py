import sys, os, argparse
from PyGmsh import executeGmsh, write_geo, isGeoFileThere

if __name__=='__main__':

    sys.path.append(os.path.abspath('.'))
    from simulation_parameters import *
    filename = 'GMSHcommand.sh'
    if params_simu.meshToMake:
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lc', c/params_simu.f * params_simu.lc_factor)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lx', params_simu.lx)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'ly', params_simu.ly)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lz', params_simu.lz)
        #executeGmsh(params_simu.pathToTarget, params_simu.targetName, 0)
        isGeoFileThere(params_simu.pathToTarget, params_simu.targetName)
        f = open(filename, 'w')
        GMSH_command = 'rm -f ' + os.path.join(params_simu.pathToTarget, params_simu.targetName) + '.msh* \n'
        GMSH_command += 'gmsh -2 -algo del2d ' + os.path.join(params_simu.pathToTarget, params_simu.targetName) + '.geo' + ' -string "General.ExpertMode=1;"\n'
        f.write(GMSH_command)
        f.write("exit 0" + "\n")
    else:
        f = open(filename, 'w')
        f.write("exit 0" + "\n")
    f.close()
    os.system("chmod u+x " + filename)
        

