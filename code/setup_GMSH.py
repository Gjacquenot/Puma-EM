import sys, os, argparse
from PyGmsh import executeGmsh, write_geo, isGeoFileThere

if __name__=='__main__':
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

    if 'geo' not in os.listdir(simuDirName):
        os.mkdir(os.path.join(simuDirName, 'geo'))
    geoDirName = os.path.join(simuDirName, 'geo')

    sys.path.append(os.path.abspath('.'))
    exec 'from ' + simuParams + ' import *'
    filename = 'GMSHcommand.sh'
    if params_simu.meshToMake:
        os.system("cp " + os.path.join(params_simu.pathToTarget, params_simu.targetName + '.geo') + ' ' + geoDirName)
        write_geo(geoDirName, params_simu.targetName, 'lc', c/params_simu.f * params_simu.lc_factor)
        write_geo(geoDirName, params_simu.targetName, 'lx', params_simu.lx)
        write_geo(geoDirName, params_simu.targetName, 'ly', params_simu.ly)
        write_geo(geoDirName, params_simu.targetName, 'lz', params_simu.lz)
        isGeoFileThere(geoDirName, params_simu.targetName)
        f = open(filename, 'w')
        GMSH_command = 'rm -f ' + os.path.join(geoDirName, params_simu.targetName) + '.msh* \n'
        GMSH_command += 'gmsh -2 -algo del2d -rand 1e-06 ' + os.path.join(geoDirName, params_simu.targetName) + '.geo' + ' -string "General.ExpertMode=1;"\n'
        f.write(GMSH_command)
        f.write("exit 0" + "\n")
    else:
        f = open(filename, 'w')
        f.write("exit 0" + "\n")
    f.close()
    os.system("chmod u+x " + filename)
        

