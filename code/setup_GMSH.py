import sys, os
from PyGmsh import executeGmsh, write_geo

if __name__=='__main__':

    sys.path.append(os.path.abspath('.'))
    from simulation_parameters import *
    if params_simu.meshToMake:
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lc', c/params_simu.f * params_simu.lc_factor)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lx', params_simu.lx)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'ly', params_simu.ly)
        write_geo(params_simu.pathToTarget, params_simu.targetName, 'lz', params_simu.lz)
        executeGmsh(params_simu.pathToTarget, params_simu.targetName, 0)


