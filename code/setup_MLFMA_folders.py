import os, argparse
from mpi4py import MPI

if __name__=='__main__':
    my_id = MPI.COMM_WORLD.Get_rank()
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--simudir')
    cmdline = parser.parse_args()
    simuDirName = cmdline.simudir

    if (my_id==0):
        if 'result' not in os.listdir(simuDirName):
            os.mkdir(os.path.join(simuDirName, 'result'))
    # creation of the directories
    tmpDirName = os.path.join(simuDirName, 'tmp' + str(my_id))
    os.mkdir( tmpDirName )
    os.mkdir( os.path.join(tmpDirName,'Z_tmp') )
    os.mkdir( os.path.join(tmpDirName,'Z_near') )
    os.mkdir( os.path.join(tmpDirName,'Mg_LeftFrob') )
    os.mkdir( os.path.join(tmpDirName,'mesh') )
    os.mkdir( os.path.join(tmpDirName,'octtree_data') )
    os.mkdir( os.path.join(tmpDirName,'V_CFIE') )
    os.mkdir( os.path.join(tmpDirName,'ZI') )
    os.mkdir( os.path.join(tmpDirName,'iterative_data') )
    os.mkdir( os.path.join(tmpDirName,'pickle') )
