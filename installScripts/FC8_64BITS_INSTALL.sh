#/bin/sh
###############################################################
# script for fedora 8 to compile puma-em MoM simulator. #
# script produced by simon.batchelor@wirelessconsultancy.com #
###############################################################

echo " "
echo "installing Puma-em"
echo "=================="
echo " "
echo " BEWARE: Some parts of this script need to be run as root... "
echo " "
echo " BEWARE: if, in the future, yum updates LAM/MPI or OpenMPI, you'll have to reinstall Puma-EM."
echo " Otherwise, mpi4py risks crashing."
echo " "
echo " press enter to continue, ctrl-C to stop "
read
rm -rf ~/.python25_compiled
echo " You will be asked for your root password so that the machine can install some programs as root"
# installing the main dependencies...
echo " Root password for installing general dependencies... "
su -c 'yum -y install python-devel gcc-c++ libgfortran gcc-gfortran libstdc++-devel compat-gcc-34-g77 compat-libstdc++-33 openmpi openmpi-libs openmpi-devel scipy numpy python-matplotlib-tk python-matplotlib blitz blitz-devel mesa-libGLU'
# relinking libg2c.so: bug in Fedora Core 8
echo " Root password for correcting a Fedora Core 8 link bug for libg2c.so "
su -c 'ln -s /usr/lib/gcc/x86_64-redhat-linux/3.4.6/libg2c.so /usr/lib/libg2c.so'
# installing binary GMSH -- usually more up-to-date than packaged GMSH
echo " Root password for installing GMSH... "
su -c 'python installGMSH.py'
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
cp $PUMA_EM_DIR/installScripts/g77_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing mpi4py. No rpm yet for this one...
cd $PUMA_EM_DIR/installScripts
#wget http://pypi.python.org/packages/source/m/mpi4py/mpi4py-0.6.0.tar.gz
wget http://mpi4py.googlecode.com/files/mpi4py-1.2.2.tar.gz
tar xzf mpi4py-1.2.2.tar.gz
cd mpi4py-1.2.2
echo " Root password for installing mpi4py... "
su -c 'python setup.py install'
# cleaning up...
cd $PUMA_EM_DIR/installScripts
echo " Root password for post installation clean... "
su -c 'make clean'
# actual Puma-em installation
cd $PUMA_EM_DIR
# choose the appropriate make according to your MPI library
make install_open-mpi
#make install_lam-mpi
echo " "
echo "=========================================================================="
echo "                         INSTALLATION COMPLETE! "
echo "=========================================================================="
echo " "

