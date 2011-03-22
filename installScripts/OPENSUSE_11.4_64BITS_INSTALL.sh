#!/bin/bash
###############################################################
# script for OpenSuse 11.2 to compile puma-em MoM simulator. #
# script produced by simon.batchelor@wirelessconsultancy.com #
###############################################################

echo " "
echo "installing Puma-em"
echo "=================="
echo " "
echo " BEWARE: Some parts of this script need to be run as sudo... "
echo " "
echo " BEWARE: if, in the future, apt updates LAM/MPI or OpenMPI, you'll have to reinstall Puma-EM."
echo " Otherwise, mpi4py risks crashing."
echo " "
echo " press enter to continue, ctrl-C to stop"
read
rm -rf ~/.python*_compile
echo " You will be asked for your root password so that the machine can install some programs as root"
# installing the main dependencies...
echo " sudo password for installing main dependencies... "
sudo zypper ar -c http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
sudo zypper refresh
sudo zypper install gcc-c++ libstdc++33-32bit Mesa-32bit gcc-fortran autoconf automake make libtool python-devel python-scipy python-matplotlib python-tk openmpi openmpi-devel cvs doxygen
# modifying the PATH variables
mpi-selector --set openmpi-1.3
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib:
export PATH=/usr/lib64/mpi/gcc/openmpi/bin:$PATH
# installing GMSH
VERSION="2.5.0"
wget http://www.geuz.org/gmsh/bin/Linux/gmsh-$VERSION-Linux.tgz
tar xzf gmsh-$VERSION-Linux.tgz
sudo mv gmsh-$VERSION-Linux/bin/gmsh /usr/bin
rm -rf gmsh-$VERSION-Linux.tgz
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
cp $PUMA_EM_DIR/installScripts/gfortran_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing blitz++
cd $PUMA_EM_DIR/installScripts
export CVSROOT=:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz
echo " next password is empty. Just press ENTER to continue"
cvs login
cvs -z3 checkout blitz
cd blitz
autoreconf -vif
./configure
make lib
sudo make install
# scipy-weave uses the old blitz++, so we need to replace them
sudo cp -r ./blitz /usr/lib64/python2.7/site-packages/scipy/weave/blitz/
# installing mpi4py. No package yet for this one...
cd $PUMA_EM_DIR/installScripts
wget http://pypi.python.org/packages/source/m/mpi4py/mpi4py-0.6.0.tar.gz
tar xzf mpi4py-0.6.0.tar.gz
cd mpi4py-0.6.0
python setup.py build
sudo python setup.py install
# cleaning up...
cd $PUMA_EM_DIR/installScripts
sudo make clean
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

