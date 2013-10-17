#!/bin/bash
###############################################################
# script for Ubuntu 12.10 to compile puma-em MoM simulator.   #
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
sudo apt-get update
sudo apt-get install g++ gfortran gmsh autoconf libtool python-dev python-scipy python-matplotlib python-tk openmpi-bin libopenmpi-dev dvipng cvs doxygen
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
make
sudo make install
# scipy-weave uses the old blitz++, so we need to replace them
sudo cp -r ./blitz /usr/lib/python2.7/dist-packages/scipy/weave/blitz/
# installing mpi4py. No package yet for this one...
cd $PUMA_EM_DIR/installScripts
wget http://mpi4py.googlecode.com/files/mpi4py-1.3.tar.gz
tar xzf mpi4py-1.3.tar.gz
cd mpi4py-1.3
sudo python setup.py install
# cleaning up...
cd $PUMA_EM_DIR/installScripts
sudo make clean
# actual Puma-em installation
cd $PUMA_EM_DIR
echo " "
echo " "
make install_open-mpi
echo " "
echo "=========================================================================="
echo "                         INSTALLATION COMPLETE! "
echo "=========================================================================="
echo " "

