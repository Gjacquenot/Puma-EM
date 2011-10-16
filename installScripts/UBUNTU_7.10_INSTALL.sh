#!/bin/bash
###############################################################
# script for Ubuntu 7.10 to compile puma-em MoM simulator. #
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
sudo apt-get update
sudo apt-get install g++ g77 libstdc++5 python-dev python-scipy python-matplotlib python-tk openmpi-bin openmpi-dev dvipng lapack3 gmsh doxygen
#installing blitz++
sudo python installBlitz++.py
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
cp $PUMA_EM_DIR/installScripts/g77_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing mpi4py. No package yet for this one...
cd $PUMA_EM_DIR/installScripts
#wget http://pypi.python.org/packages/source/m/mpi4py/mpi4py-0.6.0.tar.gz
wget http://mpi4py.googlecode.com/files/mpi4py-1.2.2.tar.gz
tar xzf mpi4py-1.2.2.tar.gz
cd mpi4py-1.2.2
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

