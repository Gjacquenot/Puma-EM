#!/bin/bash
###############################################################
# script for Ubuntu 14.04 to compile puma-em MoM simulator.   #
###############################################################

echo " "
echo "installing Puma-em"
echo "=================="
echo " "
echo " BEWARE: Some parts of this script need to be run as sudo... "
echo " "
echo " press enter to continue, ctrl-C to stop"
read
echo " You will be asked for your root password so that the machine can install some programs as root"
# installing the main dependencies...
echo " sudo password for installing main dependencies... "
sudo apt-get update
sudo apt-get install g++ gfortran gmsh autoconf libtool python-dev python-scipy python-matplotlib python-mpi4py python-tk openmpi-bin libopenmpi-dev dvipng cvs automake libblitz0-dev libblitz0ldbl
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
#cp $PUMA_EM_DIR/installScripts/gfortran_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing blitz++
#cd $PUMA_EM_DIR/installScripts
#export CVSROOT=:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz
#echo " next password is empty. Just press ENTER to continue"
#cvs login
#cvs -z3 checkout blitz
#cd blitz
#autoreconf -vif
#./configure
#make
#sudo make install
cd $PUMA_EM_DIR/installScripts
sudo make clean
# actual Puma-em installation
cd $PUMA_EM_DIR
echo " "
echo " "
make libs
echo " "
echo "=========================================================================="
echo "                         INSTALLATION COMPLETE! "
echo "=========================================================================="
echo " "

