#!/bin/bash
###############################################################
# script for Ubuntu 8.10 to compile puma-em MoM simulator. #
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
sudo apt-get install g++ gfortran gmsh autoconf libtool python-dev python-scipy python-matplotlib python-tk openmpi-bin libopenmpi-dev dvipng cvs doxygen ssh
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
sudo cp -r ./blitz /usr/lib/python2.6/dist-packages/scipy/weave/blitz/
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
# now correcting a little problem on Ubuntu 9.04:
echo " "
echo " "
echo "!!WARNING!!"
echo "Ubuntu 9.04 has a new default behavior for ssh: the local user cannot ssh into the localhost (this machine) without passord anymore, whic causes Open-MPI (and hence Puma-EM) to crash. I know: lame. However, here is the fix in 2 steps, copied from http://www.cs.umd.edu/~arun/misc/ssh.html"
echo " "
echo "  1. Firstly, generate your public/private keys using ssh-keygen

    % ssh-keygen -t rsa

You must use the -t option to specify that you are producing keys for SSHv2 using RSA. This will generate your id_rsa and id_rsa.pub in the .ssh directory in your home directory. I strongly suggest using a passphrase.

  2. copy the id_rsa.pub to the .ssh directory of the localhost you want to logon to as authorized_keys2.

    % cp ~/.ssh/id_rsa.pub ~/.ssh/authorized_keys2

This script can automatically do the above maneuvers for you. Do you want that? [y/N]"
read ANSWER

: ${ANSWER:="N"}
if [ "$ANSWER" = "y" ]; then
    echo " "
    echo " leave the following passphrase empty (just press enter to questions)"
    ssh-keygen -t rsa
    cp ~/.ssh/id_rsa.pub ~/.ssh/authorized_keys2
else
    echo " OK, we assume now that you can ssh on your localhost (this machine) without a password..."
fi

# choose the appropriate make according to your MPI library
make install_open-mpi
#make install_lam-mpi
echo " "
echo "=========================================================================="
echo "                         INSTALLATION COMPLETE! "
echo "=========================================================================="
echo " "

