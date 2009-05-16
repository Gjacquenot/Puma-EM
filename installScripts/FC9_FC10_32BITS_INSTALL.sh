#/bin/sh
###############################################################
# script for fedora 9-10 to compile puma-em MoM simulator. #
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
su -c 'yum -y install python-devel gcc-c++ libgfortran gcc-gfortran libstdc++-devel compat-libstdc++-33 openmpi openmpi-libs openmpi-devel scipy numpy python-matplotlib-tk python-matplotlib cvs autoconf automake sysconftool gettext libtool'
# installing binary GMSH -- usually more up-to-date than packaged GMSH
echo " Root password for installing GMSH... "
su -c 'python installGMSH.py'
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
cp $PUMA_EM_DIR/installScripts/gfortran_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing development version of blitz++: released version is too old for gcc >= 4.3.0
# packages to be installed prior to compiling blitz: cvs autoconf automake sysconftool gettext
cd $PUMA_EM_DIR/installScripts
export CVSROOT=:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz
echo " next password is empty. Just press ENTER to continue"
cvs login
cvs -z3 checkout blitz
cd blitz
autoreconf -vif
./configure
make lib
echo " Root password for installing blitz++ library... "
su -c 'make install'
# scipy-weave uses the old blitz++, so we need to replace them
echo " Root password for replacing scipy buggy blitz++ library by the newest version of blitz++... "
su -c 'cp -r ./blitz /usr/lib/python2.5/site-packages/scipy/weave/blitz'
# installing mpi4py. No rpm yet for this one...
cd $PUMA_EM_DIR/installScripts
wget http://pypi.python.org/packages/source/m/mpi4py/mpi4py-0.6.0.tar.gz
tar xzf mpi4py-0.6.0.tar.gz
cd mpi4py-0.6.0
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
