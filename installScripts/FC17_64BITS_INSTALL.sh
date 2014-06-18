#/bin/sh
###############################################################
# script for fedora 17 to compile puma-em MoM simulator.      #
# script produced by vandenbosch.idesbald@gmail.com)          #
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
rm -rf ~/.python*_compiled
echo " You will be asked for your root password so that the machine can install some programs as root"
# installing the main dependencies...
echo " Root password for installing general dependencies... "
su -c 'yum -y install python-devel gcc-c++ libgfortran gcc-gfortran libstdc++-devel compat-libstdc++-33 openmpi openmpi-devel scipy numpy python-matplotlib-tk python-matplotlib cvs autoconf automake sysconftool libtool.x86_64 gettext zlib.i686 mesa-libGLU mesa-libGLU-devel compat-libstdc++-33.i686 wget cmake fltk-devel lapack lapack-devel blas blas-devel libjpeg-devel libpng-devel'
# repairing shite introduced in FC13 in the Open-MPI packages...
su -c 'ln -s /usr/lib64/openmpi/bin/mpicc /usr/bin/mpicc; ln -s /usr/lib64/openmpi/bin/mpiCC /usr/bin/mpiCC; ln -s /usr/lib64/openmpi/bin/mpirun /usr/bin/mpirun; ln -s /usr/lib64/openmpi/lib/libmpi.so.1 /usr/lib64/libmpi.so.1; ln -s /usr/lib64/openmpi/lib/libopen-rte.so.1 /usr/lib64/libopen-rte.so.1; ln -s /usr/lib64/openmpi/lib/libopen-pal.so.1 /usr/lib64/libopen-pal.so.1; ln -s /usr/lib64/openmpi/lib/libmpi_cxx.so.1 /usr/lib64/libmpi_cxx.so.1;'
# installing GMSH from source
./installGMSH_fromSource.sh
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
#cp $PUMA_EM_DIR/installScripts/gfortran_makefile.inc $PUMA_EM_DIR/makefile.inc
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
su -c 'cp -r ./blitz /usr/lib64/python2.7/site-packages/scipy/weave/blitz'
# installing mpi4py. No rpm yet for this one...
cd $PUMA_EM_DIR/installScripts
wget http://mpi4py.googlecode.com/files/mpi4py-1.3.tar.gz
tar xzf mpi4py-1.3.tar.gz
cd mpi4py-1.3
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
