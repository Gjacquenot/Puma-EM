#/bin/sh

#-----------------------------------------------------------------------
# Script to compile PUMA-EM MoM simulator for openSUSE 13.1 (64 bit).
#
# Based on script for installation on Fedora Core 17 (64 bit) by 
# Idesbald Van den Bosch (vandenbosch.idesbald@gmail.com).
#
# Author: Maksims Abalenkovs (abalenkovs@ieee.org)
# Date:   Jun 17, 2014
# Ver:    0.2
#-----------------------------------------------------------------------

echo " "
echo "Installing PUMA-EM"
echo "------------------"
echo " "
echo " NOTE: Some parts of this script need to be run as root... "
echo " "
echo " NOTE: It will be necessary to re-install PUMA-EM in case YaST/zypper updates LAM/MPI or OpenMPI."
echo "       Otherwise, mpi4py might crash."
echo " "
echo " Press Enter to continue or Ctrl+C to stop."
read
rm -rf ~/.python*_compiled
echo " You will be asked for root password in order for the machine to install programs."
#
# Install main dependencies
#
echo " Installing main software dependencies..."
echo " Enter root password, please."
#
su -c 'zypper in autoconf automake blas-devel blas-man cmake cvs doxygen fltk-devel fltk-devel-static gcc-c++ gcc-fortran gettext-runtime gettext-tools glu glu-devel hdf5 hdf5-devel hdf5-openmpi hdf5-openmpi-devel lapack-devel lapack-devel-static lapack-man libatlas3 libatlas3-devel libatlas3-sse libatlas3-sse2 libatlas3-sse2-devel libatlas3-sse3 libatlas3-sse3-devel libatlas3-sse-devel libblas3 libfltk1 libgfortran3 libgmm++-devel libjpeg8 libjpeg8-devel libpng16-16 libpng16-devel libstdc++33 libstdc++33-devel libstdc++-devel libtool Mesa-libGL-devel openmpi openmpi-devel python python-argparse python-devel python-matplotlib python-matplotlib-tk python-numpy python-pip python-scipy texlive-latex texlive-fourier texlive-multirow wget zlib-devel'
#
# Re-create soft links for OpenMPI packages
#
echo " Re-creating soft links for OpenMPI packages..."
echo " Enter root password, please."
su -c 'ln -s /usr/lib64/openmpi/bin/mpicc /usr/bin/mpicc; ln -s /usr/lib64/openmpi/bin/mpiCC /usr/bin/mpiCC; ln -s /usr/lib64/openmpi/bin/mpirun /usr/bin/mpirun; ln -s /usr/lib64/openmpi/lib/libmpi.so.1 /usr/lib64/libmpi.so.1; ln -s /usr/lib64/openmpi/lib/libopen-rte.so.1 /usr/lib64/libopen-rte.so.1; ln -s /usr/lib64/openmpi/lib/libopen-pal.so.1 /usr/lib64/libopen-pal.so.1; ln -s /usr/lib64/openmpi/lib/libmpi_cxx.so.1 /usr/lib64/libmpi_cxx.so.1;'
#
# Install GMSH from source
#
echo " Installing GMSH..."
echo " Enter root password, please."
./installGMSH_fromSource.sh
# create makefile.inc
cd ..
PUMA_EM_DIR=$PWD
#cp $PUMA_EM_DIR/installScripts/gfortran_makefile.inc $PUMA_EM_DIR/makefile.inc
# installing development version of blitz++: released version is too old for gcc >= 4.3.0
# packages to be installed prior to compiling blitz: cvs autoconf automake sysconftool gettext
echo " "
echo " Obtaining Blitz++ from CVS server..."
cd $PUMA_EM_DIR/installScripts
export CVSROOT=:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz
echo " Next password is empty. Please, press Enter to continue."
cvs login
cvs -z3 checkout blitz
echo " Configuring Blitz++ library..."
cd blitz
autoreconf -vif
./configure
echo " Making Blitz++ library..."
make lib
echo " Installing Blitz++ library... "
echo " Enter root password, please."
su -c 'make install'
#
# Replace blitz++ used by Scipy-Weave with the newest version
#
echo " Replacing old Blitz++ library by the new version of Blitz++ in Scipy-Weave..."
echo " Enter root password, please."
su -c 'cp -r ./blitz /usr/lib64/python2.7/site-packages/scipy/weave/blitz'
#
# Install mpi4py from source since no RPM exists yet
#
echo " Installing mpi4py from source..."
cd $PUMA_EM_DIR/installScripts
echo " Obtaining mpi4py source code..."
wget http://mpi4py.googlecode.com/files/mpi4py-1.3.tar.gz
echo " Uncompressing mpi4py source code archive..."
tar xzf mpi4py-1.3.tar.gz
cd mpi4py-1.3
echo " Installing mpi4py..."
echo " Enter root password, please."
su -c 'python setup.py install'
#
# Clean up
#
echo " Cleaning up the installation..."
cd $PUMA_EM_DIR/installScripts
echo " Enter root password, please."
su -c 'make clean'
#
# Proceed with the actual PUMA-EM installation
#
echo "Installing PUMA-EM..."
cd $PUMA_EM_DIR
#
# Choose appropriate make according to your MPI library
#
make install_open-mpi
# make install_lam-mpi
echo " "
echo "------------------------------------------------------------------"
echo "               Installation completed successfully.               "
echo "------------------------------------------------------------------"
echo " "
