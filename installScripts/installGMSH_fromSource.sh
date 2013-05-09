#!/bin/bash

# installing source GMSH -- usually more up-to-date than packaged GMSH
echo " Root password for installing GMSH... "
wget http://www.geuz.org/gmsh/src/gmsh-2.7.0-source.tgz
tar xzf gmsh-2.7.0-source.tgz
cd gmsh-2.7.0-source
mkdir build
cd build
cmake ..
echo "On how many processes do you want to build gmsh? [default: 1]"
read N_PROCESSES
: ${N_PROCESSES:="1"}
make -j $N_PROCESSES
su -c 'make install'
cd ../..
rm -rf gmsh-2.7.0-source*

