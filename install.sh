#!/bin/bash

echo " "
echo "Script for compiling dependencies and Puma-em Modules"
echo "====================================================="
echo " "
echo "You have to run this script on each machine/node that you intend to use."
echo " "
echo "What is your Linux distribution?"
echo " (1) Fedora Core 8, 32 bits"
echo " (2) Fedora Core 8, 64 bits"
echo " (3) Fedora Core 9, 32 bits"
echo " (4) Fedora Core 9, 64 bits"
echo " (5) Fedora Core 10, 32 bits"
echo " (6) Fedora Core 10, 64 bits"
echo " (7) Ubuntu 7.10, 32 or 64 bits"
echo " (8) Ubuntu 8.04, 32 or 64 bits"
echo " (9) Ubuntu 8.10, 32 or 64 bits"
echo " (10) Ubuntu 9.04, 32 or 64 bits"
echo " (11) OpenSuse 11.0, 32 bits"
echo " (12) OpenSuse 11.0, 64 bits"
echo " (13) OpenSuse 11.1, 64 bits"
echo " (14) CentOS 5.2, 64 bits"
echo " (15) Other Linux distribution"
echo " "
echo "Enter the correct number here: "
read DISTRIB
DIR_INSTALL_SCRIPTS="./installScripts"

if [ $DISTRIB = "1" ]; then
    echo " OK, running install script for Fedora Core 8, 32 bits: $DIR_INSTALL_SCRIPTS/FC8_32BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC8_32BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "2" ]; then
    echo " OK, running install script for Fedora Core 8, 64 bits: $DIR_INSTALL_SCRIPTS/FC8_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC8_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "3" ]; then
    echo " OK, running install script for Fedora Core 9, 32 bits: $DIR_INSTALL_SCRIPTS/FC9_FC10_32BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC9_FC10_32BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "4" ]; then
    echo " OK, running install script for Fedora Core 9, 64 bits: $DIR_INSTALL_SCRIPTS/FC9_FC10_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC9_FC10_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "5" ]; then
    echo " OK, running install script for Fedora Core 10, 32 bits: $DIR_INSTALL_SCRIPTS/FC9_FC10_32BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC9_FC10_32BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "6" ]; then
    echo " OK, running install script for Fedora Core 10, 64 bits: $DIR_INSTALL_SCRIPTS/FC9_FC10_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./FC9_FC10_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "7" ]; then
    echo " OK, running install script for Ubuntu 7.10: $DIR_INSTALL_SCRIPTS/UBUNTU_7.10_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./UBUNTU_7.10_INSTALL.sh
    cd ..
elif [ $DISTRIB = "8" ]; then
    echo " OK, running install script for Ubuntu 8.04: $DIR_INSTALL_SCRIPTS/UBUNTU_8.04_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./UBUNTU_8.04_INSTALL.sh
    cd ..
elif [ $DISTRIB = "9" ]; then
    echo " OK, running install script for Ubuntu 8.10: $DIR_INSTALL_SCRIPTS/UBUNTU_8.10_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./UBUNTU_8.10_INSTALL.sh
    cd ..
elif [ $DISTRIB = "10" ]; then
    echo " OK, running install script for Ubuntu 9.04: $DIR_INSTALL_SCRIPTS/UBUNTU_9.04_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./UBUNTU_9.04_INSTALL.sh
    cd ..
elif [ $DISTRIB = "11" ]; then
    echo " OK, running install script for OpenSuse 11.0: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.0_32BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./OPENSUSE_11.0_32BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "12" ]; then
    echo " OK, running install script for OpenSuse 11.0: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.0_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./OPENSUSE_11.0_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "13" ]; then
    echo " OK, running install script for OpenSuse 11.1: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.1_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./OPENSUSE_11.1_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "14" ]; then
    echo " OK, running install script for CentOS 5.2: $DIR_INSTALL_SCRIPTS/CENTOS5_64BITS_INSTALL.sh"
    echo " read the file if you want more info about what will be installed..."
    cd $DIR_INSTALL_SCRIPTS
    ./CENTOS5_64BITS_INSTALL.sh
    cd ..
elif [ $DISTRIB = "15" ]; then
    echo " Sorry, no install script for Other distributions yet. See what is done in the other install scripts."
fi

