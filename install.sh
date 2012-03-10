#!/bin/bash

DIR_INSTALL_SCRIPTS="./installScripts"

echo " "
echo "Script for compiling dependencies and Puma-em Modules"
echo "====================================================="
echo " "
echo "You have to run this script on each machine/node that you intend to use."
echo " "
echo "What is your the type of your Linux distribution?"
echo "  (1) Fedora Core"
echo "  (2) Ubuntu"
echo "  (3) OpenSuse"
echo "  (4) CentOS"
echo "  (5) Other Linux distribution"
echo " "
echo "Enter the correct number here: "
read DISTRIB_TYPE
echo " "

if [ $DISTRIB_TYPE = "1" ]; then
    echo "What is the version of your Fedora distribution?"
    echo "  (1) Fedora Core 13, 64 bits"
    echo "  (2) Fedora Core 14, 64 bits"
    echo "  (3) Fedora Core 15, 64 bits"
    echo "  (4) Fedora Core 16, 64 bits"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for Fedora Core 13, 64 bits: $DIR_INSTALL_SCRIPTS/FC13_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./FC13_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "2" ]; then
        echo " OK, running install script for Fedora Core 14, 64 bits: $DIR_INSTALL_SCRIPTS/FC14_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./FC14_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "3" ]; then
        echo " OK, running install script for Fedora Core 15, 64 bits: $DIR_INSTALL_SCRIPTS/FC15_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./FC15_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "4" ]; then
        echo " OK, running install script for Fedora Core 16, 64 bits: $DIR_INSTALL_SCRIPTS/FC16_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./FC16_64BITS_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "2" ]; then
    echo "What is the version of your Ubuntu distribution?"
    echo "  (1) Ubuntu 10.04"
    echo "  (2) Ubuntu 10.10"
    echo "  (3) Ubuntu 11.04"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for Ubuntu 10.04: $DIR_INSTALL_SCRIPTS/UBUNTU_10.04_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./UBUNTU_10.04_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "2" ]; then
        echo " OK, running install script for Ubuntu 10.10: $DIR_INSTALL_SCRIPTS/UBUNTU_10.10_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./UBUNTU_10.10_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "3" ]; then
        echo " OK, running install script for Ubuntu 11.04: $DIR_INSTALL_SCRIPTS/UBUNTU_11.04_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./UBUNTU_11.04_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "3" ]; then
    echo "What is the version of your OpenSuse distribution?"
    echo "  (1) OpenSuse 11.0, 64 bits"
    echo "  (2) OpenSuse 11.1, 64 bits"
    echo "  (3) OpenSuse 11.2, 64 bits"
    echo "  (4) OpenSuse 11.3, 64 bits"
    echo "  (5) OpenSuse 11.4, 64 bits"
    echo "  (6) OpenSuse 12.1, 64 bits"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB


    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for OpenSuse 11.0: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.0_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_11.0_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "2" ]; then
        echo " OK, running install script for OpenSuse 11.1: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.1_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_11.1_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "3" ]; then
        echo " OK, running install script for OpenSuse 11.2: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.2_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_11.2_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "4" ]; then
        echo " OK, running install script for OpenSuse 11.3: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.3_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_11.3_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "5" ]; then
        echo " OK, running install script for OpenSuse 11.4: $DIR_INSTALL_SCRIPTS/OPENSUSE_11.4_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_11.4_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "6" ]; then
        echo " OK, running install script for OpenSuse 12.1: $DIR_INSTALL_SCRIPTS/OPENSUSE_12.1_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./OPENSUSE_12.1_64BITS_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "4" ]; then
    echo "What is the version of your CentOS distribution?"
    echo "  (1) CentOS 5.2"
    echo "  (2) CentOS 5.5"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for CentOS 5.2: $DIR_INSTALL_SCRIPTS/CENTOS5_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./CENTOS5_64BITS_INSTALL.sh
        cd ..
    elif [ $DISTRIB = "2" ]; then
        echo " OK, running install script for CentOS 5.5: $DIR_INSTALL_SCRIPTS/CENTOS5.5_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./CENTOS5.5_64BITS_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "5" ]; then
    echo " Sorry, no install script for other distributions yet. See what is done in the other install scripts, or check the guide.pdf."
fi

