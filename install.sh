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
echo "  (3) openSUSE"
echo "  (4) Other Linux distribution"
echo " "
echo "Enter the correct number here: "
read DISTRIB_TYPE
echo " "

if [ $DISTRIB_TYPE = "1" ]; then
    echo "What is the version of your Fedora distribution?"
    echo "  (1) Fedora Core 17, 64 bits"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for Fedora Core 17, 64 bits: $DIR_INSTALL_SCRIPTS/FC17_64BITS_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./FC17_64BITS_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "2" ]; then
    echo "What is the version of your Ubuntu distribution?"
    echo "  (1) Ubuntu 14.04"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for Ubuntu 14.04: $DIR_INSTALL_SCRIPTS/UBUNTU_14.04_INSTALL.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./UBUNTU_14.04_INSTALL.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "3" ]; then
    echo "What is the version of your openSUSE distribution?"
    echo " "
    echo "  (1) openSUSE 13.1 (64 bit)"
    echo " "
    echo "Enter the correct number here: "
    read DISTRIB

    if [ $DISTRIB = "1" ]; then
        echo " OK, running install script for openSUSE 13.1 (64 bit): $DIR_INSTALL_SCRIPTS/openSUSE-13.1-x86_64-install.sh"
        echo " read the file if you want more info about what will be installed..."
        cd $DIR_INSTALL_SCRIPTS
        ./openSUSE-13.1-x86_64-install.sh
        cd ..
    fi

elif [ $DISTRIB_TYPE = "4" ]; then
    echo " Sorry, no install script for other distributions yet. See what is done in the other install scripts, or check the guide.pdf."
fi

