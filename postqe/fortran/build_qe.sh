#!/bin/bash
#
# Script for download and build a local copy Quantum ESPRESSO suite, with
# Position Independent Code shared relocatable libraries.
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

QE_SOURCE_COMMIT="585716d6435014dfc493abcd08294cd9f02eee30"
QE_SOURCE_MD5SUM="05dede6dd93e4d69109062aa429161c9"

# Move into build/ directory
if [ -d build ]
then
    echo "Directory build/ exists ..."
else
    mkdir build
    if [ $? -eq 1 ]
    then
       exit 1
    fi
fi
cd build

if [ -e q-e.zip ]
then
   echo "QE source code already downloaded ..."
else
   echo "Download QE source code ..."
   curl -SL https://github.com/QEF/q-e/archive/$QE_SOURCE_COMMIT.zip -o q-e.zip
fi

echo -n "Checks MD5SUM of source archive ... "
echo `md5sum q-e.zip` | grep --quiet "^${QE_SOURCE_MD5SUM}[[:blank:]]"
if [ $? -eq 0 ]
then
   echo "OK"
else
   echo "no match!! exit ..."
   exit 1
fi

if [ -e q-e ] && [ -d q-e ]
then
    # Checks the optional argument
    if [[ $1 =~ ^[YyNn]$ ]]
    then
        REPLACE=$1
    else
        # No or wrong argument: arks for a choice
        read -p "Replace QE source files and rebuild? [Y/N] " -n 1 -r
        echo
        REPLACE=$REPLY
    fi

    if [[ $REPLACE =~ ^[Yy]$ ]]
    then
        echo "Remove existing source files ..."
        rm -Rf q-e
    else
        exit 0
    fi
fi

echo "Extract source files from archive ..."
unzip q-e.zip && mv q-e-$QE_SOURCE_COMMIT/ q-e/

echo "Build Quantum Espresso modules ..."
cd q-e
./configure DFLAGS="-D__FFTW" MPIF90=gfortran CFLAGS=-fPIC FFLAGS="-g -fPIC"
make pw
