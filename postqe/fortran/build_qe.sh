#!/bin/bash

QE_SOURCE_MD5SUM="c8c403fd9a9e6b3049683a94a68dc7d9"

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

# https://github.com/QEF/q-e/archive/1b1c427c1bc604d9746de2eefd6e2a96849e56fe.zip

if [ -e q-e.zip ]
then
   echo "QE source code already downloaded ..."
else
   echo "Download QE source code ..."
   curl -SL "https://github.com/QEF/q-e/archive/master.zip" -o q-e.zip
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
    read -p "Remove existing QE source files? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        echo "Remove existing source files ..."
        rm -Rf q-e
    else
        exit
    fi
fi

echo "Extract source files from archive ..."
unzip q-e.zip
mv q-e-master/ q-e/

echo "Build Quantum Espresso modules ..."
cd q-e
TOPDIR=`pwd`

#./configure DFLAGS="-D__FFTW -D__MPI -D__HDF5" CFLAGS=-fPIC FFLAGS="-g -fPIC" \
#    IFLAGS="-I\$(TOPDIR)/include -I\$(TOPDIR)/FoX/finclude -I../include/ -I/usr/lib64/gfortran/modules/openmpi/"

# Manually modify HDF5_LIB parameter into generated make.inc
#HDF5_LIB="-L/usr/lib64/openmpi/lib/ -lhdf5_fortran -lhdf5hl_fortran"
#sed -i "s,HDF5_LIB =,HDF5_LIB = $HDF5_LIB," make.inc

./configure DFLAGS="-D__FFTW" MPIF90=gfortran CFLAGS=-fPIC FFLAGS="-g -fPIC"
make pw
