#!/bin/bash

QE_SOURCE_URL="https://github.com/QEF/q-e/archive/585716d6435014dfc493abcd08294cd9f02eee30.zip"
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
   curl -SL $QE_SOURCE_URL -o q-e.zip
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
mv q-e-585716d6435014dfc493abcd08294cd9f02eee30/ q-e/

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
