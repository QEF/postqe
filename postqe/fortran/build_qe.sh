#!/bin/bash

QE_SOURCE_MD5SUM="c8c403fd9a9e6b3049683a94a68dc7d9"

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

if [ -e build/q-e ] && [ -d build/q-e ]
then
    read -p "Remove existing QE source files? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        echo "Remove existing source files ..."
        rm -Rf build/q-e
    else
        exit
    fi
fi

echo "Extract source files from archive ..."
unzip q-e.zip -d build/
mv build/q-e-master/ build/q-e/

echo "Build Quantum Espresso modules ..."
cd build/q-e
TOPDIR=`pwd`
# HDF5_LIB="-L\$(TOPDIR)/../hdf5/fortran/src/.libs -L\$(TOPDIR)/../hdf5/hl/fortran/src/.libs -lhdf5_fortran -lhdf5hl_fortran"
# HDF5_LIB="-L/usr/lib64/openmpi/lib/ -lhdf5_fortran -lhdf5hl_fortran"

#./configure DFLAGS="-D__FFTW -D__MPI -D__HDF5" --with-hdf5=/usr/lib64/mpich/lib/ MPIF90=mpif90 CFLAGS=-fPIC FFLAGS="-g -fPIC" \
#   IFLAGS="-I\$(TOPDIR)/include -I\$(TOPDIR)/FoX/finclude -I../include/ -I\$(TOPDIR)/../hdf5/fortran/src"

#./configure DFLAGS="-D__FFTW -D__MPI -D__HDF5" --with-hdf5=/usr/lib64/openmpi/lib/
#    IFLAGS="-I\$(TOPDIR)/include -I\$(TOPDIR)/FoX/finclude -I../include/ -I/usr/lib64/openmpi/lib/"

./configure DFLAGS="-D__FFTW -D__MPI" MPIF90=mpif90 CFLAGS=-fPIC FFLAGS="-g -fPIC"

#sed -i "s,HDF5_LIB =,HDF5_LIB = $HDF5_LIB," make.inc

make all
