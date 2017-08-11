#!/bin/bash

HDF5_VERSION="1.8"
HDF5_RELEASE="1.8.19"
HDF5_MD5SUM="7f568e2464d4ab0a74d16b23956d900b"

if [ -e hdf5.tar.gz ]
then
   echo "HDF5 source code already downloaded ..."
else
   echo "Download HDF5 source code ..."
   curl -SL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$HDF5_VERSION/hdf5-$HDF5_RELEASE/src/hdf5-$HDF5_RELEASE.tar.gz" -o hdf5.tar.gz
fi

echo -n "Checks MD5SUM of source archive ... "
echo `md5sum hdf5.tar.gz` | grep --quiet "^${HDF5_MD5SUM}[[:blank:]]"
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

if [ -e build/hdf5 ] && [ -d build/hdf5 ]
then
    read -p "Remove existing HDF5 source files? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        echo "Remove existing source files ..."
        rm -Rf build/hdf5
    else
        exit
    fi
fi

echo "Extract source files from archive ..."
mkdir build/hdf5 && tar -xz -C build/hdf5 -f hdf5.tar.gz --strip-components 1
cd build/hdf5

echo "Build HDF5 library ..."
./configure --enable-fortran --enable-fortran2003
make all

