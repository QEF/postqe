#!/bin/bash
#
# Script for download and build a local copy Quantum ESPRESSO suite, with
# Position Independent Code shared relocatable libraries.
#
# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

# QE repository settings
QE_REPOSITORY_URL="https://github.com/QEF/q-e"
QE_SOURCE_COMMIT="8693f23037bdee531470bcbef1db1f991ab503eb"
QE_SOURCE_MD5SUM="433df3b608b7bde214ef6ad1a4966b1e"

replace=false
rebuild=false

### Argument parsing ###
usage() {
    echo "usage: $0 [-h] [--replace][--rebuild]"
    echo
    echo "optional arguments:"
    echo "  -h, --help  show this help message and exit"
    echo "  --replace   replace an existing source and build"
    echo "  --rebuild   rebuild an existing source"
}

while [ "$1" != "" ]; do
    case $1 in
        --replace )     shift
                        replace=true
                        ;;
        --rebuild )     shift
                        rebuild=true
                        ;;
        -h | --help )   usage
                        exit
                        ;;
        * )             usage
                        exit 1
    esac
    shift
done

echo "+++ Build a Quantum Espresso installation for postqe +++"

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

if [ ! -e q-e.zip ] || [ "$replace" = true ]
then
    echo "Download QE source code ..."
    curl -SL $QE_REPOSITORY_URL/archive/$QE_SOURCE_COMMIT.zip -o q-e.zip
    rebuild=true
else
    echo "QE source code already downloaded ..."
fi

echo -n "Checks MD5SUM of source archive ... "
echo `md5sum q-e.zip` | grep --quiet "^${QE_SOURCE_MD5SUM}[[:blank:]]"
if [ $? -eq 0 ]
then
    echo "OK"
else
    echo "wrong checksum!! exit ..."
    exit 1
fi

if [ ! "$rebuild" = true ]
then
    echo "Nothing to build, exit ..."
    exit
fi

if [ -d q-e ]
then
    echo "Remove existing source files ..."
    rm -Rf q-e
fi

echo "Extract source files from archive ..."
unzip q-e.zip && mv q-e-$QE_SOURCE_COMMIT/ q-e/

echo "Build Quantum Espresso modules ..."
cd q-e
./configure DFLAGS="-D__FFTW" MPIF90=gfortran CFLAGS=-fPIC FFLAGS="-g -fPIC"
make pw
