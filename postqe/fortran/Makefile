#
# Makefile for building pyqe modules for postqe.
#
# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
QE_DIR = `pwd`/build/q-e
#PY3_LAPACK_LIB = `python -c "import numpy.linalg.lapack_lite; print(numpy.linalg.lapack_lite.__file__, end='')"`
PY3_LAPACK_LIB = -llapack
.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	gfortran -c -fPIC $<

LIBFILES=\
error_handler.o\
generate_gg_list.o\
struct_fact.o

qe:
	./build_qe.sh

copy:
	cp -pur $(QE_DIR)/POSTQE_LIBS  ./

libpostqe: $(LIBFILES)
	ar -r ./POSTQE_LIBS/$@.a $?

pyqe: wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90
	f2py3 -c wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90 -m pyqe \
	-I./POSTQE_LIBS -L./POSTQE_LIBS/ -lqemod -lqefft -lutil \
	./POSTQE_LIBS/*.o ./POSTQE_LIBS/libpw.a $(PY3_LAPACK_LIB) ./POSTQE_LIBS/libpostqe.a

module: copy libpostqe pyqe move

all: qe module

move :
	mv *.so ../

clean:
	rm -Rf ./POSTQE_LIBS
	rm -Rf ./build/
	rm -Rf ./clib/
	rm -Rf ./flib/
	rm -Rf ./Modules/
	rm -f *.o
	touch *.f90
	rm -f *.so

touch-dummy:
	$(dummy-variable)
