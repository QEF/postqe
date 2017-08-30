#
# Makefile for building pyqe modules for postqe.
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
QE_DIR = `pwd`/build/q-e
PY2_LAPACK_LIB = `python2 -c "import numpy.linalg.lapack_lite; print(numpy.linalg.lapack_lite.__file__)"`
PY3_LAPACK_LIB = `python3 -c "import numpy.linalg.lapack_lite; print(numpy.linalg.lapack_lite.__file__, end='')"`

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	gfortran -c -fPIC $<

LIBFILES=\
error_handler.o\
generate_gg_list.o\
struct_fact.o

qe:
	./build_qe.sh y

copy:
	mkdir -p clib/
	cp -pu $(QE_DIR)/clib/*.o clib/
	mkdir -p flib/
	cp -pu $(QE_DIR)/Modules/libqemod.a flib/
	cp -pu $(QE_DIR)/FFTXlib/libqefft.a flib/
	cp -pu $(QE_DIR)/PW/src/libpw.a flib/
	cp -pu $(QE_DIR)/UtilXlib/libutil.a flib/
	mkdir -p Modules/
	cp -pu $(QE_DIR)/Modules/funct.mod Modules/

libpostqe: $(LIBFILES)
	ar -r ./flib/$@.a $?

pyqe2: wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90
	f2py2 -c wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90 -m pyqe \
	-I./Modules/ -L./flib/ -lqemod -lqefft -lutil \
	./clib/*.o ./flib/libpw.a $(PY2_LAPACK_LIB) ./flib/libpostqe.a

pyqe3: wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90
	f2py3 -c wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90 -m pyqe \
	-I./Modules/ -L./flib/ -lqemod -lqefft -lutil \
	./clib/*.o ./flib/libpw.a $(PY3_LAPACK_LIB) ./flib/libpostqe.a

module2: copy libpostqe pyqe2 move

module3: copy libpostqe pyqe3 move

all: qe module3

move :
	mv *.so ../

clean:
	rm -Rf ./build/
	rm -Rf ./clib/
	rm -Rf ./flib/
	rm -Rf ./Modules/
	rm -f *.o
	touch *.f90
	rm -f *.so

touch-dummy:
	$(dummy-variable)
