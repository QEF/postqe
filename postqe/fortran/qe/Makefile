# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

include ${QE_TOPDIR}/make.inc

FPP = ${F90} -E -cpp

MODULES_SOURCES = \
    funct.f90 \
    recips.f90 \
    latgen.f90 \
    xc_vdW_DF.f90

PW_SOURCES = \
    vloc_of_g.f90

MODULES_FILES = $(addprefix $(QE_TOPDIR)/Modules/,${MODULES_SOURCES})
PW_FILES = $(addprefix $(QE_TOPDIR)/PW/src/,${PW_SOURCES})

WRAP_FILES = ${MODULES_FILES} ${PW_FILES}
# $(info $(WRAP_FILES))

wrap_files: ${WRAP_FILES}
	for f90_file in ${WRAP_FILES}; do \
  		export fpp_file=$${f90_file%.f90}.fpp ; \
	    echo ${FPP} $$f90_file -I$${QE_TOPDIR}/include > $$fpp_file ; \
	    ${FPP} $$f90_file -I$${QE_TOPDIR}/include > $$fpp_file ; \
	done

pw:
	make -C ${QE_TOPDIR} pw

pp:
	make -C ${QE_TOPDIR} pp

move:
	mv *.so ../
	mv pyqe.py ../

qe: pw pp

clean:
	rm -f .f2py_f2cmap
	rm -f pyqe.py
	rm -f _pyqe*
	rm -f *.o
	rm -f *.so
	rm -f *.f90
	rm -f ${QE_TOPDIR}/Modules/*.fpp

f90wrappers:
	f90wrap -k kind_map -m pyqe --skip=neqq $(QE_TOPDIR)/*/*.fpp

pyqe_module:
	f2py-f90wrap -c -m _pyqe f90wrap_*.f90 \
		$(QE_TOPDIR)/Modules/funct.o \
		-I$(QE_TOPDIR)/Modules/ \
		-L$(QE_TOPDIR)/Modules -l:libqemod.a \
		-L$(QE_TOPDIR)/upflib -l:libupf.a \
		-L$(QE_TOPDIR)/XClib -l:xc_lib.a \
		-L$(QE_TOPDIR)/UtilXlib -l:libutil.a \
		-L$(QE_TOPDIR)/clib -l:clib.a \
		-L$(QE_TOPDIR)/FFTXlib -l:libqefft.a \
		$(LAPACK_LIBS) $(FFT_LIBS)

pyqe: wrap_files f90wrappers pyqe_module

all: qe pyqe move
