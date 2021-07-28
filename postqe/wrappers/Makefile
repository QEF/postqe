# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

include ${QE_TOPDIR}/make.inc

FPP = ${F90} -E -cpp
# FFLAGS="-fPIC -fallow-argument-mismatch"

UtilXlib_sources = mp_base.f90
WRAP_FILES = $(addprefix $(QE_TOPDIR)/UtilXlib/,${UtilXlib_sources})

wrap_files: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do \
	    echo ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
	    ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
	done

pw:
	make -C ${QE_TOPDIR} pw

pp:
	make -C ${QE_TOPDIR} pp

move:
	mv *.so ../

clean:
	rm -f .f2py_f2cmap
	rm -f wrapper_module.py
	rm -f *.o
	rm -f *.so
	rm -f *.f90
	rm -f *.fpp

f90wrappers:
	f90wrap -k kind_map -m wrapper_module *.fpp

py3module:
	f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 -L
