# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

include ${QE_TOPDIR}/make.inc

FPP = ${F90} -E -cpp
# FFLAGS="-fPIC -fallow-argument-mismatch"

MODULES_SOURCES = funct.f90
#    xc_rVV10.f90 \
#    xc_vdW_DF.f90 \
#    recips.f90 \
#    latgen.f90

MODULES_FILES = $(addprefix $(QE_TOPDIR)/Modules/,${MODULES_SOURCES})

# MODULES_FILES = $(wildcard $(QE_TOPDIR)/Modules/*.f90)
UTIL_XLIB_FILES = $(wildcard $(QE_TOPDIR)/UtilXlib/*.f90)

# $(info $(MODULES_FILES))

WRAP_FILES = ${MODULES_FILES}
# WRAP_FILES = ${UTIL_XLIB_FILES}

wrap_files: ${WRAP_FILES}
	for f90_file in ${WRAP_FILES}; do \
  		export fpp_file=$${f90_file%.f90}.fpp ; \
	    echo ${FPP} $$f90_file -I$${QE_TOPDIR}/include > $$fpp_file ; \
	    ${FPP} $$f90_file -I$${QE_TOPDIR}/include > $$fpp_file ; \
	done

# wrap_files: ${WRAP_FILES}
#	for f in ${WRAP_FILES}; do \
#	    echo ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
#	    ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
#	done

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
	rm -f ${QE_TOPDIR}/Modules/*.fpp
	rm -f ${QE_TOPDIR}/UtilXlib/*.fpp

f90wrappers:
	f90wrap -k kind_map --only set_dft_from_name -m wrapper_module $(QE_TOPDIR)/*/*.fpp

py3module:
	f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 \
	    $(QE_TOPDIR)/Modules/funct.o \
	    -llapack \
		-I$(QE_TOPDIR)/Modules/ \
		-L$(QE_TOPDIR)/Modules -l:libqemod.a \
		-L$(QE_TOPDIR)/upflib -l:libupf.a \
		-L$(QE_TOPDIR)/XClib -l:xc_lib.o

py3module-old2:
	f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 \
	    $(QE_TOPDIR)/Modules/io_global.o \
	    $(QE_TOPDIR)/Modules/io_files.o \
	    $(QE_TOPDIR)/Modules/ions_base.o \
	    $(QE_TOPDIR)/Modules/mp_images.o \
	    $(QE_TOPDIR)/Modules/read_pseudo.o \
	    $(QE_TOPDIR)/upflib/uspp.o \
	    $(QE_TOPDIR)/upflib/upf_auxtools.o \
	    $(QE_TOPDIR)/upflib/upf_io.o \
	    $(QE_TOPDIR)/upflib/upf_spinorb.o \
		$(QE_TOPDIR)/upflib/atom.o \
		-I$(QE_TOPDIR)/Modules/ \
		-I$(QE_TOPDIR)/upflib/ \
		-L$(QE_TOPDIR)/upflib/ \
		-L$(QE_TOPDIR)/Modules/

old-py3module:
	f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 \
	    $(QE_TOPDIR)/Modules/*.o \
	    $(QE_TOPDIR)/UtilXlib/*.o \
	    $(QE_TOPDIR)/upflib/atom.o \
	    $(QE_TOPDIR)/upflib/uspp.o \
	    $(QE_TOPDIR)/upflib/upf_spinorb.o \
	    $(QE_TOPDIR)/MBD/libmbd-*/src/mbd*.o \
	    $(QE_TOPDIR)/dft-d3/*.o \
	    $(QE_TOPDIR)/XClib/qe_funct_*.o \
	    $(QE_TOPDIR)/PW/src/beef.o \
	    $(QE_TOPDIR)/PW/src/bp_mod.o \
	    $(QE_TOPDIR)/PW/src/pw_restart_new.o \
	    $(QE_TOPDIR)/PW/src/setup.o \
	    $(QE_TOPDIR)/PW/src/symm_base.o \
	    $(QE_TOPDIR)/PW/src/esm_common_mod.o \
		-I$(QE_TOPDIR)/Modules/ -L$(QE_TOPDIR)/Modules/
