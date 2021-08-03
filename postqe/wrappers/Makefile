# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

include ${QE_TOPDIR}/make.inc

FPP = ${F90} -E -cpp
# FFLAGS="-fPIC -fallow-argument-mismatch"

UtilXlib_sources =  clib_wrappers.f90 \
	clocks_handler.f90 \
	data_buffer.f90 \
	device_helper.f90 \
	divide.f90 \
	error_handler.f90 \
	export_gstart_2_solvers.f90 \
	find_free_unit.f90 \
	fletcher32_mod.f90 \
	hash.f90 \
	mem_counter.f90 \
	mp_bands_util.f90 \
	mp_base.f90 \
	mp_base_gpu.f90 \
	mp.f90 \
	nvtx_wrapper.f90 \
	parallel_include.f90 \
	print_mem.f90 \
	set_mpi_comm_4_solvers.f90 \
	thread_util.f90 \
	util_param.f90


WRAP_FILES = $(addprefix $(QE_TOPDIR)/UtilXlib/,${UtilXlib_sources})

wrap_files: ${WRAP_FILES}
	for f90_file in ${WRAP_FILES}; do \
  		fpp_file=${f90_file%.f90}.fpp ; \
	    echo ${FPP} $$f90_file > $$fpp_file ; \
	    ${FPP} $$f90_file > $$fpp_file ; \
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
	rm -f *.fpp

f90wrappers:
	f90wrap -k kind_map -m wrapper_module $(QE_TOPDIR)/*/*.fpp

py3module:
	f2py-f90wrap -c -m _wrapper_module f90wrap_*.f90 \
		-I$(QE_TOPDIR)/UtilXlib/
