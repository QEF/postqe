# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
move:
	mv *.so ../

clean:
	rm -f f90utils*
	rm -f *.o
	rm -f *.so

f90utils_module:
	f2py3 -c *.f90 -m f90utils

all: f90utils_module move
