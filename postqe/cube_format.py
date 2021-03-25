#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
A tentative collection of functions for the cube Gaussian format.
See http://paulbourke.net/dataformats/cube/ or http://gaussian.com/cubegen/ for more info.

Note that this is different from Fortran routines written by Axel Kohlmeyer
(file cube.f90 in Quantum Espresso) and still untested! Use at your own risk.
"""
from ase import Atom
from .constants import BOHR_RADIUS_SI


def cube(W, nx, ny, nz, e1, e2, e3, struct_info):
    """
    Creates and returns a string with cube Gaussian format. *struct_info* is
    a dictionary containing the structural parameters. Return a string formatted
    as the XSF format
    """

    # TODO: this format is untested with Gaussian or other reader

    text = 'Cubefile created from PWScf calculation\n' \
           'Contains the selected quantity on a FFT grid\n'
    text += '{:5d} {:12.6f} {:12.6f} {:12.6f}\n'.format(int(struct_info['nat']), 0., 0., 0.)
    text += '{:5d}'.format(nx)
    alat = struct_info['alat']
    for ecomp in e1:
        text += ' {:12.6f}'.format(ecomp * alat)
    text += '\n'
    text += '{:5d}'.format(ny)
    for ecomp in e2:
        text += ' {:12.6f}'.format(ecomp * alat)
    text += '\n'
    text += '{:5d}'.format(nz)
    for ecomp in e3:
        text += ' {:12.6f}'.format(ecomp*alat)
    text += '\n'
    fact = BOHR_RADIUS_SI*1e10
    for at in struct_info['atomic_positions']:
        atomic_number = Atom(at['@name']).number

        # !!! it is unclear how this is used in the Gaussian format,
        # here it is just the atomic number again
        charge = atomic_number
        text += '{:5d} {:12.6f} {:12.6f} {:12.6f} {:12.6f}' \
                '\n'.format(atomic_number, charge, at['$'][0] * fact,
                            at['$'][1] * fact, at['$'][2] * fact)
    # Now print the real data
    count = 0
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                text += '{:13.5E}'.format(W[i, j, k].real)
                if ((count + 1) % 6) == 0:
                    text += '\n'
                    count = 0
                else:
                    count += 1
            text += '\n'
            count = 0

    return text
