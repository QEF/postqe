#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tentative collection of functions for the cube Gaussian format. See http://paulbourke.net/dataformats/cube/ or
 http://gaussian.com/cubegen/ for more info.

Note that this is different from Fortran routines written by Axel Kohlmeyer (file cube.f90 in Quantum
Espresso) and still untested! Use at your own risk.
"""
################################################################################

import numpy as np
from ase import Atom
from .constants import BOHR_RADIUS_SI


def cube(W, nx, ny, nz, e1, e2, e3, struct_info):
    """
    Creates and returns a string with cube Gaussian format.
    :param struct_info: a dictionary containing the structural parameters
    :return: a string formatted as the XSF format
    """

    # TODO: this format is untested with Gaussian or other reader

    cube =  'Cubefile created from PWScf calculation\nContains the selected quantity on a FFT grid\n'
    cube +=  '{:5d} {:12.6f} {:12.6f} {:12.6f}\n'.format(int(struct_info['nat']),0.,0.,0.)
    cube +=  '{:5d}'.format(nx)
    alat = struct_info['alat']
    for ecomp in e1:
        cube += ' {:12.6f}'.format(ecomp*alat)
    cube += '\n'
    cube +=  '{:5d}'.format(ny)
    for ecomp in e2:
        cube += ' {:12.6f}'.format(ecomp*alat)
    cube += '\n'
    cube +=  '{:5d}'.format(nz)
    for ecomp in e3:
        cube += ' {:12.6f}'.format(ecomp*alat)
    cube += '\n'
    fact = BOHR_RADIUS_SI*1e10
    for at in struct_info['atomic_positions']:
        atomic_number = Atom(at['@name']).number
        charge = atomic_number   # !!! it is unclear how this is used in the Gaussian format, here it is just the atomic number again
        cube += '{:5d} {:12.6f} {:12.6f} {:12.6f} {:12.6f}' \
                '\n'.format(atomic_number, charge, at['$'][0]*fact, at['$'][1]*fact, at['$'][2]*fact)
    # Now print the real data
    count = 0
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                cube += '{:13.5E}'.format(W[i, j, k].real)
                if ((count + 1) % 6) == 0:
                    cube += '\n'
                    count = 0
                else:
                    count += 1
            cube += '\n'
            count = 0

    return cube

