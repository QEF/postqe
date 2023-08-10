#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
A tentative collection of functions for the XSF format (XCrySDen).
See http://www.xcrysden.org/doc/XSF.html for more info.

This is essentially a Python version of the original Fortran routines
written by Tone Kokalj (file xsf.f90 in Quantum Espresso).
"""
from .constants import BOHR_RADIUS_SI


def xsf_struct(struct_info):
    """
    Creates and returns a string with XSF format for structural info (periodic systems).

    :param struct_info: a dictionary containing the structural parameters
    :return: a string formatted as the XSF format
    """
    xsf = ' CRYSTAL\n PRIMVEC\n'
    fact = BOHR_RADIUS_SI * 1e10
    a = struct_info['a'] * fact  # convert in Ang
    for at in a:
        xsf += '{:15.9f} {:15.9f} {:15.9f}\n'.format(float(at[0]), float(at[1]), float(at[2]))
    xsf += ' PRIMCOORD\n'
    xsf += ' {:4d} 1\n'.format(int(struct_info['nat']))
    for at in struct_info['atomic_positions']:
        xsf += ' {:3s}   {:15.9f} {:15.9f} {:15.9f}\n'.format(
            at['@name'], at['$'][0]*fact, at['$'][1]*fact, at['$'][2]*fact)
    return xsf


def xsf_datagrid_2d(z, nx, ny, m1, m2, x0, e1, e2, struct_info):
    """
    Creates and returns a string with XSF 2D datablock for the charge
    (or similar quantity) on a grid.
    """
    xsf = 'BEGIN_BLOCK_DATAGRID_2D\n2D_PWSCF\nDATAGRID_2D_UNKNOWN\n'
    xsf += '{:4d}  {:4d}\n'.format(nx, ny)
    fact = struct_info['alat']*BOHR_RADIUS_SI*1e10
    for x in x0:
        xsf += '{:15.9f} '.format(x*fact)
    xsf += '\n'
    for x in e1:
        xsf += '{:15.9f} '.format(x * m1 * fact)
    xsf += '\n'
    for x in e2:
        xsf += '{:15.9f} '.format(x * m2 * fact)
    xsf += '\n'
    for j in range(0, ny):
        for i in range(0, nx):
            xsf += '{:12.4E}'.format(z[i, j].real)
            if ((j * nx + i + 1) % 6) == 0:
                xsf += '\n'
    xsf += '\n'
    xsf += 'END_DATAGRID_2D\nEND_BLOCK_DATAGRID_2D'
    return xsf


def xsf_datagrid_3d(w, nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, struct_info):
    """
    Creates and returns a string with XSF 3D datablock for the charge
    (or similar quantity) on a grid.
    """
    xsf = 'BEGIN_BLOCK_DATAGRID_3D\n3D_PWSCF\nDATAGRID_3D_UNKNOWN\n'
    xsf += '{:4d}  {:4d}  {:4d}\n'.format(nx, ny, nz)
    fact = struct_info['alat'] * BOHR_RADIUS_SI * 1e10
    for x in x0:
        xsf += '{:15.9f} '.format(x * fact)
    xsf += '\n'
    for x in e1:
        xsf += '{:15.9f} '.format(x * m1 * fact)
    xsf += '\n'
    for x in e2:
        xsf += '{:15.9f} '.format(x * m2 * fact)
    xsf += '\n'
    for x in e3:
        xsf += '{:15.9f} '.format(x * m3 * fact)
    xsf += '\n'
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                xsf += '{:12.4E}'.format(w[i, j, k].real)
                if ((k * (ny*nx) + j * nx + i + 1) % 6) == 0:
                    xsf += '\n'
    xsf += '\n'
    xsf += 'END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D'
    return xsf
