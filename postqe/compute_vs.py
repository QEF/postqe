#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import numpy as np

from .constants import pi
from .setlocal import wrap_setlocal

# f2py created modules
from . import pyqe


# Compute the volume from a1, a2, a3 vectors in direct space.  
def compute_volume(at1, at2, at3):
    """
    This function computes the cell volume per atom given the a1, a2, a3 vectors.
    As in QE.
    """
    V = abs(at1[0] * at2[1] * at3[2] + at1[1] * at2[2] * at3[0] + at1[2] * at2[0] * at3[1] -
            at3[0] * at2[1] * at1[2] - at3[1] * at2[2] * at1[0] - at3[2] * at2[0] * at3[1])
    return V 


def compute_G(b, nr):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr. G are the vectors in the reciprocal lattice vector.
    b[0], b[1], b[2] are the reciprocal cell base vectors
    """
    gx, gy, gz = [
        np.fromiter((i - nr[k] if i >= nr[k] // 2 else i for i in range(nr[k])), dtype=np.int64)
        for k in range(3)
    ]
    g = np.fromfunction(
        function=lambda x, y, z: np.array((gx[x], gy[y], gz[z]), dtype=np.float64),
        shape=nr, dtype=np.int32
    )
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector
    return G


def compute_G_squared(b, nr, ecutrho, alat):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing G^2 at each point, G is the 
    corresponding reciprocal lattice vector. Also apply a proper cut-off 
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16
    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector

    def g_square(x, y, z):
        g2 = np.linalg.norm(G, axis=3) ** 2  # compute the G^2
        g2[(g2 == 0.0) | (g2 > ecutm)] = bignum  # dummy high value so that n(G)/G^2 is 0
        return g2

    G2 = np.fromfunction(function=g_square, shape=nr, dtype=int)
    return G2


def compute_Gs(b, nr, ecutrho, alat):
    """
    This function computes both a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr and the G^2 moduli. G are the vectors in the
    reciprocal space. Also apply a proper cut-off for G^2
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16
    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector

    def g_square(x, y, z):
        g2 = np.linalg.norm(G, axis=3) ** 2  # compute the G^2
        g2[(g2 == 0.0) | (g2 > ecutm)] = bignum  # dummy high value so that n(G)/G^2 is 0
        return g2

    G2 = np.fromfunction(function=g_square, shape=nr, dtype=int)
    return G, G2


def compute_v_bare(ecutrho, alat, at1, at2, at3, nr, atomic_positions, species, pseudodir):
    """
    This function computes the bare potential. It calls the wrapper function
    wrap_setlocal from the setlocal python module which is an interface to
    call the proper fortran functions.
    """
    v_F = wrap_setlocal(alat, at1, at2, at3, nr[0], nr[1], nr[2], atomic_positions,
                        species, 2.0*ecutrho, pseudodir)
    return v_F


def get_v_h_from_hdf5(filename, nr, dataset='rhotot_g'):
    """
    computes the Hartree potential, reading the data from the hdf5 charge-density file.
    :param filename: name of the hdf5 containing the charge density in reciprocal space 
    :param nr: three-ple containing nr1, nr2, nr3 integers defining the fft mesh 
    :param dataset: the name of the dataset containing the charge density
    :return: returns the value of the hartree potential defined in a nr1,nr2,nr3 mesh
    """
    import h5py
    nrrr = nr[0]*nr[1]*nr[2]
    with h5py.File(filename) as h5f:
        h5d = h5f[dataset]
        ngm = h5f.attrs.get('ngm_g')
        rho_g = np.array(h5d).reshape([ngm, 2])
        h5mill = h5f['MillerIndices']
        b1 = h5mill.attrs.get('bg1')
        b2 = h5mill.attrs.get('bg2')
        b3 = h5mill.attrs.get('bg3')
        aux = np.zeros(list(nr), dtype=np.complex128)
        my_zip = zip(h5mill, rho_g)
        next(my_zip)
        for el in my_zip:
            m = el[0]
            rhog = el[1].dot((1.e0, 1.j))
            g = m.dot((b1, b2, b3))
            gg = g.dot(g)
            try:
                aux[m[0], m[1], m[2]] = rhog/gg
            except IndexError:
                pass

    return np.fft.ifftn(aux).real*8.e0*np.pi*nrrr


def compute_v_h(charge, ecutrho, alat, b):
    """
    This function computes the hartree potential from the charge and
    the cut off energy "ecutrho" on the charge. The charge is a numpy matrix nr1*nr2*nr3.

    """    
    # First compute the FFT of the charge          
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape
        
    # Compute G^2 values for the mesh, applying the proper cutoff from ecutrho
    # and alat
    G2 = compute_G_squared(b, nr, ecutrho, alat)
    
    conv_fact = 2.0 / pi * alat**2
    v = np.fft.ifftn(fft_charge / G2) * conv_fact

    return v.real


def compute_v_xc(charge, charge_core, functional):
    """
    This function computes the exchange-correlation potential from the charge and
    the type of functional given in input. The charge is a numpy matrix nr1*nr2*nr3.
    The functional is a string identifying the functional as in QE convention.
    """
    vanishing_charge = 1.0E-10
    nr = charge.shape
    v = np.zeros(nr)

    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):
                rhox = charge[x, y, z] + charge_core[x, y, z]
                arhox = abs(rhox)
                if arhox > vanishing_charge:
                    ex, ec, vx, vc = pyqe.pyqe_xc(arhox, functional)
                    v[x, y, z] = 2.0 * (vx+vc)   # the factor 2.0 is e2 in a.u.

    return v
