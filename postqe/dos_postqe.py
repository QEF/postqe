#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to calculate the electronic density of states (DOS).

TODO: Tetrahedra methods (linear and optimized).
"""

import numpy as np
from .constants import ev_to_ry
from .xmlfile import get_cell_data, get_calculation_data, get_band_strucure_data
from .pyqe import py_w0gauss


def dos_gaussian(e, nat, ks_energies, lsda, nbnd, nks, degauss, ngauss=0):
    """
    Calculated the electronic density of states with Gaussian broadening.

    :param e energy values (for which calculate the dos)
    :param ks_energies: eigenvalues with weights and k-points
    :param lsda: if true = magnetic calculation
    :param nbnd: number of bands
    :param nks: number of k-points
    :param degauss: value for the Gaussian broadening
    :param ngauss:  0   -> Simple Gaussian (default)
                    1   -> Methfessel-Paxton of order 1
                    -1  -> Marzari-Vanderbilt "cold smearing"
                    -99 -> Fermi-Dirac function
    :return: dos_up, dos_down (for the given value of energy)
    """

    # TODO non collinear case to be implemented
    dos_up = 0.
    dos_down = 0.

    if lsda:   # if magnetic
        for i in range(0, nks):
            weight = ks_energies[i]['k_point']['@weight']  # weight at k-point i
            for j in range(0, nbnd // 2):
                eigenvalue = ks_energies[i]['eigenvalues'][j] * 2 * nat           # eigenvalue at k-point i, band j
                dos_up += weight * py_w0gauss((e - eigenvalue) / degauss, ngauss)
            for j in range(nbnd // 2, nbnd):
                eigenvalue = ks_energies[i]['eigenvalues'][j] * 2 * nat        # eigenvalue at k-point i, band j
                dos_down += weight * py_w0gauss((e - eigenvalue) / degauss, ngauss)

    else:       # non magnetic
        for i in range(0, nks):
            weight = ks_energies[i]['k_point']['@weight']  # weight at k-point i
            for j in range(0, nbnd):
                eigenvalue = ks_energies[i]['eigenvalues'][j] * nat           # eigenvalue at k-point i, band j
                dos_up += weight * py_w0gauss((e - eigenvalue) / degauss, ngauss)

    dos_up /= degauss
    dos_down /= degauss

    return dos_up, dos_down


def compute_dos(xmlfile, filedos='filedos', e_min='', e_max='', e_step=0.01, degauss=0.02, ngauss=0):
    """
    Compute the electronic density of states between *e_min* and *e_max*, with step *e_step*, using
    Gaussian broadening type *ngauss* and value *degauss*.

    :param xmlfile: xml output file from QE
    :param filedos: output dos file containing the DOS values
    :param e_min, e_max, e_step: calculate the DOS in the range (*e_min*,*e_max*) with step *e_step*
    :param degauss: value for the Gaussian broadening
    :param ngauss:  0   -> Simple Gaussian (default)
                    1   -> Methfessel-Paxton of order 1
                    -1  -> Marzari-Vanderbilt "cold smearing"
                    -99 -> Fermi-Dirac function
    :return: np.array(e), np.array(dos_up), np.array(dos_down)
    """

    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    nks, nbnd, ks_energies = get_band_strucure_data(xmlfile)

    # TODO determine E_min, E_max automatically from ks_energies if not set in input parameters

    # Convert to rydberg
    e_min = e_min * ev_to_ry
    e_max = e_max * ev_to_ry
    e_step = e_step * ev_to_ry

    fout = open(filedos, "w")
    fout.write("#  E"+18*' '+" dos up "+13*' '+" dos down \n")
    fout.write("#  (eV)"+16*' '+"(states/eV/cell)"+4*' '+" (states/eV/cell)\n")
    te = e_min
    e = []
    dos_up = []
    dos_down = []
    while te < e_max:
        tdos_up, tdos_down = dos_gaussian(te, nat, ks_energies, lsda, nbnd, nks, degauss, ngauss)
        fout.write( "  {:.9E}".format(te / ev_to_ry)+4*' '+" {:.9E}".format(tdos_up * ev_to_ry) +
                    4*' '+"  {:.9E}\n".format(tdos_down * ev_to_ry) )
        e.append(te / ev_to_ry)
        dos_up.append(tdos_up * ev_to_ry)
        dos_down.append(tdos_down * ev_to_ry)
        te += e_step

    return np.array(e), np.array(dos_up), np.array(dos_down)