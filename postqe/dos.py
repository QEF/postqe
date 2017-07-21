#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to calculate the electronic density of states (DOS).
"""

from postqe.xmlfile import get_cell_data, get_calculation_data, get_band_strucure_data
from postqe.pyqe import py_w0gauss


def dos_gaussian(E, nat, ks_energies, lsda, nbnd, nks, degauss, ngauss=0):
    """
    Calculated the electronic density of states with Gaussian broadening.

    :param E energy values (for which calculate the dos)
    :param ks_energies: eigenvalues with weights and k-points
    :param lsda: if true = magnetic calculation
    :param nbnd: number of bands
    :param nks: number of k-points
    :param degauss: value for the Gaussian smearing
    :param ngauss:  0   -> Simple Gaussian (default)
                    1   -> Methfessel-Paxton of order 1
                    -1  -> Marzari-Vanderbilt "cold smearing"
                    -99 -> Fermi-Dirac function
    :return: dos_up, dos_down
    """

    # TODO non collinear case to be implemented
    nk = nks
    if lsda == 'true':
        nk = nks // 2

    dos_up = 0.
    dos_down = 0.

    for i in range(0,nk):
        for j in range(0,nbnd):
            weight = ks_energies[i]['k_point']['@weight']                  # weight at k-point i
            eigenvalue = ks_energies[i]['eigenvalues'][j]  * nat           # eigenvalue at k-point i, band j
            dos_up += weight * py_w0gauss( (E-eigenvalue)/degauss, ngauss )

    if lsda == 'true':
        for i in range(nk,nks):
            for j in range(0,nbnd):
                weight = ks_energies[i]['k_point']['@weight']              # weight at k-point i
                eigenvalue = ks_energies[i]['eigenvalues'][j] * nat        # eigenvalue at k-point i, band j
                dos_down += weight * py_w0gauss( (E-eigenvalue)/degauss, ngauss )

    dos_up /= degauss
    dos_down /= degauss

    return dos_up, dos_down


def compute_dos(xmlfile,filedos='filedos',E_min='',E_max='',E_step=0.01, degauss=0.02, ngauss=0):

    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    nks, nbnd, ks_energies = get_band_strucure_data(xmlfile)

    # TODO determine E_min, E_max automatically from ks_energies if not set in input parameters

    # Convert to rydberg
    ev_to_ry = 0.073498618
    E_min = E_min * ev_to_ry
    E_max = E_max * ev_to_ry
    E_step = E_step * ev_to_ry


    fout = open(filedos, "w")
    E = E_min
    while E < E_max:
        dos_up, dos_down = dos_gaussian(E, nat, ks_energies, lsda, nbnd, nks, degauss, ngauss)
        fout.write( "{:.9E}".format(E / ev_to_ry)+"  {:.9E}\n".format(dos_up * ev_to_ry) )
        E += E_step




compute_dos('silicon.xml','testdos.dat',-10,20,0.01)