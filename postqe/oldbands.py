#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to calculate the electronic band structure.

Note: no symmetry recognition is implemented yet.
"""

import numpy as np
from math import fabs, sqrt
from postqe.xmlfile import get_cell_data, get_calculation_data, get_band_strucure_data
from postqe.constants import ev_to_ry


def compute_bands(xmlfile, filebands='filebands', spin_component=''):
    """
    """

    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    nks, nbnd, ks_energies = get_band_strucure_data(xmlfile)

    # open output file
    fout = open(filebands, "w")
    fout.write("& plot  nbnd = "+str(nbnd)+" nks = "+str(nks)+" /\n")

    kpoints = np.zeros((nks, 3))
    bands = np.zeros((nks, nbnd))
    if lsda:   # magnetic
        for i in range(0, nks):
            kpoints[i] = ks_energies[i]['k_point']['$']
            fout.write(12 * ' ' + ' {:.6E}'.format(kpoints[i,0]) + ' {:.6E}'.format(kpoints[i,1]) + ' {:.6E}\n'.format(kpoints[i,2]))
            if (spin_component==1):     # get bands for spin up
                for j in range(0, nbnd // 2):
                    bands[i,j] = ks_energies[i]['eigenvalues'][j] * 2 * nat / ev_to_ry          # eigenvalue at k-point i, band j
                    fout.write('   {:.3E}'.format(bands[i,j]))
            else:                       # get bands for spin down
                for j in range(nbnd // 2, nbnd):
                    bands[i, j] = ks_energies[i]['eigenvalues'][j] * 2 * nat / ev_to_ry          # eigenvalue at k-point i, band j
                    fout.write('   {:.3E}'.format(bands[i,j]))
            fout.write('\n')

    else:       # non magnetic
        for i in range(0, nks):
            kpoints[i] = ks_energies[i]['k_point']['$']
            fout.write(12 * ' ' + ' {:.6E}'.format(kpoints[i,0]) + ' {:.6E}'.format(kpoints[i,1]) + ' {:.6E}\n'.format(kpoints[i,2]))
            for j in range(0, nbnd):
                bands[i, j] = ks_energies[i]['eigenvalues'][j] * nat / ev_to_ry           # eigenvalue at k-point i, band j
                fout.write('   {:.3E}'.format(bands[i,j]))
            fout.write('\n')

    return kpoints, bands


def set_high_symmetry_points(kpoints):
    """
    Determines which k-points have "high simmetry" and are at the boundaries of the Brillouin zone.

    :param kpoints: a matrix (nks,3) with the k-points coordinates. nks is the number of k-points.
    :return high_sym: an array of nks booleans, True if the kpoint is a high symmetry one
    """
    nks = kpoints.shape[0]
    high_sym = np.full(nks,False,dtype=bool)
    high_sym[0] = True
    high_sym[nks-1] = True

    k1 = np.zeros(3)
    k2 = np.zeros(3)
    for i in range(1,nks-1):
        if np.dot(kpoints[i,:],kpoints[i,:]) < 1.e-9:   # the Gamma point is always a high symmetry one
            high_sym[i] = True
        else:
            k1 = kpoints[i,:] - kpoints[i-1,:]
            k2 = kpoints[i+1,:] - kpoints[i,:]
            ps = np.dot(k1,k2) / sqrt(np.dot(k1,k1)) / sqrt(np.dot(k2,k2))
            if fabs(ps-1.0) > 1.0e-4 :
                high_sym[i] = True

    return high_sym


def compute_kx(kpoints):
    """
    This functions "linearize" the path along the k-points list in input and calculate
    the linear x variable kx for the plot.

    :param kpoints: a matrix (nks,3) with the k-points coordinates. nks is the number of k-points.
    :return kx : linear x variable for the plot determined as the k-points path
    """

    nks = kpoints.shape[0]
    kx = np.zeros(nks)

    ktemp = kpoints[2, :] - kpoints[1, :]
    dxmod_save = sqrt (np.dot(ktemp, ktemp))

    for i in range(1,nks):
        ktemp = kpoints[i, :] - kpoints[i-1, :]
        dxmod = sqrt(np.dot(ktemp, ktemp))

        if dxmod > 5*dxmod_save:    # a big jump in dxmod is a sign the points kpoints[i] and kpoints[i]
                                    # are quite distant and belong to two different lines. We put them on
                                    # the same point in the graph
            kx[i] = kx[i-1]
        elif dxmod > 1.e-5:         # this is the usual case. The two points kpoints[i] and kpoints[i] are in the
                                    # same path.
            kx[i] = kx[i-1] + dxmod
            dxmod_save = dxmod
        else:                       # !  This is the case in which dxmod is almost zero. The two points coincide
                                    # in the graph, but we do not save dxmod.
            kx[i] = kx[i-1] + dxmod

    return kx

