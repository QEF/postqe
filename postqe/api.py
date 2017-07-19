#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tentative collection of functions to be part of postqe API and exposed to the user.

"""

from .pp import get_from_xml
from .readutils import read_charge_file_hdf5, write_charge, create_header
from .plot import plotcharge1D, plotcharge2D
from .compute_vs import compute_G
from .pyqe import pyqe_getcelldms


def plot_charge1D(xmlfile, plot_file='plotout', x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50):
    """

    :param xmlfile: File xml produced by the QE calculation
    :param plot_file: Output text file with the charge
    :param x0: 3D vector, origin of the line
    :param e1: 3D vector which determines the plotting line
    :param nx: number of points in the line
    :return:
    """
    prefix, outdir, ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, \
    nat, ntyp, lsda, noncolin, pseudodir, nr, nr_smooth = get_from_xml(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

    charge_file = outdir+prefix+".save/charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    if (lsda != 'true'):    # non magnetic calculation
        write_charge(plot_file, charge, header)

        # Plot a 1D section
        G = compute_G(b, charge.shape)
        fig = plotcharge1D(charge, G, a, x0, e1, nx)
        fig.show()
    else:                   # magnetic calculation, also plot charge_up and charge_down
        write_charge(plot_file, charge, header)
        charge_up = (charge + chargediff) / 2.0
        charge_down = (charge - chargediff) / 2.0
        write_charge(plot_file+'_up', charge_up, header)
        write_charge(plot_file+'_down', charge_down, header)

        # Plot 2D sections
        G = compute_G(b, charge.shape)
        fig1 = plotcharge1D(charge, G, a, x0, e1, nx)
        fig1.show()
        fig2 = plotcharge1D(charge_up, G, a, x0, e1, nx)
        fig2.show()
        fig3 = plotcharge1D(charge_down, G, a, x0, e1, nx)
        fig3.show()


def plot_charge2D(xmlfile, plot_file='plotout', x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50,
                  e2=(0., 1., 0.), ny=50):
    """

    :param xmlfile: File xml produced by the QE calculation
    :param plot_file: Output text file with the charge
    :param x0: 3D vector, origin of the plane
    :param e1, e2: 3D vectors which determine the plotting plane
    :param nx: number of points along the e1, e2 directions of the plotting plane
    :return:
    """
    prefix, outdir, ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, \
    nat, ntyp, lsda, noncolin, pseudodir, nr, nr_smooth = get_from_xml(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

    charge_file = outdir+"charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    if (lsda != 'true'):  # non magnetic calculation
        write_charge(plot_file, charge, header)
        # Plot a 2D section
        G = compute_G(b, charge.shape)
        fig = plotcharge2D(charge, G, a, x0, e1, e2, nx, ny)
        fig.show()
    else:  # magnetic calculation, also plot charge_up and charge_down
        write_charge(plot_file, charge, header)
        charge_up = (charge + chargediff) / 2.0
        charge_down = (charge - chargediff) / 2.0
        write_charge(plot_file + '_up', charge_up, header)
        write_charge(plot_file + '_down', charge_down, header)

        # Plot 2D sections
        G = compute_G(b, charge.shape)
        fig1 = plotcharge2D(charge, G, a, x0, e1, e2, nx, ny)
        fig1.show()
        fig2 = plotcharge2D(charge_up, G, a, x0, e1, e2, nx, ny)
        fig2.show()
        fig3 = plotcharge2D(charge_down, G, a, x0, e1, e2, nx, ny)
        fig3.show()

