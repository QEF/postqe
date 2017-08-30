#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
A tentative collection of functions to be part of postqe API and exposed to the user.
"""

import numpy as np
from .xmlfile import get_cell_data, get_calculation_data
from .readutils import read_charge_file_hdf5, write_charge, create_header
from .plot import plot1D_FFTinterp, plot2D_FFTinterp
from .compute_vs import compute_G, compute_v_bare, compute_v_h, compute_v_xc
from .pyqe import pyqe_getcelldms


def get_eos(label, eos='murnaghan'):
    """
    This function returns an EOS object from a text input file containing the volumes and corresponding calculated energies.
    Different equation of states are available: Murnaghan, Birch, Vinet, etc.

    :param label: input file for volumes and energies (possibly including the full path)
    :param eos_type: type of Equation of State (Murnaghan, Birch, Vinet, etc.)
    :return: an EOS object
    """
    from .readutils import read_EtotV
    from postqe.ase.calculator import PostqeCalculator
    from ase.eos import EquationOfState

    # set a simple calculator, only to get calcul.prefix from the label
    calcul = PostqeCalculator(atoms=None, label=label)
    # Extract volumes and energies from the input file:
    volumes, energies = read_EtotV(calcul.prefix)

    # Create an object EquationOfState and fit with Murnaghan (or other) EOS
    eos = EquationOfState(volumes, energies, eos=eos)

    return eos

def get_band_structure(label, schema, reference_energy=0):
    from postqe.ase.io import get_atoms_from_xml_output
    from postqe.ase.calculator import PostqeCalculator

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    Si = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    Si.set_calculator(calcul)
    # read the results
    Si.calc.read_results()

    bs = Si.calc.band_structure(reference=reference_energy)

    return bs


def get_dos(label, schema, width=0.01, npts=100):
    """
    This function returns an DOS object from an output xml Espresso file containing the results of a DOS calculation.

    :param label: defines the system and the xml file containing the results (possibly including the full path)
    :param schema: the xml schema to be used to read and validate the xml output file
    :param width: width of the gaussian to be used for the DOS (in eV)
    :param npts:  number of points of the DOS
    :return: a DOS object
    """
    from ase.dft import DOS
    from postqe.ase.io import get_atoms_from_xml_output
    from postqe.ase.calculator import PostqeCalculator

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    Si = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    Si.set_calculator(calcul)
    # read the results
    Si.calc.read_results()

    # Create a DOS object with width= eV and npts points
    dos = DOS(calcul, width=width, npts=npts)

    return dos


def get_charge(xmlfile, outfile='postqe.out'):
    """
    This function reads the *xmlfile* and the HDF5 charge file produced by Quantum Espresso and writes the charge values
     in a text file.

    :param xmlfile: file xml produced by the QE calculation
    :param outfile: output text file for the total charge. For magnetic systems, also writes the difference between
     charge density spin up - charge density spin down in outfile_diff.
    :return charge, chargediff : the electronic charge density and the difference between charge density spin up -
      charge density spin down
    """

    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)
    charge_file = outdir+prefix+".save/charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header(prefix, nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    write_charge(outfile, charge, header)

    if (lsda == 'true'):    # non magnetic calculation
        write_charge(outfile+'_diff', charge, header)

    return charge, chargediff

def get_potential(xmlfile, outfile='postqe.out', pot_type='vtot'):
    """
    This function reads the *xmlfile* and the HDF5 charge file produced by Quantum Espresso and writes the pot_type
     potential values in a text file. The available potentials are the bare (pot_type='vbare'), Hartree
      (pot_type='vh'), exchange-correlation (pot_type='vxc') and total (pot_type='vtot').

    :param xmlfile: file xml produced by the QE calculation
    :param outfile: output text file for the total charge. Default = 'postqe.out'
    :param pot_type: which type of potential ('vtot', 'vbare', 'vh', 'vxc'). Default = 'vtot'
    :return charge, chargediff : the electronic charge density and the difference between charge density spin up -
      charge density spin down
    """
    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

    charge_file = outdir + prefix + ".save/charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    if (pot_type=='vbare'):
        v = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir)
    elif (pot_type=='vh'):
        v = compute_v_h(charge, ecutrho, alat, b)
    elif (pot_type=='vxc'):
        charge_core = np.zeros(charge.shape)    # only for now, later in input
        v = compute_v_xc(charge, charge_core, str(functional))
    elif (pot_type=='vtot'):
        v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir)
        v_h =  compute_v_h(charge,ecutrho, alat, b)
        charge_core = np.zeros(charge.shape)    # only for now, later in input
        v_xc = compute_v_xc(charge, charge_core, str(functional))
        v = v_bare + v_h + v_xc

    write_charge(outfile, v, header)

    return v


def plot_charge1D(xmlfile, plot_file='plotout', x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50):
    """

    :param xmlfile: File xml produced by the QE calculation
    :param plot_file: Output text file with the charge
    :param x0: 3D vector, origin of the line
    :param e1: 3D vector which determines the plotting line
    :param nx: number of points in the line
    :return:
    """
    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

    charge_file = outdir+prefix+".save/charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header(prefix, nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    if (lsda != 'true'):    # non magnetic calculation
        write_charge(plot_file, charge, header)

        # Plot a 1D section
        G = compute_G(b, charge.shape)
        fig = plot1D_FFTinterp(charge, G, a, x0, e1, nx)
        fig.show()
    else:                   # magnetic calculation, also plot charge_up and charge_down
        write_charge(plot_file, charge, header)
        charge_up = (charge + chargediff) / 2.0
        charge_down = (charge - chargediff) / 2.0
        write_charge(plot_file+'_up', charge_up, header)
        write_charge(plot_file+'_down', charge_down, header)

        # Plot 2D sections
        G = compute_G(b, charge.shape)
        fig1 = plot1D_FFTinterp(charge, G, a, x0, e1, nx)
        fig1.show()
        fig2 = plot1D_FFTinterp(charge_up, G, a, x0, e1, nx)
        fig2.show()
        fig3 = plot1D_FFTinterp(charge_down, G, a, x0, e1, nx)
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
    ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species = get_cell_data(xmlfile)
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
        get_calculation_data(xmlfile)
    celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

    charge_file = outdir+"charge-density.hdf5"

    charge, chargediff = read_charge_file_hdf5(charge_file, nr)
    header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

    if (lsda != 'true'):  # non magnetic calculation
        write_charge(plot_file, charge, header)
        # Plot a 2D section
        G = compute_G(b, charge.shape)
        fig = plot2D_FFTinterp(charge, G, a, x0, e1, e2, nx, ny)
        fig.show()
    else:  # magnetic calculation, also plot charge_up and charge_down
        write_charge(plot_file, charge, header)
        charge_up = (charge + chargediff) / 2.0
        charge_down = (charge - chargediff) / 2.0
        write_charge(plot_file + '_up', charge_up, header)
        write_charge(plot_file + '_down', charge_down, header)

        # Plot 2D sections
        G = compute_G(b, charge.shape)
        fig1 = plot2D_FFTinterp(charge, G, a, x0, e1, e2, nx, ny)
        fig1.show()
        fig2 = plot2D_FFTinterp(charge_up, G, a, x0, e1, e2, nx, ny)
        fig2.show()
        fig3 = plot2D_FFTinterp(charge_down, G, a, x0, e1, e2, nx, ny)
        fig3.show()

