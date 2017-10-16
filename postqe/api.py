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
A collection of functions to be part of postqe API and exposed to the user.
"""
import os
from postqe.eos import QEEquationOfState
from postqe.dos import QEDOS

from .charge import Charge, Potential
from .readutils import read_EtotV
from .ase.io import get_atoms_from_xml_output
from .ase.calculator import PostqeCalculator


def get_label(prefix, outdir=None):
    if outdir is None:
        try:
            outdir = os.environ['ESPRESSO_TMPDIR']
        except KeyError:
            outdir = os.curdir

    label = os.path.join(outdir, prefix)
    return label


def get_eos(prefix, outdir=None, eos_type='murnaghan'):
    """
    Returns an EOS object from a text input file containing the volumes
    and corresponding calculated energies.
    Different equation of states are available: Murnaghan, Birch, Vinet, etc.

    :param prefix: prefix of saved files for volumes and energies
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set, or current directory ('.') otherwise
    :param eos_type: type of equation of state (EOS) for fitting. Available types are:\n
            'murnaghan' (default) -> Murnaghan EOS, PRB 28, 5480 (1983)\n
            'sjeos' -> A third order inverse polynomial fit, PhysRevB.67.026103\n
            \t\tE(V) = c_0 + c_1 t + c_2 t^2  + c_3 t^3 ,  t = V^(-1/3)\n
            'taylor' -> A third order Taylor series expansion around the minimum volume\n
            'vinet' -> Vinet EOS, PRB 70, 224107 \n
            'birch' -> Birch EOS, Intermetallic compounds: Principles and Practice, Vol I: Principles, p. 195\n
            'birchmurnaghan' -> Birch-Murnaghan EOS, PRB 70, 224107\n
            'pouriertarantola' -> Pourier-Tarantola EOS, PRB 70, 224107\n
            'antonschmidt' -> Anton-Schmidt EOS, Intermetallics 11, 23 - 32(2003)\n
            'p3' -> A third order inverse polynomial fit\n

    :return: an EOS object
    """
    label = get_label(prefix, outdir)
    # Extract volumes and energies from the input file:
    volumes, energies = read_EtotV(label)
    # Create an object EquationOfState and fit with Murnaghan (or other) EOS
    eos = QEEquationOfState(volumes, energies, eos=eos_type)
    return eos


def compute_eos(prefix, outdir=None, eos_type='murnaghan', fileout='', fileplot='EOSplot', show=True, ax=None):
    """
    This function fits an Equation of state of type *eos* and writes the results into *filename*.
    Different equation of states are available: Murnaghan, Birch, Vinet, etc.

    :param prefix: prefix of saved files for volumes and energies
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set, or current directory ('.') otherwise
    :param eos_type: type of equation of state (EOS) for fitting. Available types are:\n
            'murnaghan' (default) -> Murnaghan EOS, PRB 28, 5480 (1983)\n
            'sjeos' -> A third order inverse polynomial fit, PhysRevB.67.026103\n
            \t\tE(V) = c_0 + c_1 t + c_2 t^2  + c_3 t^3 ,  t = V^(-1/3)\n
            'taylor' -> A third order Taylor series expansion around the minimum volume\n
            'vinet' -> Vinet EOS, PRB 70, 224107 \n
            'birch' -> Birch EOS, Intermetallic compounds: Principles and Practice, Vol I: Principles, p. 195\n
            'birchmurnaghan' -> Birch-Murnaghan EOS, PRB 70, 224107\n
            'pouriertarantola' -> Pourier-Tarantola EOS, PRB 70, 224107\n
            'antonschmidt' -> Anton-Schmidt EOS, Intermetallics 11, 23 - 32(2003)\n
            'p3' -> A third order inverse polynomial fit
    :param fileout: output file with fitting data and results (default='', not written).
    :param fileplot: output plot file (default='EOSplot') in png format.
    :param show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    :param ax: a Matplotlib "Axes" instance (see Matplotlib documentation for details). If ax=None (default), creates
            a new one
    :return: an EOS object and a Matplotlib figure object
    """

    eos = get_eos(prefix, outdir, eos_type)
    v0, e0, B = eos.fit()
    if fileout !='':
        eos.write(fileout)
    fig = eos.plot(fileplot, show=show, ax=ax)

    return eos, fig

def get_band_structure(prefix, outdir=None, schema=None, reference_energy=0):
    """
    This function returns a "band structure" object from an output xml Espresso file
    containing the results of a proper calculation along a path in the Brilluoin zone.

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param reference_energy: the Fermi level
    :return: a band structure object
    """
    label = get_label(prefix, outdir)

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.label + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    bs = atoms.calc.band_structure(reference=reference_energy)

    return bs


def get_dos(prefix, outdir=None, schema=None, width=0.01, window= None, npts=100):
    """
    This function returns an DOS object from an output xml Espresso file containing the
    results of a DOS calculation.

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param width: width of the gaussian to be used for the DOS (in eV)
    :param window = emin, emax: defines the minimun and maximun energies for the DOS
    :param npts:  number of points of the DOS
    :return: a DOS object
    """
    label = get_label(prefix, outdir)

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.label + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    # Create a DOS object with width= eV and npts points
    dos = QEDOS(calcul, width=width, window=window, npts=npts)

    return dos


def get_charge(prefix, outdir=None, schema=None):
    """
    This function returns an Charge object from an output xml Espresso file and the
    corresponding HDF5 charge file containing the results of a calculation.

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :return: a Charge object
    """
    label = get_label(prefix, outdir)

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.label + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    nr = calcul.get_nr()
    charge_file = calcul.label + ".save/charge-density.hdf5"

    charge = Charge(nr)
    charge.read(charge_file)
    charge.set_calculator(calcul)

    return charge


def get_potential(prefix, outdir=None, schema=None, pot_type='vtot'):
    """
    This function returns an Potential object from an output xml Espresso file and
    the corresponding HDF5 charge file containing the results of a calculation.
    The available potentials are the bare (pot_type='v_bare'), Hartree (pot_type='v_h'),
    exchange-correlation (pot_type='v_xc') and total (pot_type='v_tot').

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param pot_type: type of the Potential ('vtot', ....)
    :return: a Potential object
    """
    label = get_label(prefix, outdir)

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.label + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    nr = calcul.get_nr()
    charge_file = calcul.label + ".save/charge-density.hdf5"

    potential = Potential(nr)
    potential.read(charge_file)
    potential.set_calculator(calcul)
    potential.compute_potential(pot_type=pot_type)

    return potential

