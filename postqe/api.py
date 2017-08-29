#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A tentative collection of functions to be part of postqe API and exposed to the user.

"""

import numpy as np
from ase.eos import EquationOfState
from ase.dft import DOS
from .charge import Charge, Potential
from .readutils import read_EtotV
from .ase.io import get_atoms_from_xml_output
from .ase.calculator import PostqeCalculator



def get_eos(label, eos='murnaghan'):
    """
    This function returns an EOS object from a text input file containing the volumes and corresponding calculated energies.
    Different equation of states are available: Murnaghan, Birch, Vinet, etc.

    :param label: input file for volumes and energies (possibly including the full path)
    :param eos_type: type of Equation of State (Murnaghan, Birch, Vinet, etc.)
    :return: an EOS object
    """

    # set a simple calculator, only to get calcul.prefix from the label
    calcul = PostqeCalculator(atoms=None, label=label)
    # Extract volumes and energies from the input file:
    volumes, energies = read_EtotV(calcul.prefix)

    # Create an object EquationOfState and fit with Murnaghan (or other) EOS
    eos = EquationOfState(volumes, energies, eos=eos)

    return eos

def get_band_structure(label, schema, reference_energy=0):

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    bs = atoms.calc.band_structure(reference=reference_energy)

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

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    # Create a DOS object with width= eV and npts points
    dos = DOS(calcul, width=width, npts=npts)

    return dos

def get_charge(label, schema):
    """
    This function returns an Charge object from an output xml Espresso file and the corresponding HDF5 charge file
    containing the results of a calculation.

    :param label: defines the system and the xml file containing the results (possibly including the full path)
    :param schema: the xml schema to be used to read and validate the xml output file
    :return a Charge object
    """

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    nr = calcul.get_nr()
    charge_file = calcul.prefix+".save/charge-density.hdf5"

    charge = Charge(nr)
    charge.read(charge_file)
    charge.set_calculator(calcul)

    return charge


def get_potential(label, schema, pot_type='vtot'):
    """
    This function returns an Potential object from an output xml Espresso file and the corresponding HDF5 charge file
    containing the results of a calculation.  The available potentials are the bare (pot_type='v_bare'), Hartree
      (pot_type='v_h'), exchange-correlation (pot_type='v_xc') and total (pot_type='v_tot').

    :param label: defines the system and the xml file containing the results (possibly including the full path)
    :param schema: the xml schema to be used to read and validate the xml output file
    :return a Potential object
    """

    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label=label, schema=schema)
    # define the Atoms structure reading the xml file
    atoms = get_atoms_from_xml_output(calcul.prefix + ".xml", schema=schema)
    atoms.set_calculator(calcul)
    # read the results
    atoms.calc.read_results()

    nr = calcul.get_nr()
    charge_file = calcul.prefix+".save/charge-density.hdf5"

    potential = Potential(nr)
    potential.read(charge_file)
    potential.set_calculator(calcul)
    potential.compute_potential(pot_type=pot_type)

    return potential
