#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A set of utility functions to extract data from the xml file produced by QE.
"""

import numpy as np
import xmlschema


def get_dict(xmlfile):

    ##########################################################
    # TODO for whatever reason this is not working now
    # schemaLoc = xmlschema.fetch_schema(xmlfile)
    # xs = xmlschema.XMLSchema(schemaLoc)
    #
    # temporary local solution
    xs = xmlschema.XMLSchema('schemas/qes.xsd')
    ##########################################################

    d = xs.to_dict(xmlfile)

    return d


def get_cell_data(xmlfile):
    """
    Gets some data about the unit cell from the xmlfile.

    :param xmlfile:
    :return ibrav, alat, a, b:
    """
    d = get_dict(xmlfile)
    dout = d["output"]
    alat = (dout["atomic_structure"]["@alat"])
    a1 = np.array(dout["atomic_structure"]["cell"]["a1"])
    a2 = np.array(dout["atomic_structure"]["cell"]["a2"])
    a3 = np.array(dout["atomic_structure"]["cell"]["a3"])
    ibrav = (dout["atomic_structure"]["@bravais_index"])
    b1 = np.array(dout["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(dout["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(dout["basis_set"]["reciprocal_lattice"]["b3"])
    a = np.array([a1, a2, a3])
    b = np.array([b1, b2, b3])
    a_p = (dout["atomic_structure"]["atomic_positions"]["atom"])
    a_s = (dout["atomic_species"]["species"])
    nat = (dout["atomic_structure"]["@nat"])
    ntyp = (dout["atomic_species"]["@ntyp"])

    # for subsequent loops it is important to have always lists for atomic_positions
    # and atomic_species. If this is not, convert
    if (type(a_s) == type([])):
        atomic_species = a_s
    else:
        atomic_species = [a_s]

    if (type(a_p) == type([])):
        atomic_positions = a_p
    else:
        atomic_positions = [a_p]

    return ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species


def get_band_strucure_data(xmlfile):
    """
    Get some useful values from xml file
    """

    d = get_dict(xmlfile)
    dout = d["output"]
    nks = (dout["band_structure"]["nks"])
    nbnd = (dout["band_structure"]["nbnd"])
    ks_energies = (dout["band_structure"]["ks_energies"])

    return nks, nbnd, ks_energies


def get_calculation_data(xmlfile):
    """
    Get some calculation data from xml file
    """

    d = get_dict(xmlfile)
    dout = d["output"]
    prefix = (d["input"]["control_variables"]["prefix"])
    outdir = (d["input"]["control_variables"]["outdir"])
    ecutwfc = (dout["basis_set"]["ecutwfc"])
    ecutrho = (dout["basis_set"]["ecutrho"])

    functional = np.array(dout["dft"]["functional"])

    lsda = (dout["magnetization"]["lsda"])
    noncolin = (dout["magnetization"]["noncolin"])
    nr = np.array([dout["basis_set"]["fft_grid"]["@nr1"], dout["basis_set"]["fft_grid"]["@nr2"],
                   dout["basis_set"]["fft_grid"]["@nr3"]])
    nr_smooth = np.array([dout["basis_set"]["fft_smooth"]["@nr1"], dout["basis_set"]["fft_smooth"]["@nr2"],
                          dout["basis_set"]["fft_smooth"]["@nr3"]])
    pseudodir = d["input"]["control_variables"]["pseudo_dir"]


    return (prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth)
