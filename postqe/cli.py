#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
Command line interface for postqe.
"""
import sys
import time

from .api import get_charge, get_potential


def get_cli_parser_pp():
    import argparse

    parser = argparse.ArgumentParser(description='QE post processing')
    parser.add_argument('-plot_num', type=int, nargs='?', default=0, choices=range(0, 21),
                        help="""selects what to save in filplot:\n
    0  = electron (pseudo-)charge density\n
    1  = total potential V_bare + V_H + V_xc\n
    2  = local ionic potential V_bare\n
    3  = local density of states at e_fermi (number of states per volume, in bohr^3,\
         per energy unit, in Ry)\n\
    4  = local density of electronic entropy\n\
    5  = STM images\n\
        Tersoff and Hamann, PRB 31, 805 (1985)\n\
    6  = spin polarization (rho(up)-rho(down))\n

    7  = contribution of a selected wavefunction to the
        (pseudo-)charge density. For norm-conserving PPs,
        |psi|^2 (psi=selected wavefunction). Noncollinear case:
        contribution of the given state to the charge or
        to the magnetization along the direction indicated
        by spin_component (0 = charge, 1 = x, 2 = y, 3 = z )

    8  = electron localization function (ELF)

    9  = charge density minus superposition of atomic densities

    10 = integrated local density of states (ILDOS)
        from emin to emax (emin, emax in eV)
        if emax is not specified, emax=E_fermi

    11 = the V_bare + V_H potential

    12 = the sawtooth electric field potential (if present)

    13 = the noncollinear magnetization.

    17 = all-electron valence charge density
        can be performed for PAW calculations only
        requires a very dense real-space grid!

    18 = The exchange and correlation magnetic field in
        the noncollinear case

    19 = Reduced density gradient
        (J. Chem. Theory Comput. 7, 625 (2011))
        Set the isosurface between 0.3 and 0.6 to plot the
        non-covalent interactions (see also plot_num = 20)

    20 = Product of the electron density (charge) and the second
        eigenvalue of the electron-density Hessian matrix;
        used to colorize the RDG plot (plot_num = 19)

    21 = all-electron charge density (valence+core).
        For PAW calculations only; requires a very dense
        real-space grid.
    """)

    # default_prefix = "Si"
    default_prefix = "SiO2"
    # default_prefix = "SrTiO3"
    parser.add_argument('-prefix', type=str, nargs='?', default=default_prefix,
                        help='prefix of files saved by program pw.x')
    default_outdir = "../tests/" + default_prefix
    parser.add_argument('-outdir', type=str, nargs='?', default=default_outdir,
                        help='directory containing the input data, i.e. the same as in pw.x')
    parser.add_argument('-filplot', type=str, nargs='?', default="filplot",
                        help='file \"filplot\" contains the quantity selected by plot_num\
                        (can be saved for further processing)')

    parser.add_argument('-spin_component', type=int, nargs='?', default=0,
                        help='if plot_num==0: 0 = total charge (default value), 1 = spin up charge, 2 = spin down charge.\
                          if plot_num==1: 0 = spin averaged potential (default value), 1 = spin up potential,\
                          2 = spin down potential.')

    return parser


def get_cli_parser():
    import argparse
    PREFIX_HELP = 'prefix of files saved by program pw.x'
    OUTDIR_HELP = 'directory containing the input data, i.e. the same as in pw.x'
    SCHEMA_HELP = 'the XSD schema file for QE XML output file. If not provided the schema' \
                  'information is taken from xsi:schemaLocation attributes.'

    parser = argparse.ArgumentParser(description='QE post processing')
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "charge" command
    charge_parser = subparsers.add_parser(
        'charge', help='Get the charge from Espresso XML and HDF5 output.')
    charge_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    charge_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    charge_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)

    # create the parser for the "potential" command
    charge_parser = subparsers.add_parser(
        'potential', help='Get the potential from Espresso XML and HDF5 output.')
    charge_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    charge_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    charge_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    charge_parser.add_argument(
        '-pot_type', type=str, default=None, help="Type of the potential ('vtot', ...).")

    return parser


def main():
    if sys.version_info < (2, 7, 0):
        sys.stderr.write("You need python 2.7 or later to run this program\n")
        sys.exit(1)

    start_time = time.time()
    cli_parser = get_cli_parser()
    pars = cli_parser.parse_args()

    print(pars)

    from . import pp
    pp.run(pars)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: " + str(elapsed_time) + " s.")
