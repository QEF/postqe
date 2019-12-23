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
import argparse


SPIN_CHOICES = [0, 1, 2]
FERMI_LEVEL = 1.0  # FIXME: put in constants.py with an appropriate value


def vector(s):
    """Parses a 3-dimension vector argument."""
    try:
        x, y, z = map(int, s.split(','))
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("Vectors must be x,y,z")
    else:
        return x, y, z


def get_cli_parser():
    parser = argparse.ArgumentParser(description='QE post processing')
    subparsers = parser.add_subparsers(metavar="plot_num", dest='plot_num', required=True,
                                       help='Selects what to save in filplot')

    plot_num_parser = subparsers.add_parser(
        '0', aliases=['charge'], help="electron (pseudo-)charge density")
    plot_num_parser.add_argument(
        '-spin', type=int, default=0, choices=SPIN_CHOICES, help="spin component of charge")

    plot_num_parser = subparsers.add_parser(
        '1', help="total potential V_bare + V_H + V_xc")
    plot_num_parser.add_argument(
        '-spin', type=int, default=0, choices=SPIN_CHOICES, help="spin component of potential")

    subparsers.add_parser('2', help="local ionic potential V_bare")

    plot_num_parser = subparsers.add_parser(
        '3', help="local density of states at specific energy or grid of energies")
    plot_num_parser.add_argument(
        '-emin', type=float, default=FERMI_LEVEL, help="lower boundary of energy grid (in eV)")
    plot_num_parser.add_argument(
        '-emax', type=float, default=None, help="upper boundary of energy grid (in eV)")
    plot_num_parser.add_argument(
        '-delta_e', type=float, default=0.1, help="spacing of energy grid (in eV)")
    plot_num_parser.add_argument(
        '-degauss_ldos', type=float, default=None, help="broadening of energy levels for LDOS (in eV)")

    subparsers.add_parser('4', help="local density of electronic entropy")

    plot_num_parser = subparsers.add_parser(
        '5', help="STM images, Tersoff and Hamann, PRB 31, 805 (1985)")
    plot_num_parser.add_argument(
        '-sample_bias', type=float, help="the bias of the sample (Ry) in stm images")

    subparsers.add_parser('6', help="spin polarization (rho(up)-rho(down))")

    plot_num_parser = subparsers.add_parser(
        '7', help="contribution of selected wavefunction(s) to the (pseudo-)charge density")
    plot_num_parser.add_argument(
        '-kpoint', type=int, choices=[1, 2], action='append', help="k-point(s) to be plotted")
    plot_num_parser.add_argument(
        '-kband', type=int, choices=[1, 2], action='append', help="band(s) to be plotted")
    plot_num_parser.add_argument(
        '-lsign', type=bool, default=False, help="if true and k point is Gamma, plot |psi|^2 sign(psi)")
    plot_num_parser.add_argument(
        '-spin', type=int, choices=SPIN_CHOICES, action='append', help="spin component of potential")

    subparsers.add_parser('8', help="electron localization function (ELF)")
    subparsers.add_parser('9', help="charge density minus superposition of atomic densities")
    subparsers.add_parser('10', help="integrated local density of states (ILDOS)")
    subparsers.add_parser('11', help="the V_bare + V_H potential")
    subparsers.add_parser('12', help="the sawtooth electric field potential (if present)")
    subparsers.add_parser('13', help="the noncollinear magnetization")

    plot_num_parser = subparsers.add_parser(
        '17', help="all-electron valence charge density??")
    plot_num_parser.add_argument(
        '-spin', type=int, default=0, choices=SPIN_CHOICES, help="spin component of charge")

    subparsers.add_parser('18', help="the exchange and correlation magnetic field in the noncollinear case")
    subparsers.add_parser('19', help="Reduced density gradient")
    subparsers.add_parser(
        '20', help="Product of the electron density (charge) and the second "
                   "eigenvalue of the electron-density Hessian matrix")
    subparsers.add_parser('21', help="all-electron charge density (valence+core)")

    plot_num_parser = subparsers.add_parser(
        '22', help="kinetic energy density")
    plot_num_parser.add_argument(
        '-spin', type=int, default=0, choices=SPIN_CHOICES, help="spin component of density")

    parser.add_argument('-prefix', type=str, default='pwscf',
                        help="prefix of files saved by program pw.x")
    parser.add_argument('-outdir', type=str, default=None,
                        help="directory containing the input data, i.e. the same as in pw.x")
    parser.add_argument('-fileplot', type=str, default=None,
                        help="file to save the quantity selected by plot_num (use stdout instead)")

    return parser


def main():
    from . import api

    if sys.version_info < (3, 4, 0):
        sys.stderr.write("You need python 3.4 or later to run this program\n")
        sys.exit(1)

    start_time = time.time()
    cli_parser = get_cli_parser()
    cli_args = cli_parser.parse_args()

    kwargs = {k: getattr(cli_args, k) for k in dir(cli_args) if not k.startswith('_')}
    result = api.get_plot(**kwargs)
    if result is not None:
        print(result)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished. Elapsed time: " + str(elapsed_time) + " s.")
