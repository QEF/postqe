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







##
# CLI parser help messages
#
EOS_HELP = "Fit energy vs volume data with an equation of state."
BANDS_HELP = "Calculate energy bands."
DOS_HELP = "Calculate the electronic density of states"
CHARGE_HELP = """
Extract the charge from an output xml Espresso file and the corresponding HDF5 charge file 
containing the results of a calculation. Create also a Matplotlib figure object from a 1D or 2D
section of the charge. (optional) Export the charge (1, 2 or 3D section) in a text file 
according to different formats (XSF, cube, Gnuplot, etc.).
"""
POTENTIAL_HELP = """
Compute the potential specified in 'pot_type' from an output xml Espresso file and the 
corresponding HDF5 charge file containing the results of a calculation.
Create also a Matplotlib figure object from a 1D or 2D' section of the charge. 
(optional) Export the charge (1, 2 or 3D section) in a text file
according to different formats (XSF, cube, Gnuplot, etc.).
"""
PREFIX_HELP = "prefix of files saved by program pw.x"
OUTDIR_HELP = "directory containing the input data, i.e. the same as in pw.x"
SCHEMA_HELP = """the XSD schema file for QE XML output file. If not provided the schema
information is taken from xsi:schemaLocation attributes in the xml espresso file."""

EOS_PREFIX_HELP = "file containing the energy/volume data."
EOS_TYPE_HELP = """
type of equation of state (EOS) for fitting. Available types are:
murnaghan (default) -> Murnaghan EOS, PRB 28, 5480 (1983);
sjeos -> A third order inverse polynomial fit, PhysRevB.67.026103;
E(V) = c_0 + c_1 t + c_2 t^2  + c_3 t^3 ,  t = V^(-1/3);
taylor -> A third order Taylor series expansion around the minimum volume;
vinet -> Vinet EOS, PRB 70, 224107;
birch -> Birch EOS, Intermetallic compounds: Principles and Practice, Vol I: Principles, p.195;
birchmurnaghan -> Birch-Murnaghan EOS, PRB 70, 224107;
pouriertarantola -> Pourier-Tarantola EOS, PRB 70, 224107;
antonschmidt -> Anton-Schmidt EOS, Intermetallics 11, 23 - 32(2003);
p3 -> A third order inverse polynomial fit.
"""
EOS_FILEOUT_HELP = "text output file with fitting data and results (default="", not written)."
EOS_FILEPLOT_HELP = """
output plot file in png format (default='EOSplot'). Other formats are available from the Matplotlib GUI.
"""
EOS_SHOW_HELP = "True -> plot results with Matplotlib; None or False -> do nothing. Default = True."

BANDS_REFERENCE_ENERGY_HELP = "the Fermi level, defines the zero of the plot along y axis (default=0)."
BANDS_EMIN_HELP = "the minimum energy for the band plot (default=-50)."
BANDS_EMAX_HELP = "the maximum energy for the band plot (default=50)."
BANDS_FILEPLOT_HELP = "output plot file (default='bandsplot') in png format."

DOS_EMIN_HELP = 'the minimum energy for the dos plot (default=-50).'
DOS_EMAX_HELP = 'the maximum energy for the dos plot (default=50).'
DOS_NPTS_HELP = 'number of points of the DOS.'
DOS_FILEOUT_HELP = 'text output file with dos data (default='', not written).'
DOS_FILEPLOT_HELP = 'output plot file (default=\'dosplot\') in png format.'

VECTOR_HELP = " Enter the vector as components separated by commas, eg. 1,0,0 ."
CHARGE_FILEOUT_HELP = "text file with the full charge data as in the HDF5 file. Default='', nothing is written."
CHARGE_X0_HELP = "3D vector (a tuple), origin of the line or plane of the section." + VECTOR_HELP
CHARGE_E1_HELP = "1st 3D vector (a tuple) which determines the plotting section." + VECTOR_HELP
CHARGE_E2_HELP = "2nd 3D vector (a tuple) which determines the plotting section." + VECTOR_HELP
CHARGE_E3_HELP = "3rd 3D vector (a tuple) which determines the plotting section." + VECTOR_HELP
CHARGE_NX_HELP = "number of points along e1 direction."
CHARGE_NY_HELP = "number of points along e2 direction."
CHARGE_NZ_HELP = "number of points along e3 direction."
CHARGE_RADIUS_HELP = "radius of the sphere in the polar average method."
CHARGE_DIM_HELP = "1, 2, 3 for a 1D, 2D or 3D section respectively."
CHARGE_IFMAGN_HELP = """for a magnetic calculation, 'total' plot the total charge, 'up'
plot the charge with spin up, 'down' for spin down."""
CHARGE_EXPORTFILE_HELP = "file where plot data are exported in the chosen format (Gnuplot, XSF, cube Gaussian, etc.)."
CHARGE_METHOD_HELP = """
interpolation method. Available choices are:
'FFT' -> Fourier interpolation (default);
'polar' -> 2D polar plot on a sphere;
'spherical' -> 1D plot of the spherical average;
'splines' -> not implemented.
"""
CHARGE_FORMAT_HELP = """format of the (optional) exported file. Available choices are:
'gnuplot' -> plain text format for Gnuplot (default). Available for 1D and 2D sections.
'xsf' -> XSF format for the XCrySDen program. Available for 2D and 3D sections.
'cube' -> cube Gaussian format. Available for 3D sections.
'contour' -> format for the contour.x code of Quantum Espresso.
'plotrho' -> format for the plotrho.x code of Quantum Espresso.
"""
CHARGE_SHOW_HELP = "if True, show the Matplotlib plot (only for 1D and 2D sections)."

POT_TYPE_HELP = """type of the potential to calculate. Available types are:
'v_tot' (default) -> the total potential (v_bare+v_hartree+v_xc).
'v_bare' -> the bare potential.
'v_hartree' = the Hartree potential.
'v_xc' -> the exchange-correlation potential.
"""


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
    subparsers = parser.add_subparsers(help='sub-command help', dest='commands')

    # create the parser for the "eos" command
    eos_parser = subparsers.add_parser('eos', help=EOS_HELP)
    eos_parser.add_argument('-prefix', type=str, required=True, help=EOS_PREFIX_HELP)
    eos_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    eos_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    eos_parser.add_argument('-eos_type', type=str, default='murnaghan', help=EOS_TYPE_HELP)
    eos_parser.add_argument('-fileout', type=str, default='', help=EOS_FILEOUT_HELP)
    eos_parser.add_argument('-fileplot', type=str, default='EOSplot', help=EOS_FILEPLOT_HELP)
    eos_parser.add_argument('-show', type=bool, default=True, help=EOS_SHOW_HELP)

    # create the parser for the "bands" command
    bands_parser = subparsers.add_parser('bands', help=BANDS_HELP)
    bands_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    bands_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    bands_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    bands_parser.add_argument('-reference_energy', type=float, default=0, help=BANDS_REFERENCE_ENERGY_HELP)
    bands_parser.add_argument('-emin', type=float, default=-50, help=BANDS_EMIN_HELP)
    bands_parser.add_argument('-emax', type=float, default=50, help=BANDS_EMAX_HELP)
    bands_parser.add_argument('-fileplot', type=str, default='bandsplot', help=BANDS_FILEPLOT_HELP)
    bands_parser.add_argument('-show', type=bool, default=True, help=EOS_SHOW_HELP)

    # create the parser for the "dos" command
    dos_parser = subparsers.add_parser('dos', help=DOS_HELP)
    dos_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    dos_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    dos_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    dos_parser.add_argument('-width', type=float, default=0.01, help=SCHEMA_HELP)
    dos_parser.add_argument('-emin', type=float, default=None, help=DOS_EMIN_HELP)
    dos_parser.add_argument('-emax', type=float, default=None, help=DOS_EMAX_HELP)
    dos_parser.add_argument('-npts', type=int, default=100, help=DOS_NPTS_HELP)
    dos_parser.add_argument('-fileout', type=str, default='', help=DOS_FILEOUT_HELP)
    dos_parser.add_argument('-fileplot', type=str, default='dosplot', help=DOS_FILEPLOT_HELP)
    dos_parser.add_argument('-show', type=bool, default=True, help=EOS_SHOW_HELP)

    # create the parser for the "charge" command
    charge_parser = subparsers.add_parser('charge', help=CHARGE_HELP)
    charge_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    charge_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    charge_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    charge_parser.add_argument('-fileout', type=str, default='', help=CHARGE_FILEOUT_HELP)
    charge_parser.add_argument('-x0', type=vector, default=(0., 0., 0.), help=CHARGE_X0_HELP)
    charge_parser.add_argument('-e1', type=vector, default=(1., 0., 0.), help=CHARGE_E1_HELP)
    charge_parser.add_argument('-e2', type=vector, default=(0., 1., 0.), help=CHARGE_E2_HELP)
    charge_parser.add_argument('-e3', type=vector, default=(0., 0., 1.), help=CHARGE_E3_HELP)
    charge_parser.add_argument('-nx', type=int, default=20, help=CHARGE_NX_HELP)
    charge_parser.add_argument('-ny', type=int, default=20, help=CHARGE_NY_HELP)
    charge_parser.add_argument('-nz', type=int, default=20, help=CHARGE_NZ_HELP)
    charge_parser.add_argument('-radius', type=float, default=1, help=CHARGE_RADIUS_HELP)
    charge_parser.add_argument('-dim', type=int, default=1, help=CHARGE_DIM_HELP)
    charge_parser.add_argument('-ifmagn', type=str, default='total', help=CHARGE_IFMAGN_HELP)
    charge_parser.add_argument('-exportfile', type=str, default='', help=CHARGE_EXPORTFILE_HELP)
    charge_parser.add_argument('-method', type=str, default='FFT', help=CHARGE_METHOD_HELP)
    charge_parser.add_argument('-format', type=str, default='gnuplot', help=CHARGE_FORMAT_HELP)
    charge_parser.add_argument('-show', type=bool, default=True, help=CHARGE_SHOW_HELP)

    # create the parser for the "potential" command
    potential_parser = subparsers.add_parser('potential', help=POTENTIAL_HELP)
    potential_parser.add_argument('-prefix', type=str, required=True, help=PREFIX_HELP)
    potential_parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    potential_parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    potential_parser.add_argument('-pot_type', type=str, default='v_tot', help=POT_TYPE_HELP)
    potential_parser.add_argument('-fileout', type=str, default='', help=CHARGE_FILEOUT_HELP)
    potential_parser.add_argument('-x0', type=vector, default=(0., 0., 0.), help=CHARGE_X0_HELP)
    potential_parser.add_argument('-e1', type=vector, default=(1., 0., 0.), help=CHARGE_E1_HELP)
    potential_parser.add_argument('-e2', type=vector, default=(0., 1., 0.), help=CHARGE_E2_HELP)
    potential_parser.add_argument('-e3', type=vector, default=(0., 0., 1.), help=CHARGE_E3_HELP)
    potential_parser.add_argument('-nx', type=int, default=20, help=CHARGE_NX_HELP)
    potential_parser.add_argument('-ny', type=int, default=20, help=CHARGE_NY_HELP)
    potential_parser.add_argument('-nz', type=int, default=20, help=CHARGE_NZ_HELP)
    potential_parser.add_argument('-radius', type=float, default=1, help=CHARGE_RADIUS_HELP)
    potential_parser.add_argument('-dim', type=int, default=1, help=CHARGE_DIM_HELP)
    potential_parser.add_argument('-exportfile', type=str, default='', help=CHARGE_EXPORTFILE_HELP)
    potential_parser.add_argument('-method', type=str, default='FFT', help=CHARGE_METHOD_HELP)
    potential_parser.add_argument('-format', type=str, default='gnuplot', help=CHARGE_FORMAT_HELP)
    potential_parser.add_argument('-show', type=bool, default=True, help=CHARGE_SHOW_HELP)

    return parser


def main():
    from . import api

    if sys.version_info < (3, 4, 0):
        sys.stderr.write("You need python 3.4 or later to run this program\n")
        sys.exit(1)

    start_time = time.time()
    cli_parser = get_cli_parser()
    pars = cli_parser.parse_args()

    if pars.commands == 'eos':
        api.compute_eos(pars.prefix, pars.outdir, pars.eos_type, pars.fileout, pars.fileplot, pars.show)
    elif pars.commands == 'bands':
        api.compute_band_structure(
            pars.prefix, pars.outdir, pars.schema, pars.reference_energy, pars.emin,
            pars.emax, pars.fileplot, pars.show
        )
    elif pars.commands == 'dos':
        if pars.emin is None or pars.max is None:
            window = None
        else:
            window = (pars.emin, pars.emax)
        api.compute_dos(
            pars.prefix, pars.outdir, pars.schema, pars.width, window, pars.npts, pars.fileout,
            pars.fileplot, pars.show
        )
    elif pars.commands == 'charge':
        api.compute_charge(
            pars.prefix, pars.outdir, pars.schema, pars.fileout, pars.x0, pars.e1, pars.nx,
            pars.e2, pars.ny, pars.e3, pars.nz, pars.radius, pars.dim, pars.ifmagn, pars.exportfile,
            pars.method, pars.format, pars.show
        )
        api.new_get_charge(prefix, outdir, filplot)
    elif pars.commands == 'potential':
        api.compute_potential(
            pars.prefix, pars.outdir, pars.schema, pars.pot_type, pars.fileout, pars.x0,
            pars.e1, pars.nx, pars.e2, pars.ny, pars.e3, pars.nz, pars.radius, pars.dim,
            pars.exportfile, pars.method, pars.format, pars.show
        )
    else:
        print('Command not implemented! Exiting...')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished. Elapsed time: " + str(elapsed_time) + " s.")
