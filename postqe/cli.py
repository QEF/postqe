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
DOS_HELP = "Calculate the electronic density of states."
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
PREFIX_HELP = "Prefix of files saved by program pw.x."
OUTDIR_HELP = "Directory containing the input data, i.e. the same as in pw.x"
SCHEMA_HELP = """The XSD schema file for QE XML output file. If not provided the schema
information is taken from xsi:schemaLocation attributes in the xml espresso file."""
SHOW_HELP = """
True -> plot results with Matplotlib; None or False -> do nothing (default = True).
For DOS only for 1D and 2D sections.
"""
FILEPLOT_HELP = "Output plot file in png format (default = 'plot'). Other formats are available from the Matplotlib GUI."

EOS_PREFIX_HELP = "Name of the file containing the energy/volume data."
EOS_TYPE_HELP = """
Type of equation of state (EOS) for fitting. Available types are:
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
EOS_FILEOUT_HELP = "Text output file with fitting data and results (default = 'eos_data.dat')."
EOS_AX_HELP = "A Matplotlib 'Axes' instance (see Matplotlib documentation for details (default = None, creates a new one)."

BANDS_REFERENCE_ENERGY_HELP = "The Fermi level, defines the zero of the plot along y axis (default = 0)."
BANDS_EMIN_HELP = "The minimum energy for the band plot (default = -50)."
BANDS_EMAX_HELP = "The maximum energy for the band plot (default = 50)."

DOS_WINDOW_HELP = "A tuple (emin, emax) that defines the minimum and maximum energies for the DOS (default = None)."
DOS_WIDTH_HELP = "Width of the gaussian to be used for the DOS (in eV, default = 0.5)."
DOS_NPTS_HELP = 'Number of points of the DOS (default = 100).'
DOS_FILEOUT_HELP = "Text output file with dos data (default = 'dos')."

VECTOR_HELP = "Enter the vector as components separated by commas, eg. 1,0,0."
CHARGE_FILEOUT_HELP = "Text file with the full charge data as in the HDF5 file (default='charge.dat')."
CHARGE_X0_HELP = "3D vector (a tuple), origin of the line or plane of the section (default = 0,0,0)." + VECTOR_HELP
CHARGE_E1_HELP = "1st 3D vector (a tuple) which determines the plotting section (default = 1,0,0)." + VECTOR_HELP
CHARGE_E2_HELP = "2nd 3D vector (a tuple) which determines the plotting section (default = 0,1,0)." + VECTOR_HELP
CHARGE_E3_HELP = "3rd 3D vector (a tuple) which determines the plotting section (default = 0,0,1)." + VECTOR_HELP
CHARGE_NX_HELP = "Number of points along e1 direction (default = 50)."
CHARGE_NY_HELP = "Number of points along e2 direction (default = 50)."
CHARGE_NZ_HELP = "Number of points along e3 direction (default = 50)."
CHARGE_RADIUS_HELP = "Radius of the sphere in the polar average method (default = 1)"
CHARGE_DIM_HELP = "1, 2, 3 for a 1D, 2D or 3D section respectively (default = 1)"
CHARGE_IFMAGN_HELP = """
For a magnetic calculation, 'total' plot the total charge, 'up'
plot the charge with spin up, 'down' for spin down (default=None).
"""
CHARGE_EXPORTFILE_HELP = """
Name of the file where plot data are exported in the chosen format (Gnuplot, XSF, cube Gaussian, etc.)
(default = 'plot_data.dat').
"""
CHARGE_METHOD_HELP = """
Interpolation method. Available choices are:
'FFT' -> Fourier interpolation (default);
'polar' -> 2D polar plot on a sphere;
'spherical' -> 1D plot of the spherical average;
'splines' -> not implemented.
"""
CHARGE_FORMAT_HELP = """
Format of the (optional) exported file. Available choices are:
'gnuplot' -> plain text format for Gnuplot (default). Available for 1D and 2D sections.
'xsf' -> XSF format for the XCrySDen program. Available for 2D and 3D sections.
'cube' -> cube Gaussian format. Available for 3D sections.
'contour' -> format for the contour.x code of Quantum Espresso.
'plotrho' -> format for the plotrho.x code of Quantum Espresso.
"""

# SPIN_CHOICES = [0, 1, 2]
SPIN_CHOICES = ['total', 'up', 'down']
METHOD_CHOICES = ['FFT', 'polar', 'spherical']
EOS_CHOICES = ['murnaghan', 'sjeos', 'taylor', 'vinet', 'birch', 'birchmurnaghan', 'pouriertarantola', 'antonschmidt', 'p3']
FORMAT_CHOICES = ['gnuplot', 'xsf', 'cube', 'contour', 'plotrho']
BOOL_CHOICES = [True, False]
DIM_CHOICES = [1, 2, 3]
FERMI_LEVEL = 1.0  # FIXME: put in constants.py with an appropriate value


def vector(s):
    """Parses a 3-dimension vector argument."""
    try:
        x, y, z = map(float, s.split(','))
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("Vectors must be x,y,z")
    else:
        return x, y, z
#
def window(s):
    """Parses a 2-dimension argument."""
    try:
        emin, emax = map(float, s.split(','))
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("Tuple must be: emin, emax")
    else:
        return emin, emax

def get_cli_parser():
    parser = argparse.ArgumentParser(description='QE post processing')
    subparsers = parser.add_subparsers(metavar="command", dest='command', required=True,
                                       help='Selects what to save in filplot')

    #COMPUTE EOS
    command_parser = subparsers.add_parser(
        'eos', help=EOS_HELP)
    #SPECIFIC OPTIONS
    command_parser.add_argument(
        '-eos_type', type=str, default='murnaghan', choices=EOS_CHOICES, help=EOS_TYPE_HELP)
    command_parser.add_argument(
        '-fileout', type=str, default='eos_data.dat', help=EOS_FILEOUT_HELP)
    command_parser.add_argument(
        '-ax', type=str, default=None, help=EOS_AX_HELP)

    #COMPUTE BANDS
    command_parser = subparsers.add_parser(
        'bands', help=BANDS_HELP)
    #SPECIFIC OPTIONS
    command_parser.add_argument(
        '-reference_energy', type=float, default=0, help=BANDS_REFERENCE_ENERGY_HELP)
    command_parser.add_argument(
        '-emin', type=float, default=-50, help=BANDS_EMIN_HELP)
    command_parser.add_argument(
        '-emax', type=float, default=50, help=BANDS_EMAX_HELP)

    #COMPUTE DOS
    command_parser = subparsers.add_parser(
        'dos', help=DOS_HELP)
    #SPECIFIC OPTIONS
    command_parser.add_argument(
        '-window', type=window, default=None, help=DOS_WINDOW_HELP)
    command_parser.add_argument(
        '-width', type=float, default=0.5, help=DOS_WIDTH_HELP)
    command_parser.add_argument(
        '-npts', type=int, default=100, help=DOS_NPTS_HELP)
    command_parser.add_argument(
        '-fileout', type=str, default='dos', help=DOS_FILEOUT_HELP)

    #COMPUTE CHARGE
    command_parser = subparsers.add_parser(
        'charge', help=CHARGE_HELP)
    #SPECIFIC OPTIONS
    command_parser.add_argument(
        '-fileout', type=str, default='charge.dat', help=CHARGE_FILEOUT_HELP)
    command_parser.add_argument(
        '-ifmagn', type=str, default=None, choices=SPIN_CHOICES, help=CHARGE_IFMAGN_HELP)
    command_parser.add_argument(
        '-x0', type=vector, default=(0.,0.,0.), help=CHARGE_X0_HELP)
    command_parser.add_argument(
        '-e1', type=vector, default=(1.,0.,0.), help=CHARGE_E1_HELP)
    command_parser.add_argument(
        '-nx', type=int, default=50, help=CHARGE_NX_HELP)
    command_parser.add_argument(
        '-e2', type=vector, default=(0.,1.,0.), help=CHARGE_E2_HELP)
    command_parser.add_argument(
        '-ny', type=int, default=50, help=CHARGE_NY_HELP)
    command_parser.add_argument(
        '-e3', type=vector, default=(0.,0.,1.), help=CHARGE_E3_HELP)
    command_parser.add_argument(
        '-nz', type=int, default=50, help=CHARGE_NZ_HELP)
    command_parser.add_argument(
        '-radius', type=int, default=1, help=CHARGE_RADIUS_HELP)
    command_parser.add_argument(
        '-dim', type=int, default=1, choices=DIM_CHOICES, help=CHARGE_DIM_HELP)
    command_parser.add_argument(
        '-method', type=str, default='FFT', choices=METHOD_CHOICES, help=CHARGE_METHOD_HELP)
    command_parser.add_argument(
        '-format', type=str, default='gnuplot', choices=FORMAT_CHOICES, help=CHARGE_FORMAT_HELP)
    command_parser.add_argument(
        '-plot_file', type=str, default='plot_data.dat', help=CHARGE_EXPORTFILE_HELP)


    parser.add_argument('-prefix', type=str, default='pwscf', help=PREFIX_HELP)
    parser.add_argument('-outdir', type=str, default=None, help=OUTDIR_HELP)
    parser.add_argument('-fileplot', type=str, default='plot', help=FILEPLOT_HELP)
    parser.add_argument('-schema', type=str, default=None, help=SCHEMA_HELP)
    parser.add_argument('-show', type=bool, default=True, choices=BOOL_CHOICES, help=SHOW_HELP)

    #TODO:
    #COMPUTE POTENTIAL

    return parser


def main():
    from . import api

    if sys.version_info < (3, 6):
        sys.stderr.write("You need python 3.6 or later to run this program\n")
        sys.exit(1)

    start_time = time.time()
    cli_parser = get_cli_parser()
    args = cli_parser.parse_args()

    if args.command == 'eos':
        api.compute_eos(args.prefix, args.outdir, args.eos_type,
                        args.fileout, args.fileplot, args.show)
    elif args.command == 'bands':
        api.compute_band_structure(
            args.prefix, args.outdir, args.schema, args.reference_energy, args.emin,
            args.emax, args.fileplot, args.show
        )
    elif args.command == 'dos':
        if args.window is None:
            window = None
        else:
            window = args.window
        api.compute_dos(
            args.prefix, args.outdir, args.schema, args.width, window, args.npts, args.fileout,
            args.fileplot, args.show
        )
    elif args.command == 'charge':
        api.compute_charge(
            args.prefix, args.outdir, args.schema, args.fileout, args.x0, args.e1, args.nx,
            args.e2, args.ny, args.e3, args.nz, args.radius, args.dim, args.ifmagn, args.plot_file,
            args.method, args.format, args.fileplot, args.show
        )
        # api.new_get_charge(args.prefix, args.outdir, args.filplot) #TODO: new_get_charge method
    elif args.command == 'potential':
        api.compute_potential(
            args.prefix, args.outdir, args.schema, args.pot_type, args.fileout, args.x0,
            args.e1, args.nx, args.e2, args.ny, args.e3, args.nz, args.radius, args.dim,
            args.exportfile, args.method, args.format, args.show
        )
    else:
        print('Command not implemented! Exiting...')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished. Elapsed time: " + str(elapsed_time) + " s.")
