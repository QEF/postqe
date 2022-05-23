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

# SPIN_CHOICES = [0, 1, 2]
SPIN_CHOICES = ['total', 'up', 'down']
METHOD_CHOICES = ['FFT', 'polar', 'spherical']
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
        emin, emax = map(int, s.split(','))
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("Tuple must be: emin, emax")
    else:
        return emin, emax

def get_cli_parser():
    parser = argparse.ArgumentParser(description='QE post processing')
    subparsers = parser.add_subparsers(metavar="plot_num", dest='plot_num', required=True,
                                       help='Selects what to save in filplot')

    #COMPUTE CHARGE
    plot_num_parser = subparsers.add_parser(
        '0', aliases=['charge'], help="electron (pseudo-)charge density")
    #SPECIFIC OPTIONS
    plot_num_parser.add_argument(
        '-fileout', type=str, default='charge.dat', help="text file with the full charge data as in the HDF5 file")
    plot_num_parser.add_argument(
        '-ifmagn', type=str, default=None, choices=SPIN_CHOICES, help="spin component of charge")
    plot_num_parser.add_argument(
        '-x0', type=vector, default=(0.,0.,0.), help="vector (a tuple), origin of the line")
    plot_num_parser.add_argument(
        '-e1', type=vector, default=(1.,0.,0.), help="3D vector (tuples) which determines the plotting lines must be in  x,y,z  format")
    plot_num_parser.add_argument(
        '-nx', type=int, default=50, help="number of points along e1")
    plot_num_parser.add_argument(
        '-e2', type=vector, default=(0.,1.,0.), help="3D vector (tuples) which determines the plotting lines must be in '(x,y,z)' format")
    plot_num_parser.add_argument(
        '-ny', type=int, default=50, help="number of points along e2")
    plot_num_parser.add_argument(
        '-e3', type=vector, default=(0.,0.,1.), help="3D vector (tuples) which determines the plotting lines must be in '(x,y,z)' format")
    plot_num_parser.add_argument(
        '-nz', type=int, default=50, help="number of points along e3")
    plot_num_parser.add_argument(
        '-radius', type=int, default=1, help="radious of the sphere in the polar average method")
    plot_num_parser.add_argument(
        '-dim', type=int, default=1, choices=DIM_CHOICES, help="1, 2, or 3 for 1D, 2D or 3D section respectively")
    plot_num_parser.add_argument(
        '-plot_file', type=str, default='plot_data.dat', choices=DIM_CHOICES, help="file where plot data are exported in the chosen format \
                                                                                (Gnuplot, XSF, cube Gaussian, etc.).")
    plot_num_parser.add_argument(
        '-method', type=str, default='FFT', choices=METHOD_CHOICES, help="the interpolation method")
    plot_num_parser.add_argument(
        '-format', type=str, default='gnuplot', choices=FORMAT_CHOICES, help="format of the (optional) exported file")
    plot_num_parser.add_argument(
        '-show', type=bool, default=False, choices=BOOL_CHOICES, help="show the Matplotlib plot (only for 1D and 2D sections)")


    #COMPUTE POTENTIAL
    plot_num_parser = subparsers.add_parser(
        '1', help="total potential V_bare + V_H + V_xc")
    #SPECIFIC OPTION
    # plot_num_parser.add_argument(
    #     '-spin', type=int, default=0, choices=SPIN_CHOICES, help="spin component of potential")
    plot_num_parser = subparsers.add_parser(
        '0', aliases=['charge'], help="electron (pseudo-)charge density")

    subparsers.add_parser('2', help="local ionic potential V_bare")

    #COMPUTE DOS
    plot_num_parser = subparsers.add_parser(
        '3', help="local density of states at specific energy or grid of energies")
    #SPECIFIC OPTIONS
    plot_num_parser.add_argument(
        '-window', type=window, default=None, help="a tuple (emin, emax) that defines the minimum and maximum energies for the DOS")
    plot_num_parser.add_argument(
        '-width', type=float, default=0.5, help="width of the gaussian to be used for the DOS (in eV, default=0.5)")
    plot_num_parser.add_argument(
        '-npts', type=int, default=100, help="number of points of the DOS")
    plot_num_parser.add_argument(
        '-fileout', type=str, default='dos.dat', help="output file with DOS results (default='dos.dat', not written)")
    plot_num_parser.add_argument(
        '-fileplot', type=str, default='dosplot', help="output plot file (default='dosplot') in png format.")
    plot_num_parser.add_argument(
        '-show', type=bool, default=False, choices=BOOL_CHOICES, help="plot results with Matplotlib (True, False")


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
    parser.add_argument('-schema', type=str, default=None,
                        help="The XML schema to be used to read and validate the XML output file")

    return parser


def main():
    from . import api

    if sys.version_info < (3, 6):
        sys.stderr.write("You need python 3.6 or later to run this program\n")
        sys.exit(1)

    start_time = time.time()
    cli_parser = get_cli_parser()
    args = cli_parser.parse_args()

    # if args.commands == 'eos':
    #     api.compute_eos(args.prefix, args.outdir, args.eos_type,
    #                     args.fileout, args.fileplot, args.show)
    # elif args.commands == 'bands':
    #     api.compute_band_structure(
    #         args.prefix, args.outdir, args.schema, args.reference_energy, args.emin,
    #         args.emax, args.fileplot, args.show
    #     )
    # elif args.commands == 'dos':
    #     if args.emin is None or args.max is None:
    #         window = None
    #     else:
    #         window = (args.emin, args.emax)
    #     api.compute_dos(
    #         args.prefix, args.outdir, args.schema, args.width, window, args.npts, args.fileout,
    #         args.fileplot, args.show
    #     )
    # elif args.commands == 'charge':
    #     api.compute_charge(
    #         args.prefix, args.outdir, args.schema, args.fileout, args.x0, args.e1, args.nx,
    #         args.e2, args.ny, args.e3, args.nz, args.radius, args.dim, args.ifmagn, args.exportfile,
    #         args.method, args.format, args.show
    #     )
    #     api.new_get_charge(args.prefix, args.outdir, args.filplot)
    # elif args.commands == 'potential':
    #     api.compute_potential(
    #         args.prefix, args.outdir, args.schema, args.pot_type, args.fileout, args.x0,
    #         args.e1, args.nx, args.e2, args.ny, args.e3, args.nz, args.radius, args.dim,
    #         args.exportfile, args.method, args.format, args.show
    #     )
    if args.plot_num == '0':
        api.compute_charge(
            args.prefix, args.outdir, args.schema, args.fileout, args.x0, args.e1,
            args.nx, args.e2, args.ny, args.e3, args.nz, args.radius, args.dim, args.ifmagn,
            args.plot_file, args.method, args.format, args.show
        )

    # elif args.plot_num == '1':
    # elif args.plot_num == '2':
    elif args.plot_num == '3':
        if args.window is None:
            window = None
        else:
            window = args.window
        api.compute_dos(
            args.prefix, args.outdir, args.schema, args.width, window, args.npts,
            args.fileout, args.fileplot, args.show
        )
    # elif args.plot_num == '4':
    # elif args.plot_num == '5':
    # elif args.plot_num == '6':
    # elif args.plot_num == '7':
    # elif args.plot_num == '8':
    # elif args.plot_num == '9':
    # elif args.plot_num == '10':
    # elif args.plot_num == '11':
    # elif args.plot_num == '12':
    # elif args.plot_num == '13':
    # elif args.plot_num == '14':
    # elif args.plot_num == '15':
    # elif args.plot_num == '16':
    # elif args.plot_num == '17':
    # elif args.plot_num == '18':
    # elif args.plot_num == '19':
    # elif args.plot_num == '20':
    # elif args.plot_num == '21':
    # elif args.plot_num == '22':
    else:
        print('Command not implemented! Exiting...')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished. Elapsed time: " + str(elapsed_time) + " s.")
