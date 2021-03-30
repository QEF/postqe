#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
A collection of functions defining postqe API and exposed to the user.
"""
import os

from .eos import QEEquationOfState
from .dos import QEDOS
from .charge import Charge, Potential
from .readutils import read_EtotV
from .calculator import PostqeCalculator


def get_label(prefix='pwscf', outdir=None):
    """Gets a filepath for Quantum ESPRESSO environment."""
    if outdir is None:
        try:
            outdir = os.environ['ESPRESSO_TMPDIR']
        except KeyError:
            outdir = os.curdir

    label = os.path.join(outdir, '{}.save'.format(prefix))
    return label


## New CLI-API interfaces ###

def new_get_charge(prefix=None, output=None, filplot=None):
    plot = []
    # ...
    if filplot is not None:
        with open(filplot, 'bw') as f:
            f.write(plot)
    return plot


## New CLI-API interfaces ###

def get_plot(plot_num, filplot=None, prefix=None, output=None, **kwargs):
    """
    Get plot

    :param plot_num:
    :param filplot: if `None` returns the computed plot object, otherwise save \
    it to the specified file.
    :param prefix: name of the input file with volumes and energies
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set, or current directory ('.') otherwise
    :param kwargs: additional arguments for plot computation.
    """
    plot = []

    print(plot_num, prefix, output, filplot, kwargs)

    if plot_num == '0':
        spin = kwargs.get('spin', 0)
        pass
    elif plot_num == '1':
        pass
    elif plot_num == '2':
        pass
    elif plot_num == '3':
        pass
    elif plot_num == '4':
        pass
    elif plot_num == '5':
        pass
    elif plot_num == '6':
        pass
    elif plot_num == '7':
        pass
    elif plot_num == '8':
        pass
    elif plot_num == '9':
        pass
    elif plot_num == '10':
        pass
    elif plot_num == '11':
        pass
    elif plot_num == '12':
        pass
    elif plot_num == '13':
        pass
    elif plot_num == '17':
        pass
    elif plot_num == '18':
        pass
    elif plot_num == '19':
        pass
    elif plot_num == '20':
        pass
    elif plot_num == '21':
        pass
    elif plot_num == '22':
        pass
    else:
        raise ValueError("Wrong value {!r} for 'plot_num' argument".format(plot_num))

    if filplot is not None:
        with open(filplot, 'bw') as f:
            f.write(plot)
        return
    return plot


def get_eos(prefix, outdir=None, eos_type='murnaghan'):
    """
    Fits an Equation of state of type *eos* and returns an QEEquationOfState object.
    Different equation of states are available (see below).

    :param prefix: name of the input file with volumes and energies
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set, or current directory ('.') otherwise
    :param eos_type: type of equation of state (EOS) for fitting. Available types are:\n\
    'murnaghan' (default) -> Murnaghan EOS, PRB 28, 5480 (1983)\n\
    'sjeos' -> A third order inverse polynomial fit, PhysRevB.67.026103\n\
    \t\tE(V) = c_0 + c_1 t + c_2 t^2  + c_3 t^3 ,  t = V^(-1/3)\n\
    'taylor' -> A third order Taylor series expansion around the minimum volume\n\
    'vinet' -> Vinet EOS, PRB 70, 224107 \n\
    'birch' -> Birch EOS, Intermetallic compounds: Principles and Practice, \
    Vol I: Principles, p. 195\n\
    'birchmurnaghan' -> Birch-Murnaghan EOS, PRB 70, 224107\n\
    'pouriertarantola' -> Pourier-Tarantola EOS, PRB 70, 224107\n\
    'antonschmidt' -> Anton-Schmidt EOS, Intermetallics 11, 23 - 32(2003)\n\
    'p3' -> A third order inverse polynomial fit\n
    :return: an QEEquationOfState object
    """

    label = get_label(prefix, outdir)
    # Extract volumes and energies from the input file:
    volumes, energies = read_EtotV(label)
    # Create an object EquationOfState and fit with Murnaghan (or other) EOS
    eos = QEEquationOfState(volumes, energies, eos=eos_type)
    return eos


def compute_eos(prefix, outdir=None, eos_type='murnaghan', fileout='',
                fileplot='EOSplot', show=True, ax=None):
    """
    Fits an Equation of state of type *eos_type*, writes the results into *fileout* (optionally)
    and creates a Matplotlib figure. Different equation of states are available (see below).

    :param prefix: name of the input file with volumes and energies
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set, or current directory ('.') otherwise
    :param eos_type: type of equation of state (EOS) for fitting. Available types are: \
    'murnaghan' (default) -> Murnaghan EOS, PRB 28, 5480 (1983);\n \
    'sjeos' -> A third order inverse polynomial fit, PhysRevB.67.026103\n\
    \t\tE(V) = c_0 + c_1 t + c_2 t^2  + c_3 t^3 ,  t = V^(-1/3)\n\
    'taylor' -> A third order Taylor series expansion around the minimum volume\n\
    'vinet' -> Vinet EOS, PRB 70, 224107 \n\
    'birch' -> Birch EOS, Intermetallic compounds: Principles and Practice, \
    Vol I: Principles, p. 195\n\
    'birchmurnaghan' -> Birch-Murnaghan EOS, PRB 70, 224107\n\
    'pouriertarantola' -> Pourier-Tarantola EOS, PRB 70, 224107\n\
    'antonschmidt' -> Anton-Schmidt EOS, Intermetallics 11, 23 - 32(2003)\n\
    'p3' -> A third order inverse polynomial fit
    :param fileout: output file with fitting data and results (default='', not written).
    :param fileplot: output plot file (default='EOSplot') in png format.
    :param show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    :param ax: a Matplotlib "Axes" instance (see Matplotlib documentation for details). \
    If ax=None (default), creates a new one.
    :return: an QEEquationOfState object and a Matplotlib figure object
    """
    eos = get_eos(prefix, outdir, eos_type)
    eos.fit()
    if fileout != '':
        eos.write(fileout)

    return eos, eos.plot(fileplot, show=show, ax=ax)


def get_band_structure(prefix, outdir=None, schema=None, reference_energy=0):
    """
    This function returns a "band structure" object from an output xml Espresso file
    containing the results of a proper calculation along a path in the Brilluoin zone.

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param reference_energy: the Fermi level, defines the zero of the plot along y axis
    :return: an ASE band structure object
    """
    calc = PostqeCalculator(atoms=None, label=prefix, outdir=outdir, schema=schema)
    calc.read_results()

    atoms = calc.get_atoms_from_xml_output()
    atoms.calc = calc
    atoms.calc.read_results()

    return atoms.calc.band_structure(reference=reference_energy)


def compute_band_structure(prefix, outdir=None, schema=None, reference_energy=0,
                           emin=-50, emax=50, fileplot='bandsplot.png', show=True):
    """
    This function returns a "band structure" object from an output xml Espresso file
    containing the results of a proper calculation along a path in the Brilluoin zone.

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param reference_energy: the Fermi level, defines the zero of the plot along y axis
    :param emin: the minimum energy for the band plot (default=-50)
    :param emax: the maximum energy for the band plot (default=50)
    :param fileplot: output plot file (default='bandsplot.png') in png format.
    :param show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    :return: an ASE band structure object and a Matplotlib figure object
    """
    bs = get_band_structure(prefix, outdir, schema=schema, reference_energy=reference_energy)
    fig = bs.plot(emin=emin, emax=emax, show=show, filename=fileplot)

    return bs, fig


def get_dos(prefix, outdir=None, schema=None, width=0.01, window=None, npts=100):
    """
    This function returns an DOS object from an output xml Espresso file containing the
    results of a DOS calculation.

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise.
    :param schema: the XML schema to be used to read and validate the XML output file
    :param width: width of the gaussian to be used for the DOS (in eV)
    :param window: a tuple (emin, emax) that defines the minimum and maximum energies for the DOS
    :param npts: number of points of the DOS
    :return: a DOS object
    """
    label = get_label(prefix, outdir)
    calc = PostqeCalculator(atoms=None, label=label, schema=schema, outdir=outdir)
    calc.read_results()

    atoms = calc.get_atoms_from_xml_output()
    atoms.set_calculator(calc)
    atoms.calc.read_results()

    # Create a DOS object with width= eV and npts points
    dos = QEDOS(calc, width=width, window=window, npts=npts)

    return dos


def compute_dos(prefix, outdir=None, schema=None, width=0.01, window=None, npts=100,
                fileout='', fileplot='dosplot.png', show=True):
    """
    This function returns an DOS object from an output xml Espresso file containing the
    results of a DOS calculation.

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param width: width of the gaussian to be used for the DOS (in eV)
    :param window: a couple of values (emin, emax) that defines the minimum and maximum \
    energies for the DOS
    :param npts: number of points of the DOS
    :param fileout: output file with DOS results (default='', not written).
    :param fileplot: output plot file (default='dosplot') in png format.
    :param show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    :return: a DOS object and a Matplotlib figure object
    """
    # get a DOS object
    dos = get_dos(prefix, schema=schema, width=width, window=window, npts=npts)

    # save DOS in a file
    if fileout != '':
        dos.write('DOS.out')

    # get the dos and energies for further processing
    d = dos.get_dos()
    e = dos.get_energies()

    # Plot the DOS with Matplotlib...
    import matplotlib.pyplot as plt
    plt.plot(e, d)
    plt.xlabel('energy [eV]')
    plt.ylabel('DOS')
    plt.savefig(fileplot)
    if show:
        plt.show()

    return dos, plt


def get_charge(prefix, outdir=None, schema=None):
    """
    Returns an Charge object from an output xml Espresso file and the
    corresponding HDF5 charge file containing the results of a calculation.

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :return: a Charge object
    """
    label = get_label(prefix, outdir)
    calc = PostqeCalculator(atoms=None, label=label, schema=schema)
    calc.read_results()

    atoms = calc.get_atoms_from_xml_output()
    atoms.set_calculator(calc)
    atoms.calc.read_results()

    nr = calc.get_nr()
    charge_file = calc.label + ".save/charge-density.hdf5"

    charge = Charge(nr)
    charge.read(charge_file)
    charge.set_calculator(calc)

    return charge


def compute_charge(prefix, outdir=None, schema=None, fileout='', x0=(0., 0., 0.), e1=(1., 0., 0.),
                   nx=50, e2=(0., 1., 0.), ny=50, e3=(0., 0., 1.), nz=50, radius=1, dim=1,
                   ifmagn='total', plot_file='', method='FFT', format='gnuplot', show=True):
    """
    Returns an Charge object from an output xml Espresso file and the
    corresponding HDF5 charge file containing the results of a calculation.
    Returns also a Matplotlib figure object from a 1D or 2D section of the charge.
    It also (optionally) exports the charge (1, 2 or 3D section) in a text file according to
    different formats (XSF, cube, Gnuplot, etc.).

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param fileout: text file with the full charge data as in the HDF5 file. Default='', \
    nothing is written.
    :param x0: 3D vector (a tuple), origin of the line
    :param e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
    :param nx, ny, nz: number of points along e1, e2, e3
    :param radius: radious of the sphere in the polar average method
    :param dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
    :param ifmagn: for a magnetic calculation, 'total' plot the total charge, \
    'up' plot the charge with spin up, 'down' for spin down.
    :param plot_file: file where plot data are exported in the chosen format \
    (Gnuplot, XSF, cube Gaussian, etc.).
    :param method: interpolation method. Available choices are:\n\
                    'FFT' -> Fourier interpolation (default)\n\
                    'polar' -> 2D polar plot on a sphere\n\
                    'spherical' -> 1D plot of the spherical average\n\
                    'splines' -> not implemented
    :param format: format of the (optional) exported file. Available choices are:\n\
                    'gnuplot' -> plain text format for Gnuplot (default). Available \
                    for 1D and 2D sections.\n\
                    'xsf' -> XSF format for the XCrySDen program. Available for 2D \
                    and 3D sections.\n\
                    'cube' -> cube Gaussian format. Available for 3D sections.\n\
                    'contour' -> format for the contour.x code of Quantum Espresso.\n\
                    'plotrho' -> format for the plotrho.x code of Quantum Espresso.\n
    :param show: if True, show the Matplotlib plot (only for 1D and 2D sections)
    :return: a Charge object and a Matplotlib figure object for 1D and 2D sections, \
    a Charge object and None for 3D sections
    """

    charge = get_charge(prefix=prefix, outdir=outdir, schema=schema)
    if fileout != '':
        charge.write(fileout)

    figure = charge.plot(x0=x0, e1=e1, nx=nx, e2=e2, ny=ny, e3=e3, nz=nz,
                         radius=radius, dim=dim, ifmagn=ifmagn, plot_file=plot_file,
                         method=method, format=format, show=show)

    return charge, figure


def get_potential(prefix, outdir=None, schema=None, pot_type='v_tot'):
    """
    This function returns an Potential object from an output xml Espresso file and
    the corresponding HDF5 charge file containing the results of a calculation.
    The available potentials are the bare (pot_type='v_bare'), Hartree (pot_type='v_h'),
    exchange-correlation (pot_type='v_xc') and total (pot_type='v_tot').

    :param prefix: prefix of saved output files
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param pot_type: type of the Potential ('v_tot', ....)
    :return: a Potential object
    """
    label = get_label(prefix, outdir)
    calc = PostqeCalculator(atoms=None, label=label, schema=schema)
    calc.read_results()

    atoms = calc.get_atoms_from_xml_output()
    atoms.set_calculator(calc)
    atoms.calc.read_results()

    nr = calc.get_nr()
    charge_file = calc.label + ".save/charge-density.hdf5"

    potential = Potential(nr, pot_type=pot_type)
    potential.read(charge_file)
    potential.set_calculator(calc)
    potential.compute_potential()

    return potential


def compute_potential(prefix, outdir=None, schema=None, pot_type='v_tot', fileout='',
                      x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, e2 = (0., 1., 0.),
                      ny = 50, e3 = (0., 0., 1.), nz = 50, radius = 1, dim = 1,
                      plot_file = '', method = 'FFT', format = 'gnuplot', show = True):
    """
    Returns an Potential object from an output xml Espresso file and the
    corresponding HDF5 charge file containing the results of a calculation.
    Returns also a Matplotlib figure object from a 1D or 2D section of the charge.
    It also (optionally) exports the charge (1, 2 or 3D section) in a text file according to
    different formats (XSF, cube, Gnuplot, etc.).

    :param prefix: prefix of saved output file
    :param outdir: directory containing the input data. Default to the value of
            ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema: the XML schema to be used to read and validate the XML output file
    :param fileout: text file with the calculate potential data. Default='', nothing is written.
    :param x0: 3D vector (a tuple), origin of the line
    :param e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
    :param nx, ny, nz: number of points along e1, e2, e3
    :param radius: radious of the sphere in the polar average method
    :param dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
    :param plot_file: file where plot data are exported in the chosen format \
    (Gnuplot, XSF, cube Gaussian, etc.)
    :param method: interpolation method. Available choices are:\n\
                    'FFT' -> Fourier interpolation (default)\n\
                    'polar' -> 2D polar plot on a sphere\n\
                    'spherical' -> 1D plot of the spherical average\n\
                    'splines' -> not implemented
    :param format: format of the (optional) exported file. Available choices are:\n\
                    'gnuplot' -> plain text format for Gnuplot (default). Available \
                    for 1D and 2D sections.\n\
                    'xsf' -> XSF format for the XCrySDen program. Available for 2D\
                    and 3D sections.\n\
                    'cube' -> cube Gaussian format. Available for 3D sections.\n\
                    'contour' -> format for the contour.x code of Quantum Espresso.\n\
                    'plotrho' -> format for the plotrho.x code of Quantum Espresso.\n
    :param show: if True, show the Matplotlib plot (only for 1D and 2D sections)
    :return: a Potential object and a Matplotlib figure object for 1D and 2D sections, \
    a Potential object and None for 3D sections
    """

    potential = get_potential(prefix=prefix, outdir=outdir, schema=schema, pot_type=pot_type)
    if fileout != '':
        potential.write(fileout)

    figure = potential.plot(x0=x0, e1=e1, nx=nx, e2=e2, ny=ny, e3=e3, nz=nz, radius=radius, dim=dim,
                            plot_file=plot_file, method=method, format=format, show=show)

    return potential, figure
