#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of plotting a 1D section (dim=1) of the electronic charge density
using the 'compute_charge' API.

HDF5 charge file must be present, to enable HDF5 support in QE see:
    https://www.quantum-espresso.org/Doc/user_guide/node12.html#SECTION00035050000000000000,
    and look up the dmft option in pw.x:
    https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm525)

OPTIONS:
    prefix: prefix of saved output file
    outdir: directory containing the input data
    schema: the XML schema to be used to read and validate the XML output file
    fileout: text file with the full charge data as in the HDF5 file
    x0: 3D vector (a tuple), origin of the line
    e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
    nx, ny, nz: number of points along e1, e2, e3
    radius: radious of the sphere in the polar average method
    dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
    ifmagn: for a magnetic calculation, 'total' plot the total charge, 'up' plot the charge with spin up, 'down' for spin down.
    plot_file: file where plot data are exported in the chosen format (Gnuplot, XSF, cube, Gaussian, etc.).
    method: interpolation method. Available choices are:
                    'FFT' -> Fourier interpolation (default)
                    'polar' -> 2D polar plot on a sphere
                    'spherical' -> 1D plot of the spherical average
                    'splines' -> not implemented
    format: format of the (optional) exported file. Available choices are:
                    'gnuplot' -> plain text format for Gnuplot (default). Available for 1D and 2D sections
                    'xsf' -> XSF format for the XCrySDen program. Available for 2D and 3D sections.
                    'cube' -> cube Gaussian format. Available for 3D sections
                    'contour' -> format for the contour.x code of Quantum Espresso
                    'plotrho' -> format for the plotrho.x code of Quantum Espresso
    show: if True, show the Matplotlib plot (only for 1D and 2D sections)

Returns: a Charge object and a Matplotlib figure object for 1D and 2D sections (figure), a Charge object and None for 3D sections (charge)
"""

if __name__ == "__main__":
    from postqe import compute_charge

    charge, figure = compute_charge(prefix='Silicon', outdir='outdir', schema='../../schemas/qes.xsd', fileout='charge.out',
                            x0=(0., 0., 0.), e1=(2., 0., 0.), nx=100, radius=1, dim=1,
                                ifmagn='total', plot_file='charge_plot', method='FFT', format='gnuplot', show=False)

    figure.savefig('1D_charge_plot.png', format='png')
