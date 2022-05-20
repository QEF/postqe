#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS using the "compute_eos" API

OPTIONS:
    prefix: name of the input file with volumes and energies
    outdir: directory containing the input data
    eos_type: type of (EOS), murnaghan, sjeos, taylor, vinet, birch, birchmurnaghan, pouriertarantola, antonschmidt
    fileout: output file with fitting data and results, if not specified the file will not be written
    fileplot: output plot file (default='EOSplot') in png format
    show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    ax: a Matplotlib "Axes" instance (see Matplotlib documentation for details)

Returns: an QEEquationOfState object (eos) and a Matplotlib figure object (eos_plot)
"""

if __name__ == "__main__":
    from postqe import compute_eos

    eos, eos_plot = compute_eos(prefix='Nienergies.dat', outdir='outdir',eos_type='murnaghan',
                                        fileout='eos.out', fileplot='EOSplot', show=True, ax=None)

    fig = eos.plot()
    # Save the plot in a different formats (pdf) with Matplotlib if you like
    fig.figure.savefig('figure.jpg', format='jpg')
    fig.figure.savefig('figure.pdf', format='pdf')
