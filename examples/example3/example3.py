#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an example showing how to compute the electronic density of states (DOS or edos).

Computes the DOS for silicon (file Silicon.xml) with the "compute_dos" API

OPTIONS:
    prefix: prefix of saved output files
    outdir: directory containing the input data. Default to the value
    schema: the XML schema to be used to read and validate the XML output file
    width: width of the gaussian to be used for the DOS (in eV)
    window: a couple of values (emin, emax) that defines the minimum and maximum energies for the DOS
    npts: number of points of the DOS
    fileout: output file with DOS results (default='', not written).
    fileplot: output plot file (default='dosplot') in png format.
    show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True

Return: a DOS object (dos) and a Matplotlib figure object (plot)
"""

if __name__ == "__main__":
    from postqe import compute_dos

    # get a DOS object
    dos, plot = compute_dos(prefix='Silicon', outdir='si/out/nscf', width=0.25, npts=200,
                                fileout='DOSresults', fileplot='DOSplot.png', show=True)

    # If you want, get the dos and energies for further processing
    # d = dos.get_dos()
    # e = dos.get_energies()

    # save DOS in a file
    #dos.write('DOS.out')

    # this is an example for a magnetic system
    
    dos, plot = compute_dos(prefix='co', outdir='co/out/nscf', width=0.1, 
                            npts=200, fileout='co_DOSresults', window=(-13.0, 5.0),
                            fileplot='co_DOSplot.png', show=True)


