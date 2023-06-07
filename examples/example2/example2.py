#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an example showing how to compute the the band structure of silicon
using the "compute_bands_structure" API

OPTIONS:
    prefix: prefix of saved output files
    outdir: directory containing the input data
    schema: the XML schema to be used to read and validate the XML output file
    reference_energy: the Fermi level, defines the zero of the plot along y axis
    emin: the minimum energy for the band plot (eV)
    emax: the maximum energy for the band plot (eV)
    fileplot: output plot file (default='bandsplot.png') in png format.
    show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True

Returns: an ASE band structure object (bands_structure) and a Matplotlib figure object (figure)
"""

if __name__ == "__main__":
    from postqe import compute_band_structure
    print("Plotting bands of silicon and saving them in Si_bandsplot.png")

    si_bands_structure, si_figure = compute_band_structure(
            prefix='Si', outdir='si/out/bands', reference_energy=6.1423, emin=-10,
            emax=20, fileplot='Si_bandsplot.png', show=True)

    print("Plotting bands of Cobalt  and saving them in Co_bandsplot.png")
  
    compute_band_structure, co_figure = compute_band_structure(
            prefix='co', outdir='./co/out/bands/', reference_energy=18.3568,
            emin=8, emax=24, fileplot='Co_bandsplot.png', show=True)
