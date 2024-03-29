#!/bin/bash

### COMMENT
#     This is an example showing how to compute the bands structure plot using the postqe CLI
#
#     1) postqe general parameters (for all commands):
#         -prefix: prefix of files saved by program pw.x (default: pwscf)
#         -outdir: directory containing the input data, i.e. the same as in pw.x
#         -fileplot: output plot file (default='plot_file.png') in png format
#         -schema: the XML schema to be used to read and validate the XML output file
#         -show: plot results with Matplotlib (True, False)
#
#     2) specific parameters for the dos command:
#         -reference_energy: the Fermi level, defines the zero of the plot along y axis (default=0)
#         -emin: the minimum energy for the band plot (default=-50)
#         -emax: the maximum energy for the band plot (default=50)

postqe -prefix='Si' -outdir='si/out/bands' -fileplot='CLI_bandsPLOT' bands -reference_energy=6.1423 -emin=-10 -emax=20
#
echo "running magnetic example: Hexagonal Co" 
postqe -prefix='co' -outdir='co/out/bands' -fileplot='cobalt_CLI_bandsPLOT' bands -reference_energy=18.3568  -emin=7.4 -emax=24
