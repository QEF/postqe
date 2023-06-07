#!/bin/bash

### COMMENT
#     This is an example showing how to compute the electronic density of states (DOS or edos)
#     plot using the postqe CLI
#
#     1) postqe general parameters (for all commands):
#         -prefix: prefix of files saved by program pw.x (default: pwscf)
#         -outdir: directory containing the input data, i.e. the same as in pw.x
#         -fileplot: output plot file (default='plot_file.png') in png format
#         -schema: the XML schema to be used to read and validate the XML output file
#         -show: plot results with Matplotlib (True, False)
#
#     2) specific parameters for the dos command:
#         -window: a tuple (emin, emax) that defines the minimum and maximum energies for the DOS
#         -width: width of the gaussian to be used for the DOS (in eV, default=0.5)
#         -npts: number of points of the DOS
#         -fileout: output file with DOS results (default='dos.out')


postqe -prefix='Silicon' -outdir='si/out/nscf' -fileplot='CLI_dosPLOT' -show=True dos -fileout='CLI_dosOUT' -width=0.25 -npts=300

postqe -prefix='co' -outdir='co/out/nscf' -fileplot='CLI_dos_co_PLOT' -show=True dos -fileout='CLI_co_dosOUT' -width=0.1 -npts=300 -window=-13,6.0
