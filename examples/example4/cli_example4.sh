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
#         -fileout: text file with the full charge data as in the HDF5 file (default: charge.dat)
#         -ifmagn: for a magnetic calculation, 'total' plot the total charge,'up' plot the charge with spin up, 'down' for spin down
#         -x0: vector (a tuple), origin of the line (default: (0,0,0))
#         -e1: 3D vector (tuples) which determines the plotting lines must be in  x,y,z format (default: (1,0,0))
#         -nx: number of points along e1 (default: 50)
#         -e2: 3D vector (tuples) which determines the plotting lines must be in x,y,z format default: (0,1,0))
#         -ny: number of points along e2 (default: 50)
#         -e3: 3D vector (tuples) which determines the plotting lines must be in x,y,z format (default: (0,0,1))
#         -nz: number of points along e3 (default: 50))
#         -radius: radious of the sphere in the polar average method (default: 1)
#         -dim: 1, 2, or 3 for 1D, 2D or 3D section respectively (default: 1)
#         -method: the interpolation method (default: FFT)
#         -format: format of the (optional) exported file (default: gnuplot)
#         -plot_file: file where plot data are exported in the chosen format (default: plot_data.dat)

postqe -prefix='Silicon' -outdir='outdir' -fileplot='CLI_1D_charge_plot' charge -dim=1 -e1=2,0,0 -nx=50 \
                -ifmagn='total' -plot_file='CLI_plot_datafile.dat' -fileout='CLI_charge.dat'
