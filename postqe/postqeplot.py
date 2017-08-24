#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################

from constants import pi
from math import sin, cos
import cmath, time
import numpy as np
from postqe import get_from_xml
from celldm import calcola_celldm
from compute_vs import compute_G
from readutils import read_charge_file_iotk, read_postqe_output_file,\
write_charge, create_header
from plot import plot1D, plot2D
import settings
  
    
def get_plot_input_parameters():
    """
    Get postprocessing plot input parameters using argparse.
    """
    import argparse 
    
    plotparser = argparse.ArgumentParser(description="QE post processing plotting\n\
    This program reads an input file with charge or a potential data and computes\
    interpolated data along one direction or a plane using Fourier interpolation.")

    default_prefix = "Si"
    plotparser.add_argument('-prefix', type=str, nargs='?', default=default_prefix,
                    help='prefix of files saved by program pw.x')
                    
    plotparser.add_argument('-iplot', type=int, nargs='?', default=2,
        help="0 -> 1D plot of the spherical average (not implemented yet)\
            1 -> 1D plot\
            2 -> 2D plot\
            3 -> 3D plot (not implemented yet)\
            4 -> 2D polar plot on a sphere (not implemented yet)")
            
    plotparser.add_argument('-filein', type=str, nargs='?', default="Siout0",
                    help='name of the data file to be read')
        
    plotparser.add_argument('-fileout', type=str, nargs='?', default="fileout",
                    help='name of the file to which the plot is written')
                    
    plotparser.add_argument('-x0_1', type=float, nargs='?', default=0,
                    help='1st component of the 3D vector origin of the line or plane')
    plotparser.add_argument('-x0_2', type=float, nargs='?', default=0,
                    help='2nd component of the 3D vector origin of the line or plane')
    plotparser.add_argument('-x0_3', type=float, nargs='?', default=0,
                    help='3rd component of the 3D vector origin of the line or plane')
                    
    plotparser.add_argument('-e1_1', type=float, nargs='?', default=1,
                    help='1st component of the 3D vector determining the plotting line')
    plotparser.add_argument('-e1_2', type=float, nargs='?', default=0,
                    help='2nd component of the 3D vector determining the plotting line')
    plotparser.add_argument('-e1_3', type=float, nargs='?', default=0,
                    help='3rd component of the 3D vector determining the plotting line')
                    
    plotparser.add_argument('-nx', type=int, nargs='?', default=10,
                    help='number of points in the line:\
                        rho(i) = rho( x0 + e1 * (i-1)/(nx-1) ), i=1, nx')
                        
    plotparser.add_argument('-e2_1', type=float, nargs='?', default=0,
                    help='1st component of the 3D vector determining the plotting plane')
    plotparser.add_argument('-e2_2', type=float, nargs='?', default=1,
                    help='2nd component of the 3D vector determining the plotting plane')
    plotparser.add_argument('-e2_3', type=float, nargs='?', default=0,
                    help='3rd component of the 3D vector determining the plotting plane')
                        
    plotparser.add_argument('-ny', type=int, nargs='?', default=10,
                    help='number of points in the line:\
                        rho(i) = rho( x0 + e1 * (i-1)/(ny-1) ), i=1, ny')
                    
    args = plotparser.parse_args()
    
    return args

################################################################################
#   MAIN
################################################################################

if __name__ == "__main__":
    
    start_time = time.time()
    
    # get the input parameters
    pars = get_plot_input_parameters()
    print (pars)
    
    # get some needed values from the xml output
    ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species,\
    nat, ntyp = get_from_xml(pars.prefix+".xml",settings.schema)    
    celldms = calcola_celldm(alat,a[0],a[1],a[2],ibrav)
    
    charge = read_postqe_output_file(pars.filein)
    
    if (pars.iplot==1):   # Read the charge and write it in filpl       
        x0 = [pars.x0_1,pars.x0_2,pars.x0_3]
        e1 = [pars.e1_1,pars.e1_2,pars.e1_3]
        nx = pars.nx
        G = compute_G(b,charge.shape,ecutrho,alat)
        fig = plot1D(charge,G,a,x0,e1,nx)
        fig.show()
        
    elif (pars.iplot==2):
        x0 = [pars.x0_1,pars.x0_2,pars.x0_3]
        e1 = [pars.e1_1,pars.e1_2,pars.e1_3]
        nx = pars.nx
        e2 = [pars.e2_1,pars.e2_2,pars.e2_3]
        ny = pars.ny
        G = compute_G(b,charge.shape,ecutrho,alat)
        fig = plot2D(charge,G,a,x0,e1,e2,nx,ny)
        fig.show()
    else:
        print ("Not implemented yet.")
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")
    
    input("Enter to exit")
    

    


  



