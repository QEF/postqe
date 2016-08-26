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
import settings


def plot1D(charge,G,x0=[0,0,0],e1=[1,0,0],nx=20,ylab=0):
    """
    This function calculates a 1D plot of the input charge (or else), starting from the 
    input point x0 and along the direction given by the vector e1. The G vectors
    in the reciprocal space must also be given in input. nx is the number of 
    points where the charge is effectively calculated using Fourier interpolation.
    """
    
    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1)<1.0E-6):    # if the module is less than 1.0E-6
        e1 = a[1]
        m1 = np.linalg.norm(e1)        
    e1 = e1 / m1  
    
    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape
  
    # Steps along the e1 direction...
    deltax = m1 / (nx-1)
    toplot = np.zeros(nx,dtype=complex)
    xv = np.zeros(nx)
    
    try:
        wx
        progdlg = wx.ProgressDialog("Plotting...","Time remaining", nx,
        style=wx.PD_APP_MODAL | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
    except:
        pass
        
    for i in range(0,nx):
        xi = x0[0] + i * deltax * e1[0]
        yi = x0[1] + i * deltax * e1[1]
        zi = x0[2] + i * deltax * e1[2]
    
        # For each point, evaluate the charge by Fourier interpolation
        for x in range(0,nr[0]):
            for y in range(0,nr[1]):
                for z in range(0,nr[2]):
                    arg = 2.0 * pi * (xi*G[x,y,z,0] + yi*G[x,y,z,1]  + zi*G[x,y,z,2])
                    toplot[i] += fft_charge[x,y,z] * complex(cos(arg),sin(arg)) 
                   
        xv[i] = i*deltax 
        toplot[i] = toplot[i]/(nr[0]*nr[1]*nr[2])
        print (xv[i],toplot[i].real)
        
        try:
            res = progdlg.Update(i)
        except:
            pass
        
    try:
        progdlg.Destroy()
    except:
        pass
    
    # Plot with Matplotlib library
    try:
        wx
        from matplotlib import use
        use('WXAgg')
    except:
        pass
    import matplotlib.pyplot as plt

    ylabels = ['charge','Vbare','Vbare+VH','Vtot']
    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    fig = plt.figure()
    plt.xlabel(xlab)
    plt.ylabel(ylabels[ylab])
    plt.plot(xv, np.real(toplot), 'r')
    return fig
    
    
def plot2D(charge,G,x0=[0,0,0],e1=[1,0,0],e2=[1,0,0],nx=20,ny=20,zlab=0):
    """
    This function calculates a 2D plot of the input charge (or else), starting from the 
    input point x0 and along the directions given by the vectors e1, e2. These
    vectors define the section plane along which the plot is draw. The G vectors
    in the reciprocal space must also be given in input. nx is the number of 
    points where the charge is effectively calculated using Fourier interpolation.
    """
    
    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1)<1.0E-6):    # if the module is less than 1.0E-6
        e1 = a[1]
        m1 = np.linalg.norm(e1)        
    e1 = e1 / m1  
    
    # normalize e2
    m2 = np.linalg.norm(e2)
    if (abs(m2)<1.0E-6):    # if the module is less than 1.0E-6
        e2 = a[2]
        m2 = np.linalg.norm(e2)        
    e2 = e2 / m2  
    
    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape
  
    # Steps along the e1 and e2 directions...
    deltax = m1 / (nx-1)
    deltay = m2 / (ny-1)
    
    temp = np.zeros((nx,ny),dtype=complex)
    X = np.zeros((nx,ny))
    Y = np.zeros((nx,ny))   
    Z = np.zeros((nx,ny))
    
    # loop(s) over the G points
    for x in range(0,nr[0]):
        for y in range(0,nr[1]):
            for z in range(0,nr[2]):
   
                # eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
                # compute these factors to save CPU time
                eigx = np.zeros(nx,dtype=complex)
                for i in range(0,nx):
                    eigx[i] = cmath.exp( 2.0 * pi * complex(0.0,1.0) * ( i * deltax *\
                    (e1[0] * G[x,y,z,0] + e1[1] * G[x,y,z,1] + e1[2] * G[x,y,z,2]) +\
                    (x0[0] * G[x,y,z,0] + x0[1] * G[x,y,z,1] + x0[2] * G[x,y,z,2]))) 

                eigy = np.zeros(ny,dtype=complex)
                for j in range(0,ny):
                    eigy[j] = cmath.exp( 2.0 * pi * complex(0.0,1.0) * ( j * deltax *\
                    (e2[0] * G[x,y,z,0] + e2[1] * G[x,y,z,1] + e2[2] * G[x,y,z,2]))) 
                   
                for i in range(0,nx):
                    for j in range(0,ny):
                        temp[i,j] += fft_charge[x,y,z] * eigx[i] * eigy[j]
                        
    Z = temp.real/(nr[0]*nr[1]*nr[2])
    # loop again over nx,ny to normalize and print
    for i in range(0,nx):
        for j in range(0,ny): 
            X[i,j] = i * deltax
            Y[i,j] = j * deltay
            print (X[i,j],Y[i,j],Z[i,j])
    

    # Plot with Matplotlib library
    try:
        wx
        from matplotlib import use
        use('WXAgg')
    except:
        pass
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.3)
    cset = ax.contour(X, Y, Z, zdir='z', offset=Z.min(), cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='x', offset=X.min(), cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='y', offset=Y.max(), cmap=cm.coolwarm)

    zlabels = ['charge','Vbare','Vbare+VH','Vtot']
    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    ylab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    ylab += "y*("+str(e2[0])+","+str(e2[1])+","+str(e2[2])+")"
    ax.set_xlabel(xlab)
    ax.set_xlim(X.min(),X.max())
    ax.set_ylabel(ylab)
    ax.set_ylim(Y.min(),Y.max())
    ax.set_zlabel(zlabels[zlab])
    ax.set_zlim(Z.min(),Z.max())

    return fig
    
    
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
                    
    plotparser.add_argument('-iplot', type=int, nargs='?', default=1,
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
                        
    plotparser.add_argument('-e2_1', type=float, nargs='?', default=1,
                    help='1st component of the 3D vector determining the plotting plane')
    plotparser.add_argument('-e2_2', type=float, nargs='?', default=0,
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
    
    charge, nr = read_postqe_output_file(pars.filein)
    
    if (pars.iplot==1):   # Read the charge and write it in filpl       
        x0 = [pars.x0_1,pars.x0_2,pars.x0_3]
        e1 = [pars.e1_1,pars.e1_2,pars.e1_3]
        nx = pars.nx
        G = compute_G(b,charge.shape,ecutrho,alat)
        fig = plot1D(charge,G,x0,e1,nx)
        fig.show()
        
    elif (pars.iplot==2):
        x0 = [pars.x0_1,pars.x0_2,pars.x0_3]
        e1 = [pars.e1_1,pars.e1_2,pars.e1_3]
        nx = pars.nx
        e2 = [pars.e2_1,pars.e2_2,pars.e2_3]
        ny = pars.ny
        G = compute_G(b,charge.shape,ecutrho,alat)
        fig = plot2D(charge,G,x0,e1,e2,nx,ny)
        fig.show()
    else:
        print ("Not implemented yet.")
    
    input("Enter to exit")
    


  



