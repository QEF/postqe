#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A collection of wrappers for the *matplotlib* functions.

.. Note::
  All functions return a *matplotlib* which can be modified by the user.
"""

################################################################################

import numpy as np
from math import sin, cos
import cmath

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from .eos import calculate_fitted_points
from .constants import pi


def plotcharge1D(charge, G, a, x0=(0, 0, 0), e1=(1, 0, 0), nx=20, ylab=0):
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

    for i in range(0, nx):
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
        #print (xv[i],toplot[i].real)
        
    ylabels = ['charge','Vbare','Vbare+VH','Vtot']
    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    fig = plt.figure()
    plt.xlabel(xlab)
    plt.ylabel(ylabels[ylab])
    plt.plot(xv, np.real(toplot), 'r')
    return fig
    
    
def plotcharge2D(charge, G, a, x0=(0, 0, 0), e1=(1, 0, 0), e2=(1, 0, 0), nx=20, ny=20, zlab=0):
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


def simple_plot_xy(x, y, xlabel="", ylabel=""):
    """
    This function generates a simple xy plot with matplotlib.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    ax.plot(x, y, 'r')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()

    return fig


def multiple_plot_xy(x, y, xlabel="", ylabel="", labels=""):
    """
    This function generates a simple xy plot with matplotlib overlapping several
    lines as in the matrix y. y second index refers to a line in the plot, the first
    index is for the array to be plotted.
    """

    if (len(y[0, :]) > 7):
        print("Too many data on y axis!")
        return

    colors = ['k', 'r', 'b', 'g', 'c', 'm', 'y']

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    if (labels == ""):
        try:  # try if there are multiple data on x axis
            for i in range(0, len(y[0, :])):
                ax.plot(x[:, i], y[:, i], colors[i])
        except:  # if not use a single x axis
            for i in range(0, len(y[0, :])):
                ax.plot(x, y[:, i], colors[i])
    else:
        try:  # try if there are multiple data on x axis
            for i in range(0, len(y[0, :])):
                ax.plot(x[:, i], y[:, i], colors[i], label=labels[i])
        except:  # if not use a single x axis
            for i in range(0, len(y[0, :])):
                ax.plot(x, y[:, i], colors[i], label=labels[i])
        ax.legend()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()

    return fig


def plot_EV(V, E, a=[0., 0., 0., 0.], labely="Etot"):
    """
    This function plots with matplotlib E(V) data and if a is given it also plot
    the fitted results
    """

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    ax.plot(V, E, 'o', label=labely + " data", markersize=10)
    if (a.all() != 0):
        Vdense, Edensefitted = calculate_fitted_points(V, a)
        ax.plot(Vdense, Edensefitted, 'r', label='Fitted EOS')
    ax.legend()
    ax.set_xlabel('V (a.u.^3)')
    ax.set_ylabel('E (Ry)')
    plt.show()

    return fig

