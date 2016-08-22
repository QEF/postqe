#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################

from constants import pi
from math import sin, cos
import numpy as np
from postqe import get_from_xml
from celldm import calcola_celldm
from compute_vs import compute_G
from readutils import read_charge_file_iotk, write_charge, create_header


def plot1Dcharge(charge,G,x0=[0,0,0],e1=[1,0,0],nx=20):
    """
    This function calculates a 1D plot of the input charge, starting from the 
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
        print (xv[i],toplot[i])
    
    # Plot with Matplotlib library
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt

    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    fig = plt.figure()
    plt.xlabel(xlab)
    plt.ylabel("charge")
    plt.plot(xv, np.real(toplot), 'r')
    return fig
    
    

################################################################################
#   MAIN, only for testing
################################################################################

if __name__ == "__main__":
    
    ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species,\
    nat, ntyp = get_from_xml("scf.xml")    
    celldms = calcola_celldm(alat,a[0],a[1],a[2],ibrav)
    
    charge_file = "charge-density.dat"
    charge = read_charge_file_iotk(charge_file)
    
    import wx
    from wxPlot1DDialog import Plot1DDialog
    app = wx.PySimpleApp()
    dlg = Plot1DDialog()
    dlg.ShowModal()
    x0 = np.array([float(dlg.x0_0.GetValue()), float(dlg.x0_1.GetValue()), float(dlg.x0_2.GetValue())])
    e1 = np.array([float(dlg.e1_0.GetValue()), float(dlg.e1_1.GetValue()),float(dlg.e1_2.GetValue())])
    nx = int(dlg.nx.GetValue())
    dlg.Destroy()
    
    G = compute_G(b,charge.shape,ecutrho,alat)
    fig = plot1Dcharge(charge,G,x0,e1,nx)
    fig.show()
    
    app.MainLoop()

  



