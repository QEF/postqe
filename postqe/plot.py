#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

"""
Utility functions for Fourier and other interpolation methods and for plotting the charge.
Some wrappers for the *matplotlib* functions.

.. Notes::\n
Plotting functions return a *matplotlib* figure object which can be modified by the user.\n
Faster Cython functions are available for FFT interpolation in 'cythonextensions' directory.
"""

################################################################################

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from .writecharge import write_1Dcharge_file, write_2Dcharge_file, write_3Dcharge_file
from .constants import pi


def FFTinterp1D(charge, G, a, x0, e1, nx):

        # normalize e1
        m1 = np.linalg.norm(e1)
        if abs(m1) < 1.0E-6:  # if the module is less than 1.0E-6
            e1 = a[0]
            m1 = np.linalg.norm(e1)
        e1 = e1 / m1

        # Computes the FFT of the charge
        fft_charge = np.fft.fftn(charge)
        nr = charge.shape

        # Steps along the e1 direction...
        deltax = m1 / (nx - 1)
        X = np.zeros(nx)
        Y = np.zeros(nx, dtype=complex)

        for i in range(0, nx):
            xi = x0[0] + i * deltax * e1[0]
            yi = x0[1] + i * deltax * e1[1]
            zi = x0[2] + i * deltax * e1[2]

            # For each point, evaluate the charge by Fourier interpolation
            for x in range(0, nr[0]):
                for y in range(0, nr[1]):
                    for z in range(0, nr[2]):
                        arg = 2.0 * pi * (xi * G[x, y, z, 0] + yi * G[x, y, z, 1] + zi * G[x, y, z, 2])
                        Y[i] += fft_charge[x, y, z] * complex(np.cos(arg), np.sin(arg))

            X[i] = i * deltax
            Y[i] = Y[i] / (nr[0] * nr[1] * nr[2])
            #print(X[i], Y[i].real)

        return X, Y

def spherical1D(charge, G, a, x0, e1, nx):
    """
    This function calculates a 1D plot of the input charge (or else), starting from the
    input point x0 and along the direction given by the vector e1. The G vectors
    in the reciprocal space must also be given in input. nx is the number of
    points where the spherically averaged charge is calculated (rho0(|r|) = int rho(r) dOmega
    rho0(r) = 4pi / sum_G rho(G) j_0(|G||r|)).

    :param charge:  eletronic charge density (or other quantity) to be plotted
    :param G:  G vectors in the reciprocal space
    :param a:  basis vectors of the unit cell
    :param x0: 3D vector, origin of the line
    :param e1: 3D vector which determines the plotting line
    :param nx: number of points in the line
    """
    #TODO: this function needs some optimization...
    # normalize e1
    m1 = np.linalg.norm(e1)
    if abs(m1) < 1.0E-6:  # if the module is less than 1.0E-6
        e1 = a[0]
        m1 = np.linalg.norm(e1)
    #e1 = e1 / m1

    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    # Steps along the e1 direction...
    deltax = m1 / (nx - 1)
    X = np.zeros(nx)
    Y = np.zeros(nx, dtype=complex)

    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):
                if (x==0) and (y==0) and (z==0):    # at the Gamma point
                    if (np.linalg.norm(G[0, 0, 0])<1e-10):
                        for i in range(0, nx):
                            Y[i] += 4.0 * pi * fft_charge[0, 0, 0]
                else:                               # not at Gamma
                    arg = 2.0 * pi * np.dot(x0, G[x, y, z])
                    rho0g = fft_charge[x, y, z] * complex(np.cos(arg), np.sin(arg))  # Move the origin to x0
                    Y[0] += 4.0 * pi * rho0g        # term at r=0
                    for i in range(1, nx):          # other terms at r!=0
                        gr = 2.0 * pi * np.linalg.norm(G[x,y,z]) * i * deltax
                        Y[i] += 4.0 * pi * rho0g * np.sin(gr) / gr

    for i in range(0, nx):
        X[i] = i * deltax                           # fill
        Y[i] = Y[i] / (nr[0] * nr[1] * nr[2])       # normalize

    return X, Y


def FFTinterp2D(charge, G, a, x0, e1, e2, nx, ny):
    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1) < 1.0E-6):  # if the module is less than 1.0E-6
        e1 = a[0]
        m1 = np.linalg.norm(e1)
    e1 = e1 / m1

    # normalize e2
    m2 = np.linalg.norm(e2)
    if abs(m2) < 1.0E-6:  # if the module is less than 1.0E-6
        e2 = a[1]
        m2 = np.linalg.norm(e2)
    e2 = e2 / m2

    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    # Steps along the e1 and e2 directions...
    deltax = m1 / (nx - 1)
    deltay = m2 / (ny - 1)

    X = np.zeros((nx, ny))
    Y = np.zeros((nx, ny))
    Z = np.zeros((nx, ny), dtype=complex)

    for i in range(0, nx):
        for j in range(0, ny):
            X[i, j] = i * deltax
            Y[i, j] = j * deltay

    # loop(s) over the G points
    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):

                # eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
                # compute these factors to save CPU time
                eigx = np.zeros(nx, dtype=complex)
                for i in range(0, nx):
                    np.exp(np.complex(0, 2.0 * pi * i * deltax * np.dot(e1, G[x, y, z])) + np.dot(x0, G[x, y, z]))

                eigy = np.zeros(ny, dtype=complex)
                for j in range(0, ny):
                    eigy[j] = np.exp(np.complex(0, 2.0 * pi * j * deltay * np.dot(e2, G[x, y, z])))

                for i in range(0, nx):
                    for j in range(0, ny):
                        Z[i, j] += fft_charge[x, y, z] * eigx[i] * eigy[j]
    Z = Z / (nr[0] * nr[1] * nr[2])

    return X, Y, Z


def plot_1Dcharge(charge, G, struct_info, x0=(0, 0, 0), e1=(1, 0, 0), nx=20, ylab='charge', plot_file='',
                  method='FFT', format=''):
    """
    This function calculates a 1D plot of the input charge (or else), starting from the
    input point x0 and along the direction given by the vector e1. The G vectors
    in the reciprocal space must also be given in input. nx is the number of
    points where the charge is effectively calculated using Fourier interpolation.

    :param charge:  eletronic charge density (or other quantity) to be plotted
    :param G:  G vectors in the reciprocal space
    :param struct_info: a dictionary with structural info on the unit cell
    :param x0: 3D vector, origin of the line
    :param e1: 3D vector which determines the plotting line
    :param nx: number of points in the line
    :param ylab: y axix label in the plot ('charge', 'Vtot', etc.)
    :param plot_file: if plot_file!='', write the plotting values on a text file
    :param format: output file format (not used in 1D)
    :param method: interpolation method. 'spherical' is for averaged spherical, 'FFT' for Fourier interpolation
    :return: the matplotlib figure object
    """
    if (method == 'FFT'):
        try:
            from compute_vs_cython import FFTinterp1D_Cython
            X, Y = FFTinterp1D_Cython(charge, G, struct_info['a'], x0, e1, nx)
        except ImportError:
            X, Y = FFTinterp1D(charge, G, struct_info['a'], x0, e1, nx)
    elif (method == 'spherical'):
        try:
            #TODO: Cython function to be implemented
            from compute_vs_cython import spherical1D_Cython
            X, Y = spherical1D_Cython(charge, G, struct_info['a'], x0, e1, nx)
        except ImportError:
            X, Y = spherical1D(charge, G, struct_info['a'], x0, e1, nx)
    else:
        raise NotImplementedError

    if plot_file != '':
        write_1Dcharge_file(X, Y, nx, plot_file)

    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    fig = plt.figure()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.plot(X, np.real(Y), 'r')
    plt.show()

    return fig


def polar2D(charge, G, nx, ny, radius, alat):
    """
    Computes a polar plot on a sphere.
    """

    phi = np.zeros((nx, ny))
    theta = np.zeros((nx, ny))
    Z = np.zeros((nx, ny), dtype=complex)
    r = np.zeros((nx, ny, 3))

    radius /= alat
    deltax = 2.0 * pi / (nx - 1)
    deltay = pi / (ny - 1)

    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    # first compute r coordinates on the sphere
    for i in range(0, nx):
        for j in range(0, ny):
            phi[i,j] = i * deltax
            theta[i,j] = j * deltay
            r[i,j,0] = radius * np.sin(theta[i,j]) * np.cos(phi[i,j])
            r[i,j,1] = radius * np.sin(theta[i,j]) * np.sin(phi[i,j])
            r[i,j,2] = radius * np.cos(theta[i,j])

    # loop(s) over the G points
    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):
                for i in range(0, nx):
                    for j in range(0, ny):
                        eig =  np.exp(np.complex(0,2.0 * pi * np.dot(r[i,j], G[x,y,z])))
                        Z[i,j] += fft_charge[x, y, z] * eig

    Z = Z / (nr[0] * nr[1] * nr[2])

    return phi, theta, Z


def plot_2Dcharge(charge, G, struct_info, x0=(0, 0, 0), e1=(1, 0, 0), e2=(0, 1, 0), nx=20, ny=20, radius=1, zlab='charge', plot_file='',
                  method='FFT', format=''):
    """
    This function calculates a 2D plot of the input charge (or else). If method=='FFT', the plot starts from the
    input point x0 and along the directions given by the vectors e1, e2. These
    vectors define the section plane along which the plot is draw. If method=='polar', the plot is a polar one on a
    sphere of radius = 'radius'. The G vectors in the reciprocal space must also be given in input. nx and ny define the
    grid of points where the charge is effectively calculated using Fourier interpolation (or polar projection).

    :param charge:  eletronic charge density (or other quantity) to be plotted
    :param G:  G vectors in the reciprocal space
    :param struct_info: a dictionary with structural info on the unit cell
    :param x0: 3D vector, origin of the line
    :param e1, e2: 3D vectors which determines the plotting plane
    :param nx, ny: number of points along e1, e2 respectively
    :param zlab: y axix label in the plot
    :param plot_file: if plot_file!='', write the plotting values on a text file
    :param method: interpolation method ('FFT' or 'splines' or 'polar', only 'FFT' and 'polar' implemented)
    :param format: output file format
    :return: the matplotlib figure object
    """

    #TODO: check that e1 and e2 are orthonormal
    # if (abs(e1(1) * e2(1) + e1(2) * e2(2) + e1(3) * e2(3)) > 1e-6):
    # throw something

    if (method == 'FFT'):
        try:
            from compute_vs_cython import FFTinterp2D_Cython
            X, Y, Z = FFTinterp2D_Cython(charge, G, struct_info['a'], x0, e1, e2, nx, ny)
        except ImportError:
            X, Y, Z = FFTinterp2D(charge, G, struct_info['a'], x0, e1, e2, nx, ny)
    elif (method == 'polar'):
        try:
            #TODO: Cython function to be implemented
            from compute_vs_cython import polar2D_Cython
            X, Y, Z = polar2D_Cython(charge, G, nx, ny, radius, struct_info['alat'])
        except ImportError:
            X, Y, Z = polar2D(charge, G, nx, ny, radius, struct_info['alat'])
    elif (method == 'splines'):
        # TODO: to be implemented
        raise NotImplementedError
    else:
        raise NotImplementedError

    if plot_file != '':
        write_2Dcharge_file(X, Y, Z, struct_info, x0, e1, e2, nx, ny, plot_file, method, format)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.3)
    cset = ax.contour(X, Y, Z, zdir='z', offset=Z.min(), cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='x', offset=X.min(), cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='y', offset=Y.max(), cmap=cm.coolwarm)

    xlab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    xlab += "x*("+str(e1[0])+","+str(e1[1])+","+str(e1[2])+")"
    ylab = "("+str(x0[0])+","+str(x0[1])+","+str(x0[2])+") + "
    ylab += "y*("+str(e2[0])+","+str(e2[1])+","+str(e2[2])+")"
    ax.set_xlabel(xlab)
    ax.set_xlim(X.min(), X.max())
    ax.set_ylabel(ylab)
    ax.set_ylim(Y.min(), Y.max())
    ax.set_zlabel(zlab)
    ax.set_zlim(Z.min(), Z.max())
    plt.show()

    return fig


def FFTinterp3D(charge, G, a, x0, e1, e2, e3, nx, ny, nz):
    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1) < 1.0E-6):  # if the module is less than 1.0E-6
        e1 = a[0]
        m1 = np.linalg.norm(e1)
    e1 = e1 / m1

    # normalize e2
    m2 = np.linalg.norm(e2)
    if abs(m2) < 1.0E-6:  # if the module is less than 1.0E-6
        e2 = a[1]
        m2 = np.linalg.norm(e2)
    e2 = e2 / m2

    # normalize e3
    m3 = np.linalg.norm(e3)
    if abs(m3) < 1.0E-6:  # if the module is less than 1.0E-6
        e3 = a[2]
        m3 = np.linalg.norm(e3)
    e3 = e3 / m3

    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    # Steps along the e1, e2 and e3 directions...
    deltax = m1 / (nx - 1)
    deltay = m2 / (ny - 1)
    deltaz = m3 / (nz - 1)

    X = np.zeros((nx, ny, nz))
    Y = np.zeros((nx, ny, nz))
    Z = np.zeros((nx, ny, nz))
    W = np.zeros((nx, ny, nz), dtype=complex)

    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                X[i, j, k] = i * deltax
                Y[i, j, k] = j * deltay
                Z[i, j, k] = k * deltaz

    # loop(s) over the G points
    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):

                # eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
                # compute these factors to save CPU time
                eigx = np.zeros(nx, dtype=complex)
                for i in range(0, nx):
                    eigx[i] = np.exp( np.complex(0, 2.0 * pi * i * deltax * np.dot(e1,G[x, y, z])) + np.dot(x0,G[x, y, z]))

                eigy = np.zeros(ny, dtype=complex)
                for j in range(0, ny):
                    eigy[j] = np.exp( np.complex(0, 2.0 * pi * j * deltay * np.dot(e2,G[x, y, z])))

                eigz = np.zeros(nz, dtype=complex)
                for k in range(0, nz):
                    eigz[k] = np.exp( np.complex(0, 2.0 * pi * k * deltaz * np.dot(e3,G[x, y, z])))

                for i in range(0, nx):
                    for j in range(0, ny):
                        for k in range(0, nz):
                            W[i, j, k] += fft_charge[x, y, z] * eigx[i] * eigy[j] * eigz[k]

    W = W / (nr[0] * nr[1] * nr[2])

    return X, Y, Z, W


def plot_3Dcharge(charge, G, struct_info, x0=(0, 0, 0), e1=(1, 0, 0), e2=(0, 1, 0), e3=(0, 0, 1), nx=20, ny=20, nz=20,
                  zlab='charge', plot_file='', method='FFT', format=''):
    try:
        from compute_vs_cython import FFTinterp3D_Cython
        X, Y, Z, W = FFTinterp3D_Cython(charge, G, struct_info['a'], x0, e1, e2, e3, nx, ny, nz)
    except ImportError:
        X, Y, Z, W = FFTinterp3D(charge, G, struct_info['a'], x0, e1, e2, e3, nx, ny, nz)
    pass

    write_3Dcharge_file(X, Y, Z, W, struct_info, x0, e1, e2, e3, nx, ny, nz, plot_file, method, format)


def simple_plot_xy(x, y, xlabel="", ylabel=""):
    """
    This function generates a simple x:y plot with matplotlib.
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
    This function generates a simple x:y plot with matplotlib overlapping several
    lines as in the matrix y. y second index refers to a line in the plot, the first
    index is for the array to be plotted.
    """

    if len(y[0, :]) > 7:
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
