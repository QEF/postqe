import numpy as np
from natsort import natsorted
import glob
from numpy import linalg as la
import pylab
import matplotlib.pyplot as plt
import scipy 
from scipy import interpolate as ipl
import colormath
from colormath.color_objects import sRGBColor, XYZColor, AdobeRGBColor
from colormath.color_conversions import convert_color
from matplotlib import rc
import moviepy.editor as mpy
import eig2spectra as sp
import RGB as rgb

#variables for plots and movie

frps = 25

def plotSp(x,S,rgb):
    plt.plot(x, S,'-', color = rgb)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.fill_between(x,S,0.3,color = rgb)
    plt.xlabel('wavelength [nm]', fontsize=45, color='black')
    plt.ylabel('absorption', fontsize=45, color='black')
    plt.axis([380, 780, 0.30, 25])
    plt.tick_params(axis='x', which = 'major', length=50, width=10, labelsize='40', top = 'off', bottom = 'off')
    plt.xticks([450,600,750])
    plt.yticks([])

def plotMovie(x,S,rgb):
    plt.plot(x, S,'-', color = rgb)
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.fill_between(x,S,0.3,color = rgb)
    plt.xlabel('wavelength [nm]', fontsize=40, color='black')
    plt.ylabel('absorption', fontsize=40, color='black')
    plt.axis([380, 780, 0.30, 25])
    plt.tick_params(axis='x', which = 'major', length=50, width=10, labelsize='35', top = 'off', bottom = 'off')
    plt.xticks([450,600,750])
    plt.yticks([])

#interpolation function for movie

def nPoints(x,y,M,n):
    M_out=list()
    nDim=n*(x-1)
    for i in np.arange(1,x):
        M_out.append(np.transpose(np.array([np.linspace(a, b, n) for a, b in zip(M[i-1,:], M[i,:])])))
    M_out=np.asarray(M_out).reshape(nDim,y)
    return M_out    








