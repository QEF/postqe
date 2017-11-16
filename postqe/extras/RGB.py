#color calculator

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


#variables and files for rgb calculations

lambdas_col = np.linspace(380, 780, 81) 

SunShine_read=np.loadtxt("D65.light")
CIEx_read=np.loadtxt("CIEx")
CIEy_read=np.loadtxt("CIEy")
CIEz_read=np.loadtxt("CIEz")


SunShine=SunShine_read[:,1]
CIEx=CIEx_read[:,1]
CIEy=CIEy_read[:,1]
CIEz=CIEz_read[:,1]

#functions

#transmitted spectra
def calcTrasmittedSpec(Io, Spec, x):
        max=Spec.max(axis=1)
        Itmp = np.exp(-(Spec*x)/max[:,np.newaxis])
        I = Itmp*Io
        return I

#sun RGB 
def RGB_SS():
    X = np.sum(SunShine * CIEx)
    Y = np.sum(SunShine * CIEy)
    Z = np.sum(SunShine * CIEz)
    tmp=XYZColor(X*100/Y,Y*100/Y,Z*100/Y)
    rgb=convert_color(tmp,sRGBColor)
    return rgb   

#XYZ and RGB codes

def XYZ(Spec, p):
    I = calcTrasmittedSpec(SunShine, Spec, p)
    X = np.sum(I * CIEx,axis=1)
    Y = np.sum(I * CIEy,axis=1)
    Z = np.sum(I * CIEz,axis=1)
    XYZ=np.vstack((X,Y,Z))
    nomrY_SS=100/np.sum(SunShine * CIEy)   #normalized, sunlight Y parameter must be 100
    return np.transpose(XYZ)*nomrY_SS

   

def RGB(Spec,griglia,depth):
    S_col=ipl.interp1d(griglia,Spec,kind='linear')(lambdas_col)
    zxy=XYZ(S_col, depth)
    rgb=np.empty(zxy.shape)
    rgb_ss=RGB_SS()
    for ii, xyz in enumerate(zxy):
        tmp=XYZColor(xyz[0],xyz[1],xyz[2])
        tmp2=convert_color(tmp,sRGBColor)
        rgb[ii,:]=np.array([tmp2.rgb_r/rgb_ss.rgb_r,tmp2.rgb_g/rgb_ss.rgb_g,tmp2.rgb_b/rgb_ss.rgb_b])
    return rgb










    

    
