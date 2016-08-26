#Take eigenvalues in espresso format and print 

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


#default parameters used in optical_spectra function

G                  = 1001                #punti della griglia
Wl                 = 0.01                #Larghezza Lorenziane in Ry
Ry2nm              = 91.124664           # Ry = 91.124664/nm
Emin               = 91.124664/780       #omega min in the spec
Emax               = 91.124664/380       #omega max in the spec

points_plot        = 401


lambdas_nm = np.linspace(380, 780, points_plot)




#identify file temporal order by name
def nameFiles(folder):
    files = natsorted(glob.glob(folder))
    files=np.asarray(files)
    return files
#read eig files
def EIG(name_files):
    E=list()
    for i in np.arange(name_files.size):
        tmp=np.loadtxt(name_files[i],usecols=(0,1))
        E.append(tmp)
    E=np.array(E)
    return E   

#create a Lorenzian centered on c, spread Gamma
def Lore_f(c,Gamma,a):
    diff=a[np.newaxis, np.newaxis,:]-c[:, : ,np.newaxis]
    return (Gamma) / ( diff**2 + (Gamma)**2 )


def x_ry(e_min,e_max,N):
    xaxis=np.linspace(e_min,e_max,N)
    return xaxis



#create a spectrum for each eigenvalues files in Ry

def S_ry(omega,osc,xaxis,Gamma):
    yaxis=Lore_f(omega,Gamma,xaxis)
    L=xaxis[np.newaxis,np.newaxis,:]*yaxis
    S=L*osc[:,:,np.newaxis]
    S_all=np.sum(S,axis=1)
    return S_all


#convert Ry spectra in nm
def ry2nm(Sx, Sy, nuovox):
    griglia_nm = Ry2nm/Sx
    griglia_nm=griglia_nm[::-1]
    Spectra=Sy[:,::-1]
    S_nm=ipl.interp1d(griglia_nm,Spectra,kind='linear')(nuovox)
    return S_nm


#from espresso eigenvalues to specra in nm, default parameters used (see beginning of file)
def optical_spectra_default(omega,osc):
    x_energy = x_ry(Emin, Emax, G)
    spectra_ry = S_ry(omega, osc,x_energy, Wl)
    spectra_nm = ry2nm(x_energy, spectra_ry, lambdas_nm)
    return spectra_nm


def Saverage_nm(omega,osc):
    x_energy = x_ry(Emin, Emax, G)
    spectra_ry = S_ry(omega,osc,x_energy, Wl)
    griglia_nm = Ry2nm/x_energy
    griglia_nm=griglia_nm[::-1]
    S_av=np.mean(spectra_ry,axis=0)
    S_av_rev=S_av[::-1]
    S_av_nm=np.expand_dims(ipl.interp1d(griglia_nm,S_av_rev)(lambdas_nm),axis=0)
    return S_av_nm


#from espresso eigenvalues to specra in nm, user parameters

def optical_spectra_user(files,e_min,e_max,N, Gamma, xaxis_nm):
    energy_matrix = EIG(files)
    x_energy = x_ry(e_min, e_max, N)
    spectra_ry = S_ry(energy_matrix[:,:,0], energy_matrix[:,:,1],x_energy, Gamma)
    spectra_nm = ry2nm(x_energy, spectra_ry, xaxis_nm)
    return spectra_nm






