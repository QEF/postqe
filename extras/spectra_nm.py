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
import sys
import argparse

###################################################################################################
#take esprsso eigenvalues files and save nm spectra data files, from 380 to 780 with default parameters written in 
#eig2spectra.py. For custom parameters sp.optical_spectra_user(..) must be used.
###################################################################################################

#sys.path.insert(0, './RGBlib')
import eig2spectra as sp



parser = argparse.ArgumentParser(description='From eigenvalues to spectra data files')
parser.add_argument('-eig', default="*.eigen", type=str, help='path of eigenvalues files. Format: "path/name_file ". Default: "*.eigen"')
parser.add_argument('-av', default=0, help='name of average spectrum. Format: "path/name_file". Default=0')
parser.add_argument('-out', default="S", help='path of spectra data files. Default: "S"')


###################################################################################################


arg=parser.parse_args()


files=sp.nameFiles(arg.eig)
griglia_nm=sp.lambdas_nm
energy_matrix = sp.EIG(files)
S=sp.optical_spectra_default(energy_matrix[:,:,0],energy_matrix[:,:,1])
    
for i in np.arange(S[:,0].size):
    nameS="{}.dat"
    out=arg.out+nameS
    tmp=np.transpose(np.array((griglia_nm, S[i])))
    np.savetxt(out.format(i), tmp)



if arg.av!=0:
   ext=".dat"
   out=arg.av+ext
   energy_matrix = sp.EIG(files)
   S_av=sp.Saverage_nm(energy_matrix[:,:,0],energy_matrix[:,:,1])
   S_av=np.squeeze(S_av,axis=0)
   avtmp=np.transpose(np.array((griglia_nm, S_av)))
   np.savetxt(out, avtmp)


