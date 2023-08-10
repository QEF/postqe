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


#sys.path.insert(0, './RGBlib')
import eig2spectra as sp
import RGB as rgb
import graphic as gr

###################################################################################################
#take esprsso eigenvalues or nm data files and calculate RGB values.
###################################################################################################

parser = argparse.ArgumentParser(description='Calculate RGB of spectra, calculate average spectrum, plot spectra, make movie')
parser.add_argument('-eig', default=".eigen", help='path of eigenvalues files. Format: "path/name_file". Default: "*.eigen"')
parser.add_argument('-spec', default=0, help='path of spectra data files. Format: "path/name_file". Default:0')
parser.add_argument('-depth', default=2, type=int, help='color intensity for rgb calculation. Default=2')
parser.add_argument('-rgb', default="rgb", help='path of rgb data file. Format: "path/name_file". Default: "./rgb"')

parser.add_argument('-av', default=0, help='Compute average spectrum and its rgb code. Format: "path/name_file".  Default=0.')

parser.add_argument('-plot', default=0, choices=["spectra", "average", "all"], help='make plots of single spectra, average spectrum or both')

parser.add_argument('-outS', default="S", help='path for plot.eps files. Format: "path/name_file". Default: "./S"')

parser.add_argument('-movie', default=0, help='make movie. Only for -eig flag. Format: "path/name_file". Default=0')
parser.add_argument('-duration', type=int, default=10, help='movie duration in seconds. Default: 10')


####################################################################################################


arg=parser.parse_args()
format=".dat"

#read files in eigenvalues or .dat format

if arg.spec!=0:
     files=sp.nameFiles(arg.spec)
     griglia_nm=np.loadtxt(files[0],usecols=(0,))
     S = np.zeros(files.size*griglia_nm.size).reshape(files.size,griglia_nm.size)
     for i in np.arange(files.size):
         S[i,:]=np.loadtxt(files[i],usecols=(1,))
     if arg.av!=0:
         name_av=arg.av+format
         S_av=np.expand_dims(np.mean(S,axis=0),axis=0)
         tmp=np.transpose(np.array((griglia_nm, S_av)))
         np.savetxt(name_av,S_av)

elif arg.eig!=0:
    files=sp.nameFiles(arg.eig)
    griglia_nm=sp.lambdas_nm
    energy_matrix = sp.EIG(files)
    S=sp.optical_spectra_default(energy_matrix[:,:,0],energy_matrix[:,:,1])
    if arg.av!=0:
        name_av=arg.av+format
        S_av=sp.Saverage_nm(energy_matrix[:,:,0],energy_matrix[:,:,1])
        tmp=np.transpose(np.array((griglia_nm, S_av)))
        np.savetxt(name_av,S_av)

#calculate and save RGB

RGB=rgb.RGB(S,griglia_nm,arg.depth)
name_rgb=arg.rgb+format
np.savetxt(name_rgb,RGB)

#print spectra images in eps format 

if arg.plot=="spectra" or arg.plot=="all":
    for i in np.arange(S[:,0].size):
        fig=plt.figure(i,figsize=(12, 10),dpi=1200)
        gr.plot_sp(griglia_nm, S[i], RGB[i])
        name = "{}.eps"
        folder=arg.outS+name
        plt.savefig(folder.format(i))
        plt.close(fig)    


if arg.av!=0:
    RGB_av=rgb.RGB(S_av,griglia_nm,arg.depth)
    avrgb="_rgb"
    rgb_av_out=arg.av+avrgb+format
    np.savetxt(rgb_av_out,RGB_av)
    if arg.plot=="average" or arg.plot=="all":
        plt.figure(99999,figsize=(12, 10),dpi=1200)
        gr.plot_sp(griglia_nm, np.squeeze(S_av, axis=0), np.squeeze(RGB_av, axis=0))
        name = "_av.eps"
        folder=arg.outS+name
        plt.savefig(folder)

#make movie in mp4 format

if arg.movie!=0:
    points=int((arg.duration*gr.frps)/(files.size-1))+1
    if arg.eig!=0:
        nOmega=gr.npoints(files.size,energy_matrix[0,:,0].size,energy_matrix[:,:,0],points)
        nOsc=gr.npoints(files.size,energy_matrix[0,:,0].size,energy_matrix[:,:,1],points)
        nS_nm=sp.optical_spectra_default(nOmega, nOsc)
    
    elif arg.spec!=0:
       print ("non si puo'fare")

    nRGB=rgb.RGB(nS_nm,griglia_nm,arg.depth)
    


    def make_frame(t):
        plt.clf()
        fig=plt.figure(1,figsize=(12, 10),dpi=1200)
        gr.plot_movie(griglia_nm, nS_nm[t * gr.frps], nRGB[t * gr.frps])
        fig.canvas.draw()
        figure = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        figure = figure.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        return figure
    
    
    format_movie=".mp4"
    movie_out=arg.movie+format_movie
    animation=mpy.VideoClip(make_frame,duration=arg.duration)
    animation.write_videofile(movie_out,fps=gr.frps)


