#!/usr/bin/env python3
#encoding: UTF-8

"""
This submodule groups for fitting the energy using the Murnaghan EOS. It will be substituted by ASE EOS class.
"""

from math import pow
import numpy as np
from scipy.optimize import curve_fit

from .constants import RY_KBAR
from .readutils import read_EtotV

################################################################################
# Murnaghan EOS functions 
#
# This one is in the format ideal for fitting, not the canonical one in textbooks 
def E_MurnV(V,a0,a1,a2,a3):
    """
    This function implements the Murnaghan EOS (in a form which is best for fitting).
    Returns the energy at the volume *V* using the coefficients *a0,a1,a2,a3* 
    from the equation:
    
    .. math::
       E = a_0 - (a_2*a_1)/(a_3-1.0) V a_2/a_3 ( a_1/V^{a_3})/(a_3-1.0) +1.0 )
    
    """
    res=np.zeros(len(V))
    for i in range(0,len(V)):
        res[i]=a0 - a2*a1/(a3-1.0) + V[i]*a2/a3*( pow(a1/V[i],a3)/(a3-1.0)+1.0 )
    return res

# Other functions
def E_Murn(V,a):
    """
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]*.
    """
    return a[0] - a[2]*a[1]/(a[3]-1.0) + V*a[2]/a[3]*( pow(a[1]/V,a[3])/(a[3]-1.0)+1.0 )

def P_Murn(V,a):
    """
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]* and it returns the pressure not the energy from the EOS.
    """
    return a[2]/a[3]*(pow(a[1]/V,a[3])-1.0)

def H_Murn(V,a):
    """ 
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]* and it returns the enthalpy not the energy from the EOS.
    """
    return E_Murn(V,a)+P_Murn(V,a)*V


################################################################################

def print_eos_data(x,y,a,chi,ylabel="Etot"):   
    """
    Print the data and the fitted results using the Murnaghan EOS. It can be used for
    different fitted quantities using the proper ylabel. ylabel can be "Etot", 
    "Fvib", etc.
    """
    print ("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    print ("# "+ylabel+"min= {:.10e} Ry".format(a[0])+"\t Vmin= {:.10e} a.u.^3".format(a[1])+"\t B0= {:.10e} kbar".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    print (80*"#")
    print ("# V (a.u.^3)","\t\t",ylabel," (Ry)\t\t",ylabel+"fit"," (Ry)\t\t",ylabel+"-"+ylabel+"fit (Ry)\tP (kbar)")
    for i in range(0,len(y)):
        print ("{:.10e}".format(x[i]),"\t", "{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################

def write_Etotfitted(filename,x,y,a,chi,ylabel="E"): 
    """
    Write in filename the data and the fitted results using the Murnaghan EOS. It can be used for
    different fitted quantities using the proper ylabel. ylabel can be "Etot", 
    "Fvib", etc.
    """
    fout=open(filename, "w")
    fout.write("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    fout.write("# E0= {:.10e} Ry".format(a[1])+"\t V0= {:.10e} a.u.^3".format(a[1])+"\t B0= {:.10e} kbar".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    fout.write(80*"#")
    print ("# V *a.u.^3)","\t\t",ylabel," (Ry)\t\t",ylabel+"fit"," (Ry)\t\t",ylabel+"-"+ylabel+"fit (Ry)\tP (kbar)")
    for i in range(0,len(y)):
        fout.write("{:.10e}".format(x[i])+"\t"+"{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################
#
def calculate_fitted_points(V,a):
    """
    Calculates a denser mesh of E(V) points (1000) for plotting.
    """
    Vstep = (V[len(V)-1]-V[0])/1000
    Vdense = np.zeros(1000)
    Edensefitted = np.zeros(1000)
    for i in range(0,1000):
        Vdense[i] = V[0] + Vstep*i 
        Edensefitted[i] = E_Murn(Vdense[i],a)
        
    return Vdense, Edensefitted


################################################################################

def fit_Murn(V,E,guess=[0.0,0.0,900/RY_KBAR,1.15],lm_pars={}):
    """
    This is the function for fitting with the Murnaghan EOS as a function of volume only.

    The input variable *V* is an 1D array of volumes, *E* are the corresponding 
    energies (or other analogous quantity to be fitted with the Murnaghan EOS.
    *guess* (optional) is a list of 4 floats with initial guessess for *E_{min}*, *V_{min}*,
    *B_T* and *B_T'*. *lm_pars* is an optional dictionary of parameters for the scipy
    curve_fit routine (see related documentation).
    
    Note: volumes must be in a.u.^3 and energies in Rydberg.
    
    """
    # reasonable initial guesses for EOS parameters
    if guess[0]==0.0:
        guess[0] = E[len(E) // 2]
    if guess[1]==0.0:
        guess[1] = V[len(V) // 2]

    a, pcov = curve_fit(E_MurnV, V, E, p0=guess, **lm_pars)
    
    chi = 0
    for i in range(0,len(V)):
        chi += (E[i]-E_Murn(V[i],a))**2
    
    return a, pcov, chi


def fitEtotV(fin, fout=None):
    """
    This function reads :math:`E(V)` data from the input file *fin*, fits them with a Murnaghan EOS,
    prints the results on the *stdout* and write them in the file "fout".
    It returns the volumes and energies read from the input file, the fitted coefficients
    of the EOS and the corresponding :math:`\chi^2`.
    """

    V, E = read_EtotV(fin)
    a, cov, chi = fit_Murn(V, E)
    print_eos_data(V, E, a, chi, "Etot")
    if (fout != None):
        write_Etotfitted(fout, V, E, a, chi, "Etot")

    return V, E, a, chi