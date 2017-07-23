#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":

    from postqe import fitEtotV, plot_EV

    fin = "./EtotV.dat"  		            # file with the total energy data E(V)
    V, E, a, chi2 = fitEtotV(fin)    	    # fits the E(V) data, returns the coefficients a and
                                    	    # the chi squared chi2
    
    fig1 = plot_EV(V,E,a)                  	# plot the E(V) data and the fitting line
    fig1.savefig("figure_1.png")            # save the matplotlib figure

