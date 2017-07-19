#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":

    from postqe import plot_charge1D

    fin = "./Ni.xml"  				# file xml produce by QE
    plot_charge1D(fin)    			# plot a 1D section of the charge 
                                    	    
    #fig1 = plot_EV(V,E,a)                  	# plot the E(V) data and the fitting line
    #fig1.savefig("figure_1.png")

