#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of plotting a 2D section of the electronic charge density.
"""
    
if __name__ == "__main__":

    from postqe import plot_charge2D

    fin = "./Ni.xml"  				# file xml produce by QE
    plot_charge2D(fin)    			# plot a 1D section of the charge 
                                    	    
    #fig1 = plot_EV(V,E,a)                  	# plot the E(V) data and the fitting line
    #fig1.savefig("figure_1.png")

