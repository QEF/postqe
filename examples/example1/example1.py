#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":

    from postqe import get_charge, get_cell_data, compute_G, plot1D_FFTinterp

    fin = "./Ni.xml"  				            # file xml produce by QE
    ibrav, alat, a, b = get_cell_data(fin)  # get some data on the unit cell
    charge, chargediff = get_charge(fin)    	# get the charge (and charge diff)

    G = compute_G(b, charge.shape)
    fig1 = plot1D_FFTinterp(charge, G, a, plot_file='plotfile')
    fig1.show()
    fig1.savefig("figure_1.png")




