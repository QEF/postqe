#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of plotting a 2D section of the electronic charge density.
"""
    
if __name__ == "__main__":

    from postqe import get_charge, get_cell_data, compute_G, plot2D_FFTinterp

    fin = "./Ni.xml"  				            # file xml produce by QE
    ibrav, alat, a, b, nat, ntyp,\
    atomic_positions, atomic_species = get_cell_data(fin)  # get some data on the unit cell
    charge, chargediff = get_charge(fin)    	# get the charge (and charge diff)

    G = compute_G(b, charge.shape)
    fig1 = plot2D_FFTinterp(charge, G, a, x0=(0, 0, 0), e1=(1, 0, 0), e2=(0, 1, 0),
                            nx=100, ny=100, plot_file='plotfile')
    fig1.savefig("figure_1.png")
