#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of plotting a 1D section of the electronic charge density.
"""
    
if __name__ == "__main__":

    import numpy as np
    from postqe import get_charge, get_cell_data, get_calculation_data, \
        compute_v_bare, compute_v_h, compute_v_xc, compute_G, plot1D_FFTinterp

    fin = "./Ni.xml"  				            # file xml produce by QE
    ibrav, alat, a, b, nat, ntyp,\
    atomic_positions, atomic_species = get_cell_data(fin)  # get some data on the unit cell
    prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth =\
        get_calculation_data(fin)               # get some data on the QE calculation

    charge, chargediff = get_charge(fin)    	# get the charge (and charge diff) from the HDF5 file

    # Compute the bare potential (v_bare), the Hartree potential (v_h) and the exhange-correlation potential (v_xc)
    # Add them up to get the total potential
    v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir)
    v_h = compute_v_h(charge, ecutrho, alat, b)
    charge_core = np.zeros(charge.shape)
    v_xc = compute_v_xc(charge, charge_core, str(functional))
    v_tot = v_bare + v_h + v_xc

    G = compute_G(b, charge.shape)      # calculate the G vectors for plotting
    # Plot the potentials as 1D sections from (0,0,0) along (1,0,0)
    fig1 = plot1D_FFTinterp(v_bare, G, a, x0=(0, 0, 0), e1=(1, 0, 0), nx=50, ylab='v_bare', plot_file='plot_v_bare')
    fig1.savefig("figure_v_bare.png")

    fig2 = plot1D_FFTinterp(v_h, G, a, x0=(0, 0, 0), e1=(1, 0, 0), nx=50, ylab='v_h', plot_file='plot_v_h')
    fig2.savefig("figure_v_h.png")

    fig3 = plot1D_FFTinterp(v_xc, G, a, x0=(0, 0, 0), e1=(1, 0, 0), nx=50, ylab='v_xc', plot_file='plot_v_xc')
    fig3.savefig("figure_v_xc.png")

    fig4 = plot1D_FFTinterp(v_tot, G, a, x0=(0, 0, 0), e1=(1, 0, 0), nx=50, ylab='v_tot', plot_file='plot_v_tot')
    fig4.savefig("figure_v_tot.png")
