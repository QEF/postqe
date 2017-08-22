#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the the band structure of silicon.
"""
    
if __name__ == "__main__":

    from postqe.ase.calculator import Postqe_calc_full
    import numpy as np

    from ase import Atoms
    from ase.visualize import view

    from postqe import compute_bands, plot_bands

    fin = 'Si.xml'
    kpoints, bands = compute_bands(fin, filebands='filebands', spin_component='')

    fig1 = plot_bands(kpoints,bands)





