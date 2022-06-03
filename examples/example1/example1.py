#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS using the "compute_eos" API.

The file with volumes and enegies can be plotted using ASE and the postqe EspressoCalculator

OPTIONS compute_eos:
    prefix: name of the input file with volumes and energies
    outdir: directory containing the input data
    eos_type: type of (EOS), murnaghan, sjeos, taylor, vinet, birch, birchmurnaghan, pouriertarantola, antonschmidt
    fileout: output file with fitting data and results, if not specified the file will not be written
    fileplot: output plot file (default='EOSplot') in png format
    show: True -> plot results with Matplotlib; None or False -> do nothing. Default = True
    ax: a Matplotlib "Axes" instance (see Matplotlib documentation for details)

Returns: an QEEquationOfState object (eos) and a Matplotlib figure object (eos_plot)
"""

if __name__ == "__main__":
    import numpy as np
    from ase import Atoms
    from postqe import compute_eos, EspressoCalculator

    #Runnin 15 scf calculations with ASE and EspressoCalculator to generate a file with Volumes and Energies
    a = 3.2 # approximated lattice constant
    ni = Atoms('Ni', #Nickel bulk
              cell = [(0, a/2, a/2), (a/2, 0, a/2), (a/2, a/2, 0)], #FCC cell
              pbc = 1, #activate periodic boudary conditions
              calculator = EspressoCalculator(calculation='scf', ecutwfc=100, occupations='smearing', smearing='gaussian',
                                                degauss=0.001, kpoints=[6,6,6,0,0,0], conv_thr=1e-8, pseudo_dir='../'))
    cell = ni.get_cell()
    
    f = open('vol_and_energies.dat', 'w')
    for cell_multi in np.linspace(1.85, 2.15, 15):
        ni.set_cell(cell*cell_multi, scale_atoms=True)
        v = ni.get_volume()
        e = ni.get_potential_energy()
        f.write(('\t{:.14f}\t{:.14f}\n').format(v,e))
    f.close()

    #call the compute_dos API
    eos, eos_plot = compute_eos(prefix='vol_and_energies.dat', outdir='.',eos_type='murnaghan',
                                        fileout='eos.out', fileplot='EOSplot', show=True, ax=None)

    ## Save the plot in a different formats (pdf) with Matplotlib if you like

    # fig = eos.plot()
    # fig.figure.savefig('figure.jpg', format='jpg')
    # fig.figure.savefig('figure.pdf', format='pdf')
