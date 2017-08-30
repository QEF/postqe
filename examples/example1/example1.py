#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":
    from postqe import units, get_eos

    eos = get_eos(label="./Nienergies.dat", eos='murnaghan')
    v0, e0, B = eos.fit()
    # Print some data and plot
    print('Equilibrium volume = '+str(v0)+' Ang^3')
    print('Equilibrium energy = '+str(e0)+' eV')
    print('Equilibrium Bulk modulus = '+str(B / units.kJ * 1.0e24)+' GPa')
    fig = eos.plot('Ni-eos.png', show=True)

    # Save the plot in a different format (pdf) with Matplotlib if you like
    fig.figure.savefig("figure.pdf", format='pdf')
