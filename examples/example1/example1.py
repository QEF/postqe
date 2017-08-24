#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":
    import os.path
    from postqe import read_EtotV

    # file with the total energy data E(V)
    fin = os.path.join(os.path.dirname(__file__), "EtotV.dat")

    from ase.units import kJ
    from ase.eos import EquationOfState

    # Extract volumes and energies from the input file:
    volumes, energies = read_EtotV(fin)

    # Create an object EquationOfState and fit with Murnaghan (or other) EOS
    eos = EquationOfState(volumes, energies, eos='murnaghan')
    v0, e0, B = eos.fit()
    # Print some data and plot
    print('Equilibrium volume = '+str(v0)+' Ang^3')
    print('Equilibrium energy = '+str(e0)+' eV')
    print('Equilibrium Bulk modulus = '+str(B / kJ * 1.0e24)+' GPa')
    eos.plot('Ni-eos.png')

