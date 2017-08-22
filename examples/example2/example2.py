#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the the band structure of silicon.
"""
    
if __name__ == "__main__":

    from ase import Atoms
    import numpy as np
    from postqe.ase.io import read_espresso_xml
    from postqe.ase.calculator import Postqe_calc_full

    # Set the calculator and parameters for Quantum Espresso
    QEparameters = {'outdir': 'temp', 'smearing': 'mp', 'occupations': 'smearing', 'degauss': 0.02,
                    'pp_dict': { 'Ni': 'Ni.pz-n-rrkjus_psl.1.0.0.UPF', 'Ag': 'Ag.pz-n-rrkjus_psl.1.0.0.UPF'},
                    }

    calcul = Postqe_calc_full(atoms=None, label='./Ni', restart=None, ibrav=0, ecutwfc=50,
                              kpts={'path': 'GXWLGK', 'npoints': 200}, tstress=True, tprnfor=True,
                              command='/home/mauropalumbo/q-e/bin/pw.x < PREFIX.in > PREFIX.out',
                              pseudo_dir='../PSEUDOPOTENTIALS', **QEparameters)

    Ni = read_espresso_xml('Ni.xml')
    Ni.set_calculator(calcul)

    Ni.get_potential_energy()
    bs = Ni.calc.band_structure()
    print(bs.energies)
    bs.plot(emax=10, filename='ni.png')





