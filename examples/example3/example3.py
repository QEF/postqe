#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the electronic density of states (DOS or edos).

First it computes the DOS for silicon (input file Si.xml), then for magnetic Ni (lsda) showing how to
obtain the dos for spin up and down and how to plot it.
"""
    
if __name__ == "__main__":
    from ase.dft import DOS
    from postqe.ase.io import read_espresso_xml
    from postqe.ase.calculator import Postqe_calc_full

    # Set the calculator and parameters for Quantum Espresso
    QEparameters = {'outdir': 'temp', 'smearing': 'mp', 'occupations': 'smearing', 'degauss': 0.02,
                    'pp_dict': { 'Ni': 'Ni.pz-n-rrkjus_psl.1.0.0.UPF', 'Ag': 'Ag.pz-n-rrkjus_psl.1.0.0.UPF'},
                    }

    calcul = Postqe_calc_full(atoms=None, label='./Ni', restart=None, ibrav=0, ecutwfc=50,
                              kpts=[8, 8, 8], tstress=True, tprnfor=True,
                              command='/home/mauropalumbo/q-e/bin/pw.x < PREFIX.in > PREFIX.out',
                              pseudo_dir='../PSEUDOPOTENTIALS', **QEparameters)

    Ni = read_espresso_xml('Ni.xml')
    Ni.set_calculator(calcul)
    Ni.get_potential_energy()

    dos = DOS(calcul, width=0.02)
    d = dos.get_dos()
    e = dos.get_energies()

    import matplotlib.pyplot as plt

    plt.plot(e, d)
    plt.xlabel('energy [eV]')
    plt.ylabel('DOS')
    plt.show()


