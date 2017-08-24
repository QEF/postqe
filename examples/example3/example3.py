#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the electronic density of states (DOS or edos).

First it computes the DOS for silicon (input file Si.xml), then for magnetic Ni (lsda) showing how to
obtain the dos for spin up and down and how to plot it.
"""
    
if __name__ == "__main__":
    from ase.dft import DOS
    from postqe.ase.io import read_espresso_output
    from postqe.ase.calculator import EspressoCalculator

    # Set the calculator and parameters for Quantum Espresso
    QEparameters = {'outdir': 'temp', 'smearing': 'mp', 'occupations': 'smearing', 'degauss': 0.02,
                    'pp_dict': { 'Ni': 'Ni.pz-n-rrkjus_psl.1.0.0.UPF', 'Ag': 'Ag.pz-n-rrkjus_psl.1.0.0.UPF'},
                    }

    calcul = EspressoCalculator(atoms=None, label='./Ni',
                                schema='../../schemas/qes.xsd',
                                restart=None, ibrav=0, ecutwfc=50,
                                kpts=[8, 8, 8], tstress=True, tprnfor=True,
                                command='../../postqe/fortran/build/q-e/bin/pw.x < PREFIX.in > PREFIX.out',
                                pseudo_dir='../PSEUDOPOTENTIALS', **QEparameters)

    Ni = read_espresso_output('Ni.xml', schema='../../schemas/qes.xsd')
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


