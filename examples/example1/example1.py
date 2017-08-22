#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the volume
with the Murnaghan EOS.
"""
    
if __name__ == "__main__":
    from postqe.ase.calculator import Postqe_calc_full
    import numpy as np

    from ase import Atoms
    from ase.io import read
    from ase.units import kJ
    from ase.eos import EquationOfState
    from ase.io.trajectory import Trajectory
    from ase.visualize import view

    # Set the calculator and parameters for Quantum Espresso
    QEparameters = {'outdir': 'temp', 'smearing': 'mp', 'occupations': 'smearing', 'degauss': 0.02,
                    'pp_dict': { 'Ni': 'Ni.pz-n-rrkjus_psl.1.0.0.UPF', 'Ag': 'Ag.pz-n-rrkjus_psl.1.0.0.UPF'},
                    }

    calcul = Postqe_calc_full(atoms=None, label='./Ni', restart=None, ibrav=0, ecutwfc=50,
                              kpoints=[3, 3, 3, 0, 0, 0], tstress=True, tprnfor=True,
                              command='/home/mauropalumbo/q-e/bin/pw.x < PREFIX.in > PREFIX.out',
                              pseudo_dir='../PSEUDOPOTENTIALS', **QEparameters)

    a = 6.5  # approximate lattice constant
    b = a / 2
    Ni = Atoms('Ni',
               cell=[(0, b, b), (b, 0, b), (b, b, 0)],
               pbc=1,
               calculator=calcul)
    cell = Ni.get_cell()
    traj = Trajectory('Ni.traj', 'w')
    for x in np.linspace(0.95, 1.05, 9):
        Ni.set_cell(cell * x, scale_atoms=True)
        Ni.get_potential_energy()
        traj.write(Ni)

    configs = read('Ni.traj@0:9')  # read 9 configurations
    # Extract volumes and energies:
    volumes = [Ni.get_volume() for Ni in configs]
    energies = [Ni.get_potential_energy() for Ni in configs]
    eos = EquationOfState(volumes, energies, eos='murnaghan')
    v0, e0, B = eos.fit()
    print('Equilibrium volume = '+str(v0)+' Ang^3')
    print('Equilibrium energy = '+str(e0)+' eV')
    print('Equilibrium Bulk modulus = '+str(B / kJ * 1.0e24)+' GPa')
    eos.plot('Ni-eos.png')

