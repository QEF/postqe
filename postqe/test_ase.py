#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################

from ase import Atoms
from postqe.aseio import read_espresso_xml

if __name__ == "__main__":
    from postqe.asecalc import Espresso
    from ase.visualize import view

    Ni2 = read_espresso_xml('Ni.xml')
    #Ni2 = Atoms()
    Ni2.set_calculator(Espresso(label = 'Ni'))
    print (Ni2.get_atomic_numbers())
    print (Ni2.get_cell(True))
    print (Ni2.get_positions())
    print (Ni2.get_volume())
    view(Ni2)


    #Ni2.calc.read_results()
    print(Ni2.get_potential_energy())
    print(Ni2.calc.get_xc_functional())
    print(Ni2.calc.get_fermi_level())


