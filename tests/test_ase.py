#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
#
from ase import Atoms
from postqe.ase.io import read_espresso_output

if __name__ == "__main__":
    from ase.calculators.calculator import FileIOCalculator, Calculator, kpts2ndarray

    Ni = read_espresso_output('Ni.xml')
    #test= kpts2ndarray({'path': 'GXG', 'npoints': 200},Ni)
    test= kpts2ndarray([2, 2, 2],Ni)
    print (test)
    exit()
    from postqe.ase.calculator import PostqeCalculator
    from ase.visualize import view

    system = 'Ni'
    test = read_espresso_output(system + '.xml')
    test.set_calculator(PostqeCalculator(label = system))
    print (test.get_atomic_numbers())
    print (test.get_cell(True))
    print (test.get_positions())
    print (test.get_volume())
    #view(test)

    #Ni2.calc.read_results()
    print(test.get_potential_energy())
    print(test.calc.get_xc_functional())
    print(test.calc.get_number_of_spins())
    print(test.calc.get_spin_polarized())
    print(test.calc.get_fermi_level())
    print(test.calc.get_number_of_bands())
    print(test.calc.get_bz_k_points())
    print(test.calc.get_k_point_weights())
    print(test.calc.get_eigenvalues(0,0))
    print(test.calc.get_occupation_numbers(0,0))
    print(test.calc.get_eigenvalues(0,1))
    print(test.calc.get_occupation_numbers(0,1))



