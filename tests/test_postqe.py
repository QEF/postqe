#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2018, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
import unittest
import os
import numpy as np

from ase.calculators.calculator import kpts2ndarray
from postqe import PostqeCalculator


def compare_data(datafile1, datafile2, header=0, tolerance=0.0001):
    """Compare two data files, discarding the first 'header' lines."""
    lr = 1 - tolerance
    hr = 1 + tolerance
    if lr < hr:
        compare_values = lambda x, y: lr < x / y < hr
    else:
        compare_values = lambda x, y: x == y

    with open(datafile1) as f1, open(datafile2) as f2:
        for line1, line2 in zip(f1, f2):
            if header > 0:
                header -= 1
            else:
                values1 = [float(v) for v in line1.split()]
                values2 = [float(v) for v in line2.split()]
                if len(values1) != len(values2):
                    return False
                if not all([compare_values(v1, v2) for v1, v2 in zip(values1, values2)]):
                    for v1, v2 in zip(values1, values2):
                        if not compare_values(v1, v2):
                            print("\nValues %s and %s differs (tolerance=%s)." % (v1, v2, tolerance))
                            break
                    return False
    return True


class TestPostQE(unittest.TestCase):

    def abspath(self, rel_path):
        return os.path.join(os.path.abspath(os.path.dirname(__file__)), rel_path)

    def test_get_band_structure(self):
        from ase.dft.band_structure import BandStructure
        from postqe import get_band_structure

        bs = get_band_structure(
            prefix=self.abspath('examples/Si'),
            schema="releases/qes-20180510.xsd",
            reference_energy=0
        )
        self.assertIsInstance(bs, BandStructure)

    def test_get_atoms_from_xml_output(self):
        system = 'Ni'
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples/Ni_pz_nc')
        calc = PostqeCalculator(label=system, outdir=outdir, schema='qes-20180510.xsd')
        calc.read_results()

        ni_atoms = calc.get_atoms_from_xml_output()
        result = [[-0.25, -0.25, -0.25],
                  [-0.25, -0.25,  0.25],
                  [-0.25,  0.25, -0.25],
                  [-0.25,  0.25,  0.25],
                  [+0.25, -0.25, -0.25],
                  [+0.25, -0.25,  0.25],
                  [+0.25,  0.25, -0.25],
                  [+0.25,  0.25,  0.25]]
        self.assertListEqual(kpts2ndarray([2, 2, 2], ni_atoms).tolist(), result)

    def test_get_atomic_numbers(self):
        system = 'Ni'
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'examples/Ni_pbe_us')
        calc = PostqeCalculator(label=system, outdir=outdir, schema='qes-20180510.xsd')
        calc.read_results()

        ni_atoms = calc.get_atoms_from_xml_output()
        print(ni_atoms.get_atomic_numbers())

        print(ni_atoms.get_cell(True))
        print(ni_atoms.get_positions())
        print(ni_atoms.get_volume())

        print(calc.get_potential_energy())
        print(calc.get_xc_functional())
        print(calc.get_number_of_spins())
        print(calc.get_spin_polarized())
        print(calc.get_fermi_level())
        print(calc.get_number_of_bands())
        print(calc.get_bz_k_points())
        print(calc.get_k_point_weights())
        print(calc.get_eigenvalues(0, 0))
        print(calc.get_occupation_numbers(0, 0))
        print(calc.get_eigenvalues(0, 1))
        print(calc.get_occupation_numbers(0, 1))


if __name__ == '__main__':
    unittest.main()
