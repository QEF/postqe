#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
import unittest
import pathlib
import numpy as np

from xmlschema import XMLSchema10

from ase.cell import Cell
from ase.calculators.calculator import kpts2ndarray
from ase.spectrum.band_structure import BandStructure
from postqe import EspressoCalculator, PostqeCalculator, get_band_structure


def abspath(rel_path):
    return str(pathlib.Path(__file__).parent.absolute().joinpath(rel_path))


# Setup directory paths for testing

current_workdir = pathlib.Path('.').resolve()
examples_dir = pathlib.Path(__file__).parent.joinpath('examples')

try:
    rel_examples_dir = examples_dir.relative_to(current_workdir)
except ValueError:
    rel_examples_dir = None


def compare_data(datafile1, datafile2, header=0, tolerance=0.0001):
    """Compare two data files, discarding the first 'header' lines."""
    lr = 1 - tolerance
    hr = 1 + tolerance
    if lr < hr:
        def compare_values(x, y):
            return lr < x / y < hr
    else:
        def compare_values(x, y):
            return x == y

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
                            print("\nValues %s and %s differs "
                                  "(tolerance=%s)." % (v1, v2, tolerance))
                            break
                    return False
    return True


class TestPostqeCalculators(unittest.TestCase):

    @unittest.skipIf(rel_examples_dir is None,
                     'a relative path to examples/ dir is not available!')
    def test_calculator_paths(self):
        si_dir = rel_examples_dir.joinpath('Si')
        label = str(si_dir.joinpath('Si'))

        calculator = EspressoCalculator(label=label)
        self.assertEqual(calculator.prefix, 'Si')
        self.assertEqual(calculator.directory, str(si_dir))
        self.assertEqual(calculator.label, label)

        calculator = EspressoCalculator(label='Si', directory=str(si_dir))
        self.assertEqual(calculator.prefix, 'Si')
        self.assertEqual(calculator.directory, str(si_dir))
        self.assertEqual(calculator.label, label)

        with self.assertRaises(ValueError):
            EspressoCalculator(label='./Si', directory=str(si_dir.absolute()))

        with self.assertRaises(ValueError):
            EspressoCalculator(label=label, directory=str(si_dir))

        directory = str(rel_examples_dir.joinpath('Si'))
        calc = EspressoCalculator(directory=directory)

        self.assertEqual(calc.prefix, 'pwscf')
        self.assertEqual(calc.directory, directory)
        self.assertEqual(calc.label, directory + '/pwscf')

    def test_calculator_schema(self):
        calculator = EspressoCalculator()
        self.assertIsInstance(calculator.schema, XMLSchema10)
        self.assertEqual(calculator.schema.name, 'qes.xsd')
        self.assertTrue(calculator.schema.url.endswith('schemas/qes.xsd'))

        calculator = EspressoCalculator(schema='qes_200420.xsd')
        self.assertIsInstance(calculator.schema, XMLSchema10)
        self.assertEqual(calculator.schema.name, 'qes_200420.xsd')
        self.assertTrue(calculator.schema.url.endswith('schemas/releases/qes_200420.xsd'))

    def test_calculator_command(self):
        calculator = EspressoCalculator(label='examples/Si')
        self.assertIsInstance(calculator.command, str)
        self.assertTrue(calculator.command.endswith('pw.x < PREFIX.in > PREFIX.out'))

        calculator = PostqeCalculator(label='examples/Si')
        self.assertIsNone(calculator.command)

        with self.assertRaises(RuntimeError) as ec:
            calculator.calculate()
        self.assertIn("only for QE results post-processing", str(ec.exception))


class TestNiCase01(unittest.TestCase):

    label = 'Ni_pz_nc'
    calc: PostqeCalculator

    @classmethod
    def setUpClass(cls):
        cls.calc = PostqeCalculator(
            label=cls.label,
            directory=str(examples_dir),
            schema='qes-20180510.xsd'
        )
        cls.calc.read_results()

    def test_get_atoms_from_xml_output(self):
        ni_atoms = self.calc.get_atoms_from_xml_output()
        result = [[-0.25, -0.25, -0.25],
                  [-0.25, -0.25,  0.25],
                  [-0.25,  0.25, -0.25],
                  [-0.25,  0.25,  0.25],
                  [+0.25, -0.25, -0.25],
                  [+0.25, -0.25,  0.25],
                  [+0.25,  0.25, -0.25],
                  [+0.25,  0.25,  0.25]]
        self.assertListEqual(kpts2ndarray([2, 2, 2], ni_atoms).tolist(), result)

    def test_get_atomic_numbers_from_atoms(self):
        ni_atoms = self.calc.get_atoms_from_xml_output()
        atomic_numbers = ni_atoms.get_atomic_numbers()
        self.assertIsInstance(atomic_numbers, np.ndarray)
        self.assertEqual(atomic_numbers, [28])

    def test_get_cell_from_atoms(self):
        ni_atoms = self.calc.get_atoms_from_xml_output()
        cell = ni_atoms.get_cell(complete=True)
        self.assertIsInstance(cell, Cell)
        self.assertIsInstance(cell.array, np.ndarray)
        self.assertListEqual(
            cell.array.tolist(),
            [[-3.325, 0.0, 3.325], [0.0, 3.325, 3.325], [-3.325, 3.325, 0.0]]
        )

    def test_get_positions_from_atoms(self):
        ni_atoms = self.calc.get_atoms_from_xml_output()
        positions = ni_atoms.get_positions()
        self.assertIsInstance(positions, np.ndarray)
        self.assertListEqual(positions.tolist(), [[0.0, 0.0, 0.0]])

    def test_get_volume_from_atoms(self):
        ni_atoms = self.calc.get_atoms_from_xml_output()
        volume = ni_atoms.get_volume()
        self.assertIsInstance(volume, float)
        self.assertEqual(volume, 73.51990624999998)

    def test_get_potential_energy(self):
        potential_energy = self.calc.get_potential_energy()
        self.assertEqual(potential_energy, -414.23324992553125)

    def test_get_xc_functional(self):
        self.assertEqual(self.calc.get_xc_functional(), 'PZ')

    def test_get_number_of_spins(self):
        self.assertEqual(self.calc.get_number_of_spins(), 2)

    def test_get_spin_polarized(self):
        self.assertTrue(self.calc.get_spin_polarized())

    def test_get_fermi_level(self):
        self.assertEqual(self.calc.get_fermi_level(), 20.105067177719004)

    def test_get_number_of_bands(self):
        self.assertEqual(self.calc.get_number_of_bands(), 18)

    def test_get_k_point_weights(self):
        weights = self.calc.get_k_point_weights()
        self.assertIsInstance(weights, np.ndarray)
        self.assertListEqual(
            weights.tolist(),
            [0.07407407407407407, 0.2222222222222222,
             0.2222222222222222, 0.2222222222222222,
             0.2222222222222222, 0.03703703703703703]
        )

    def test_get_eigenvalues(self):
        eigenvalues = self.calc.get_eigenvalues(0, 0)
        self.assertIsInstance(eigenvalues, np.ndarray)
        self.assertListEqual(
            eigenvalues.tolist(),
            [6.0931957834356165, 16.22730531613425, 16.54646013921364,
             16.54646013921364, 18.543772203139113, 18.543772203139113,
             33.95164971822426, 38.847247548667376, 38.8472475767997]
        )

        eigenvalues = self.calc.get_eigenvalues(0, 1)
        self.assertIsInstance(eigenvalues, np.ndarray)
        self.assertListEqual(
            eigenvalues.tolist(),
            [6.093195875780507, 16.22730680526378, 16.546461820745535,
             16.546461820745535, 18.543773841967827, 18.543773841967827,
             33.95165023183349, 38.8472478753786, 38.84724789445674]
        )

    def test_get_occupation_numbers(self):
        occupation_numbers = self.calc.get_occupation_numbers(0, 0)
        self.assertIsInstance(occupation_numbers, np.ndarray)
        self.assertListEqual(
            occupation_numbers.tolist(),
            [1.0, 1.0, 1.0, 1.0, 1.000000000000008, 1.000000000000008, -1.986509392945334e-86,
             -2.688859669326466e-86, -2.68885967336249e-86]
        )

        occupation_numbers = self.calc.get_occupation_numbers(0, 2)
        self.assertIsInstance(occupation_numbers, np.ndarray)
        self.assertListEqual(
            occupation_numbers.tolist(),
            [1.0, 1.0, 1.0, 1.0, 1.000000000000008, 1.000000000000008, -1.98650946663063e-86,
             -2.688859716198313e-86, -2.688859718935374e-86]
        )

    def test_get_bz_k_points(self):
        result = self.calc.get_bz_k_points()

    def test_get_ibz_k_points(self):
        with self.assertRaises(NotImplementedError):
            self.calc.get_ibz_k_points()

    def test_get_pseudo_density(self):
        result = self.calc.get_pseudo_density()

    def test_get_effective_potential(self):
        with self.assertRaises(NotImplementedError):
            self.calc.get_effective_potential()

    def get_pseudo_wave_function(self):
        result = self.calc.get_pseudo_wave_function()


class TestNiCase02(TestNiCase01):

    label = 'Ni_pbe_us'

    def test_get_potential_energy(self):
        potential_energy = self.calc.get_potential_energy()
        self.assertEqual(potential_energy, -584.3311788465203)

    def test_get_xc_functional(self):
        self.assertEqual(self.calc.get_xc_functional(), 'PBE')

    def test_get_fermi_level(self):
        self.assertEqual(self.calc.get_fermi_level(), 14.119346198988783 )

    def test_get_eigenvalues(self):
        eigenvalues = self.calc.get_eigenvalues(0, 0)
        self.assertIsInstance(eigenvalues, np.ndarray)
        self.assertListEqual(
            eigenvalues.tolist(),
            [6.348975841094524, 11.467210784393725, 11.864196654943136, 11.86419665494325,
             12.820196483280753, 12.820196483281327, 33.07284580473358, 38.87662604267879,
             40.65476430191017]
        )

        eigenvalues = self.calc.get_eigenvalues(0, 1)
        self.assertIsInstance(eigenvalues, np.ndarray)
        self.assertListEqual(
            eigenvalues.tolist(),
            [6.244436910933738, 12.271436044529775, 12.710188610698562, 12.710188610698733,
             13.530331317147054, 13.530331317147443, 33.08582596886793, 38.66354541182688,
             40.65219032336821]
        )

    def test_get_occupation_numbers(self):
        occupation_numbers = self.calc.get_occupation_numbers(0, 0)
        self.assertIsInstance(occupation_numbers, np.ndarray)
        self.assertListEqual(
            occupation_numbers.tolist(),
            [1.0, 1.0, 1.0, 1.0, 1.000000000162555, 1.000000000162555, -2.719176727243473e-86,
             -3.551820011139861e-86, -3.806921827316206e-86]
        )

        occupation_numbers = self.calc.get_occupation_numbers(0, 2)
        self.assertIsInstance(occupation_numbers, np.ndarray)
        self.assertListEqual(
            occupation_numbers.tolist(),
            [1.0, 1.0, 1.000000000003175, 1.000000000003175, 1.004532869927512, 1.004532869927536,
             -2.721038935329794e-86, -3.521250253338337e-86, -3.806552549731944e-86]
        )


class TestPostqeAPI(unittest.TestCase):

    def test_get_band_structure(self):
        bs = get_band_structure(
            prefix=abspath('examples/Si'),
            schema="releases/qes-20180510.xsd",
            reference_energy=0
        )
        self.assertIsInstance(bs, BandStructure)


if __name__ == '__main__':
    unittest.main()
