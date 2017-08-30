#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
Tests for postqe.
"""
import unittest
import sys
import os

# Adds the the package directory to sys.path, in order to make
# the development module loadable also without set PYTHONPATH.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
if sys.path[0] != PACKAGE_DIR:
    sys.path.insert(0, PACKAGE_DIR)


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

    def test_get_band_structure(self):
        from postqe import get_band_structure
        import pdb
        pdb.set_trace()
        bs = get_band_structure(label="./examples/Si", schema='../schemas/qes.xsd', reference_energy=0)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
