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

import numpy as np
from postqe.readutils import read_charge_file_hdf5, write_charge, create_header
from postqe.compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from postqe.pyqe import pyqe_getcelldms
from postqe.xmldata import PWData


def run_test_case(pseudo_dir):
    # get some needed values from the XML file
    data = PWData(os.path.join(pseudo_dir, "Ni.xml"), os.path.join(PACKAGE_DIR, 'postqe/schemas/qes.xsd'))
    celldms = pyqe_getcelldms(data.alat, data.a[0], data.a[1], data.a[2], data.ibrav)
    charge_file = os.path.join(pseudo_dir, "charge-density.hdf5")

    charge, chargediff = read_charge_file_hdf5(charge_file, data.nr)
    header = create_header("Ni", data.nr, data.nr_smooth, data.ibrav, celldms, data.nat,
                           data.ntyp, data.atomic_species, data.atomic_positions)

    # plot_num = 0
    write_charge(os.path.join(pseudo_dir, 'postqeout0'), charge, header)

    # plot_num = 6
    write_charge(os.path.join(pseudo_dir, 'postqeout6'), chargediff, header)

    # plot_num = 2
    v_bare = compute_v_bare(
        data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
        data.atomic_positions, data.atomic_species, pseudo_dir
    )
    write_charge(os.path.join(pseudo_dir, 'postqeout2'), v_bare, header)

    # plot_num = 11
    v_bare = compute_v_bare(data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
                            data.atomic_positions, data.atomic_species, pseudo_dir)
    v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
    v_tot = v_bare + v_h
    write_charge(os.path.join(pseudo_dir, 'postqeout11'), v_tot, header)

    # plot_num = 1
    v_bare = compute_v_bare(
        data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
        data.atomic_positions, data.atomic_species, pseudo_dir
    )
    v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
    charge_core = np.zeros(charge.shape)
    v_xc = compute_v_xc(charge, charge_core, str(data.functional))
    v_tot = v_bare + v_h + v_xc
    write_charge(os.path.join(pseudo_dir, 'postqeout1'), v_tot, header)


def compare_data(datafile1, datafile2, header=0):
    """Compare two data files, discarding the first 'header' lines."""
    with open(datafile1) as f1, open(datafile2) as f2:
        for line1, line2 in zip(f1, f2):
            if header > 0:
                header -= 1
            else:
                values1 = line1.split()
                values2 = line2.split()
                if len(values1) != len(values2):
                    return False
                if any([float(v1) != float(v2) for v1, v2 in zip(values1, values2)]):
                    print(values1)
                    print(values2)
                    return False
    return True

class TestPostQE(unittest.TestCase):

    def test_001(self):
        """
        Test case: diamond-like Si LDA norm-conserving
        """
        test_dir = "./Si"

    def test_002(self):
        """
        Test case: magnetic Ni fcc with norm-conserving LDA pseudo
        """
        run_test_case(os.path.join(TEST_DIR, "Ni_pz_nc"))
        for k in (0, 1, 2, 6, 11):
            self.assertTrue(compare_data(
                os.path.join(TEST_DIR, 'Ni_pz_nc/reference/ppout%d' % k),
                os.path.join(TEST_DIR, 'Ni_pz_nc/postqeout%d' % k),
                header=6),
                msg="Data file postqeout%d differs!" % k
            )

    def test_003(self):
        """
        Test case: magnetic Ni fcc with ultrasoft PBE pseudo
        """
        test_dir = "./Ni_pbe_us"


if __name__ == '__main__':
    import postqe
    print("module = ", postqe)
    unittest.main()
