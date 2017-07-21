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

# Adds the parent oth the test directory to sys.path,
# in order to make the development module loadable.
os.chdir(os.path.dirname(__file__) or '.')
pkg_search_path = os.path.abspath('..')
if sys.path[0] != pkg_search_path:
    sys.path.insert(0, pkg_search_path)


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
        test_dir = "./Ni_pz_nc"

    def test_003(self):
        """
        Test case: magnetic Ni fcc with ultrasoft PBE pseudo
        """
        test_dir = "./Ni_pbe_us"


if __name__ == '__main__':
    import postqe
    print("module = ", postqe)
    unittest.main()
