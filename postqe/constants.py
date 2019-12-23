#
# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
"""
Some useful constants
"""

from numpy import pi

tpi = 2.0 * pi
fpi = 4.0 * pi

C_SI = 2.99792458E+8                # m sec^-1
H_PLANCK_SI = 6.62606896E-34        # J s
K_BOLTZMANN_SI = 1.3806504E-23      # J K^-1
HARTREE_SI = 4.35974394E-18         # J
BOHR_RADIUS_SI = 0.52917720859E-10  # m
RYDBERG_SI = HARTREE_SI/2.0         # J
AU_SEC = H_PLANCK_SI/tpi/HARTREE_SI
AU_PS = AU_SEC * 1.0E+12
AU_TERAHERTZ = AU_PS
AU_GPA = HARTREE_SI / BOHR_RADIUS_SI ** 3 / 1.0E+9
K_BOLTZMANN_RY = K_BOLTZMANN_SI / RYDBERG_SI
RY_TO_THZ = 1.0 / AU_TERAHERTZ / fpi
RY_TO_GHZ = RY_TO_THZ * 1000.0
RY_TO_CMM1 = 1.0E+10 * RY_TO_THZ / C_SI
RY_KBAR = 10.0 * AU_GPA / 2.0

kb1 = 1.0 / K_BOLTZMANN_RY / RY_TO_CMM1  # inverse Boltzmann constant in cm^{-1}/K
ev_to_ry = 0.073498618
