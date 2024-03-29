#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
# noinspection PyPackageRequirements
from . import pyqe

from .calculator import EspressoCalculator, PostqeCalculator

from .compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from .api import get_eos, get_band_structure, get_dos, get_charge, get_potential, \
    compute_eos, compute_band_structure, compute_dos, compute_charge, compute_potential
from .plot import plot_1d_charge, plot_2d_charge, plot_3d_charge, \
    simple_plot_xy, multiple_plot_xy


__version__ = '1.0.0'
__all__ = [
    'pyqe', 'EspressoCalculator', 'PostqeCalculator', 'compute_v_bare',
    'compute_v_h', 'compute_v_xc', 'get_eos', 'get_band_structure', 'get_dos',
    'get_charge', 'get_potential', 'compute_eos', 'compute_band_structure',
    'compute_dos', 'compute_charge', 'compute_potential', 'plot_1d_charge',
    'plot_2d_charge', 'plot_3d_charge', 'simple_plot_xy', 'multiple_plot_xy'
]
