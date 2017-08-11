#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from ase.calculators.calculator import FileIOCalculator, Calculator, ReadError
import ase.units as units

# Fix python3 types
try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)


class Espresso(Calculator):
    """
    This is a limited implementation of an ASE calculator for postqe.

    It does NOT generate an input file for espresso, NOR does it run pw.
    It reads some properties from an espresso xml output file (written
    according to espresso xml schema). The implemented properties are
    those needed for postprocessing.
    """

    implemented_properties = ['energy', 'forces']
    command = None

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, command=None, **kwargs):
        """File-IO calculator.

        command: str
            Command used to start calculation.
        """

        self.species = None
        #self.pp_dict = pp_dict

        Calculator.__init__(self, restart, ignore_bad_restart_file, label,
                            atoms, **kwargs)

        # Read here the xml file and store a dictionary for the output part (dout)
        import xmlschema
        filename = self.label + '.xml'

        ##########################################################
        # TODO for whatever reason this is not working now
        # schemaLoc = xmlschema.fetch_schema(filename)
        # xs = xmlschema.XMLSchema(schemaLoc)
        #
        # temporary local solution
        xs = xmlschema.XMLSchema('schemas/qes.xsd')
        ##########################################################

        print("Reading xml file: ", filename)
        d = xs.to_dict(filename)
        self.dout = d["output"]


    def calculate(self, atoms=None, properties=['energy'], system_changes=[]):
        """
        This calculate method only reads the results from an xml output file.

        :param atoms:
        :param properties:
        :param system_changes:
        :return:
        """
        Calculator.calculate(self, atoms, properties, system_changes)
        self.read_results()


    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def get_number_of_bands(self):
        """Return the number of bands."""
        return int(self.dout["band_structure"]["nbnd"])

    def get_xc_functional(self):
        """Return the XC-functional identifier.

        'LDA', 'PBE', ..."""
        return self.dout["dft"]["functional"]

    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return np.zeros((1, 3))

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return 1

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return False

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return np.zeros((1, 3))

    def get_k_point_weights(self):
        """Weights of the k-points.

        The sum of all weights is one."""
        return np.ones(1)

    def get_pseudo_density(self, spin=None, pad=True):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        return np.zeros((40, 40, 40))

    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        return np.zeros((40, 40, 40))

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        return np.zeros((40, 40, 40))

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.arange(42, float)

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return np.ones(42)


    def get_fermi_level(self):
        """Return the Fermi level."""
        return float(self.dout["band_structure"]["fermi_energy"]) * units.Ry

    def read_results(self):
        self.results['energy'] = float(self.dout["total_energy"]["etot"]) * units.Ry


