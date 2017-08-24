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


    def read_results(self):
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
        dout = d["output"]
        self.results['energy'] = float(dout["total_energy"]["etot"]) * units.Ry


