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
XML data access classes for postqe.
"""
from collections import MutableMapping
import numpy as np
import xmlschema


class XMLData(MutableMapping):
    """
    Dictionary-like class for mapping data from an XML file.
    """
    def __init__(self, xmlfile=None, schema=None):
        if xmlfile is not None:
            self.read(xmlfile, schema)
        else:
            self._data = {}

    def __repr__(self):
        return '<%s %s at %#x>' % (self.__class__.__name__, str(self._data), id(self))

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __getattr__(self, key):
        try:
            return super(XMLData, self).__getattribute__(key)
        except AttributeError:
            if key not in self._data:
                raise
            return self._data[key]

    def __setattr__(self, key, value):
        if key == '_data':
            super(XMLData, self).__setattr__(key, value)
        else:
            self._data[key] = value

    def __delattr__(self, key):
        if key == '_data':
            super(XMLData, self).__delattr__(key)
        else:
            del self._data[key]

    def read(self, xmlfile, schema=None):
        self._data = xmlschema.to_dict(xmlfile, schema)


class PWData(XMLData):

    def read(self, xmlfile, schema=None):
        self._data = {}
        data = xmlschema.to_dict(xmlfile, schema)
        try:
            self['pseudodir'] = data["input"]["control_variables"]["pseudo_dir"]
        except KeyError:
            self['pseudodir'] = './'  # DB: ['./', os.path.dirname(filename)] ... alternatives?

        self["prefix"] = data["input"]["control_variables"]["prefix"]
        self["outdir"] = data["input"]["control_variables"]["outdir"]

        dout = data["output"]
        self["ecutwfc"] = dout["basis_set"]["ecutwfc"]
        self["ecutrho"] = dout["basis_set"]["ecutrho"]
        self["alat"] = dout["atomic_structure"]["@alat"]
        self["ibrav"] = dout["atomic_structure"]["@bravais_index"]

        self['a'] = np.array([
            np.array(dout["atomic_structure"]["cell"]["a1"]),
            np.array(dout["atomic_structure"]["cell"]["a2"]),
            np.array(dout["atomic_structure"]["cell"]["a3"])
        ])
        self['b'] = np.array([
            np.array(dout["basis_set"]["reciprocal_lattice"]["b1"]),
            np.array(dout["basis_set"]["reciprocal_lattice"]["b2"]),
            np.array(dout["basis_set"]["reciprocal_lattice"]["b3"])
        ])
        self["functional"] = np.array(dout["dft"]["functional"])

        self["nat"] = (dout["atomic_structure"]["@nat"])

        # for subsequent loops it is important to have always lists for atomic_positions
        # and atomic_species. If this is not, convert
        a_p = dout["atomic_structure"]["atomic_positions"]["atom"]
        if isinstance(a_p, list):
            self["atomic_positions"] = a_p
        else:
            self["atomic_positions"] = [a_p]

        a_s = dout["atomic_species"]["species"]
        if isinstance(a_s, list):
            self["atomic_species"] = a_s
        else:
            self["atomic_species"] = [a_s]

        self["ntyp"] = dout["atomic_species"]["@ntyp"]
        self["lsda"] = dout["magnetization"]["lsda"]
        self["noncolin"] = dout["magnetization"]["noncolin"]
        self['nr'] = np.array([
            dout["basis_set"]["fft_grid"]["@nr1"],
            dout["basis_set"]["fft_grid"]["@nr2"],
            dout["basis_set"]["fft_grid"]["@nr3"]
        ])
        self['nr_smooth'] = np.array([
            dout["basis_set"]["fft_smooth"]["@nr1"],
            dout["basis_set"]["fft_smooth"]["@nr2"],
            dout["basis_set"]["fft_smooth"]["@nr3"]
        ])
