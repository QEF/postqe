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
A specialization of ASE EquationOfState class with a modified default EOS type and a write method.
"""

from ase.eos import EquationOfState

def create_header(eos, v0, e0, B):
    header = '# EOS fitting, equation: '+eos+'\n'
    header += '# Equilibrium volume: '+str(v0)+'\n'
    header += '# Equilibrium energy: '+str(e0)+'\n'
    header += '# Bulk modulus: '+str(B)+'\n'
    return header

class QEEquationOfState(EquationOfState):

    def __init__(self, volumes, energies, eos='murnaghan'):
        EquationOfState.__init__(self, volumes, energies, eos=eos)

    def write(self,filename='eos.out'):
        self.filename = filename
        fout = open(self.filename, "w")

        # The header contains some information on the system, the grid nr, etc.
        header = create_header(self.eos_string, self.v0, self.e0, self.B)
        fout.write(header)
        fout.write(80*'#'+'\n')
        for i in range(0,len(self.v)):
            fout.write(str(self.v[i])+'\t'+str(self.e[i]))
            if self.eos_string == 'sj':
                y = self.fit0(self.v[i] ** -(1 / 3))
            else:
                y = self.func(self.v[i], *self.eos_parameters)
            fout.write('\t'+str(y)+'\n')

        fout.close()
