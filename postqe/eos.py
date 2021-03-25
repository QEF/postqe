#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
from ase.eos import EquationOfState

HEADER_TMPL = ("# EOS fitting, equation: {}\n"
               "# Equilibrium volume: {}\n"
               "# Equilibrium energy: {}\n"
               "# Bulk modulus: {}\n")


class QEEquationOfState(EquationOfState):
    """
    A specialization of ASE EquationOfState class with a modified default
    EOS type and a write method.
    """
    def __init__(self, volumes, energies, eos='murnaghan'):
        EquationOfState.__init__(self, volumes, energies, eos=eos)
        self.filename = None

    def write(self, filename='eos.out'):
        self.filename = filename

        with open(self.filename, "w") as fout:
            # The header contains some information on the system, the grid nr, etc.
            fout.write(HEADER_TMPL.format(self.eos_string, self.v0, self.e0, self.B))
            fout.write(80 * '#' + '\n')

            for i in range(len(self.v)):
                fout.write('{}\t{}'.format(self.v[i], self.e[i]))
                if self.eos_string == 'sj':
                    y = self.fit0(self.v[i] ** -(1 / 3))
                else:
                    y = self.func(self.v[i], *self.eos_parameters)
                fout.write('\t{}\n'.format(y))
