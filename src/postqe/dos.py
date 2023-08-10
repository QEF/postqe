#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
import numpy as np
from ase.dft.dos import DOS


class QEDOS(DOS):
    """
    A specialization of ASE DOS class with a new get_dos_int method for computing
    the integral of the DOS and a modified write method to write it properly.
    """
    def __init__(self, calc, width=0.1, window=None, npts=201):
        DOS.__init__(self, calc, width=width, window=window, npts=npts)
        self.filename = None
        if window is None:
            self.emin = self.energies.min()
            self.emax = self.energies.max()
        else:
            self.emin, self.emax = window

    def get_dos_int(self):
        """Get array of integral DOS values (always for the total DOS)."""
        dos = self.get_dos()
        dos_int = np.zeros(self.npts)
        temp = 0
        for i in range(len(dos)):
            temp += dos[i]
            dos_int[i] = temp
        return dos_int

    def write(self, filename='dos.out'):
        self.filename = filename

        with open(self.filename, "w") as fout:
            # The header contains some information on the system, the grid nr, etc.
            header = '# DOS file, width={} window=({}, {}) npts={}\n'.format(
                self.width, self.emin, self.emax, self.npts
            )
            fout.write(header)
            fout.write(80 * '#' + '\n')
            dos_int = self.get_dos_int()
            if self.nspins == 2:    # non magnetic
                dos_up = self.get_dos(spin=0)
                dos_down = self.get_dos(spin=1)
                for i in range(len(self.energies)):
                    fout.write('{:.14f}\t{:.14f}\t{:.14f}\t{:.14f}\n'.format(
                        self.energies[i], dos_up[i], dos_down[i], dos_int[i])
                    )
            else:
                dos = self.get_dos()
                for i in range(len(self.energies)):
                    fout.write('{:.14f}\t{:.14f}\t{:.14f}\n'.format(self.energies[i], dos[i], dos_int[i]))
