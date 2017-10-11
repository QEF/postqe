#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
from ase.dft.dos import DOS

class QEDOS(DOS):

    def __init__(self, calc, width=0.1, window=None, npts=201):
        DOS.__init__(self, calc, width=width, window=window, npts=npts)
        if window is None:
            self.emin = self.energies.min()
            self.emax = self.energies.max()
        else:
            self.emin, self.emax = window

    def get_dos_int(self):
        """
        Get array of integral DOS values (always for the total DOS).

        """
        dos = self.get_dos()
        dos_int = np.zeros(self.npts)
        temp = 0
        for i in range(0,len(dos)):
            temp += dos[i]
            dos_int[i] = temp
        return dos_int

    def write(self,filename='dos.out'):


        self.filename = filename
        fout = open(self.filename, "w")

        # The header contains some information on the system, the grid nr, etc.
        header = '# DOS file, width= '+str(self.width)+' window= '+str(self.emin)+', '+str(self.emax)+\
                 ' npts='+str(self.npts)+'\n'
        fout.write(header)
        fout.write(80*'#'+'\n')
        dos_int = self.get_dos_int()
        if self.nspins == 2:    # non magnetic
            dos_up = self.get_dos(spin=0)
            dos_down = self.get_dos(spin=1)
            for i in range(0,len(self.energies)):
                fout.write(str(self.energies[i])+'\t'+str(dos_up[i])+'\t'+str(dos_down[i])+'\t'+str(dos_int[i])+'\n')
        else:
            dos = self.get_dos()
            for i in range(0,len(self.energies)):
                fout.write(str(self.energies[i])+'\t'+str(dos[i])+'\t'+str(dos_int[i])+'\n')

        fout.close()
