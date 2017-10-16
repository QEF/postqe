#!/usr/bin/env python3
#encoding: UTF-8

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

    def __init__(self, volumes, energies, eos_type='murnaghan'):
        EquationOfState.__init__(self, volumes, energies, eos=eos_type)

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
