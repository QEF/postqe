#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the electronic density of states (DOS or edos).

Computes the DOS for silicon (file Si.xml)
"""
    
if __name__ == "__main__":
    from ase.dft import DOS
    from postqe.ase.io import get_atoms_from_xml_output
    from postqe.ase.calculator import PostqeCalculator

    calcul = PostqeCalculator(atoms=None, label='./Si', schema='../../schemas/qes.xsd')

    Si = get_atoms_from_xml_output('Si.xml', schema='../../schemas/qes.xsd')
    Si.set_calculator(calcul)
    Si.calc.read_results()

    dos = DOS(calcul, width=0.02*13.605698066*2, npts=3000)
    d = dos.get_dos()
    e = dos.get_energies()

    # if you want, you can plot the DOS with Matplotlib...
    import matplotlib.pyplot as plt
    plt.plot(e, d)
    plt.xlabel('energy [eV]')
    plt.ylabel('DOS')
    plt.show()
    plt.savefig("figure_Sidos.png")

