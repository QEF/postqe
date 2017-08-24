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

    # define the Atoms structure reading the xml file
    Si = get_atoms_from_xml_output('Si.xml', schema='../../schemas/qes.xsd')
    # set a simple calculator, only to read the xml file results
    calcul = PostqeCalculator(atoms=None, label='./Si', schema='../../schemas/qes.xsd')
    Si.set_calculator(calcul)
    # read the results
    Si.calc.read_results()
    # TODO: can the above lines be moved into a function?

    # Create a DOS object with width= eV and npts points
    dos = DOS(calcul, width=0.02*13.605698066*2, npts=3000)
    # get the dos and energies for further processing
    d = dos.get_dos()
    e = dos.get_energies()

    # Plot the DOS with Matplotlib...
    import matplotlib.pyplot as plt
    plt.plot(e, d)
    plt.xlabel('energy [eV]')
    plt.ylabel('DOS')
    plt.show()
    plt.savefig("figure_Sidos.png")

