#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the the band structure of silicon.
"""
    
if __name__ == "__main__":
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

    print(Si.calc.get_bz_k_points())
    bs = Si.calc.band_structure(0)
    #print(bs.energies)
    bs.plot(emin=-20, emax=50, filename='si.png')





