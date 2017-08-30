#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the electronic density of states (DOS or edos).

Computes the DOS for silicon (file Si.xml)
"""
    
if __name__ == "__main__":
    from postqe import get_dos

    # get a DOS object
    dos = get_dos(label="./Si", schema='../../schemas/qes.xsd', width=0.5, npts=200)

    # get the dos and energies for further processing
    d = dos.get_dos()
    e = dos.get_energies()

    # Plot the DOS with Matplotlib...
    import matplotlib.pyplot as plt
    plt.plot(e, d)
    plt.xlabel('energy [eV]')
    plt.ylabel('DOS')
    plt.savefig("figure.png")
    plt.show()


