#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of plotting a 1D section of the electronic charge density.
"""
    
if __name__ == "__main__":

    from postqe import get_charge

    charge = get_charge(label="./Ni", schema='../../schemas/qes.xsd')
    charge.write('outputcharge.dat')

    figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=100)
    figure.savefig("figure_1.png")
