#!/usr/bin/env python3
#encoding: UTF-8

"""
This is a simple example of plotting a 2D section of the electronic charge density.
"""
    
if __name__ == "__main__":

    from postqe import get_charge

    charge = get_charge(label="./Ni", schema='../../schemas/qes.xsd')

    figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), e2=(0, 1, 0), nx=50, ny=50, dim=2)
    figure.savefig("figure_1.png")
    figure.savefig("figure_1.pdf", format='pdf')
