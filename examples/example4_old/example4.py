#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of plotting a 1D section of the electronic charge density.
"""

if __name__ == "__main__":

    # import sys
    # from os import path
    # sys.path.insert(0,'/home/pietro/postqe')

    from postqe import get_charge

    charge = get_charge(prefix='Ni', schema='../../schemas/qes.xsd')
    charge.write('outputcharge.dat')
    figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=100)
    figure.savefig('figure_1.png')
