#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of plotting a 2D section (dim=2) of the electronic charge density
using the 'compute_charge' API (see example4.py for the list of this API parameters options)
"""

if __name__ == "__main__":
    from postqe import compute_charge

    charge, figure = compute_charge(prefix='Silicon', outdir='outdir', schema='../../schemas/qes.xsd', fileout='charge.out',
                            x0=(0., 0., 0.), e1=(1., 0., 0.), nx=50, e2=(0., 1., 0.), ny=50, radius=10, dim=2,
                                ifmagn='total', plot_file='charge_plot', method='FFT', format='gnuplot', show=False)

    figure.savefig('2D_charge_plot.png', format='png')
