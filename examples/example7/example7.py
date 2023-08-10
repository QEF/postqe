#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example of plotting a 2D section of the electronic charge density.
"""

if __name__ == "__main__":

    from postqe import get_charge

    charge = get_charge(prefix='Silicon', outdir='outdir', schema='../../schemas/qes.xsd')

    figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=10, dim=1, ifmagn='total',
                 plot_file='postqe1Dsphericalaveragegnuplot',
                         method='spherical', format='gnuplot', show=False)

    figure = charge.plot(x0=(0, 0, 0), e1=(1, 0, 0), nx=10, dim=1, ifmagn='total',
                 plot_file='postqe1DFFTgnuplot', method='FFT',
                         format='gnuplot', show = False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
                plot_file='postqe2DFFTgnuplot', method='FFT', format='gnuplot',
                        show=False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
                plot_file='postqe2DFFTcontour', method='FFT',
                         format='contour.x', show = False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
                plot_file='postqe2DFFTplotrho', method='FFT',
                         format='plotrho.x', show=False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, dim=2, ifmagn='total',
                plot_file='postqe2DFFTxsf', method='FFT', format='xsf',
                         show=False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, radius = 10, dim=2, ifmagn='total',
                plot_file='postqe2Dpolargnuplot', method='polar',
                         format='gnuplot', show = False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, radius = 10, dim=2, ifmagn='total',
                plot_file='postqe2Dpolarcontour', method='polar',
                         format='contour.x', show = False)

    figure = charge.plot(x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 10, e2 = (0., 1., 0.), ny=10, e3 = (0., 0., 1.), nz=10,
                   dim=3, ifmagn='total', plot_file='postqe3DFFTxsf', method='FFT', format='xsf', show=False)

    figure = charge.plot(x0=(0., 0., 0.), e1=(1., 0., 0.), nx=10, e2=(0., 1., 0.), ny=10, e3=(0., 0., 1.), nz=10,
                         dim=3, ifmagn='total', plot_file='postqe3DFFTcube',
                         method='FFT', format='cube', show = False)
