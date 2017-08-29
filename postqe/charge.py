#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
import h5py
from .compute_vs import compute_G
from .plot import plot1D_FFTinterp, plot2D_FFTinterp

def read_charge_file_hdf5(filename, nr):
    """
    Reads a charge file written with QE in HDF5 format. *nr = [nr1,nr2,nr3]* (the dimensions of
    the charge k-points grid) are given as parameter (taken for the xml output file by the caller).

    Notes: In the new format, the values of the charge in the reciprocal space are stored.
    Besides, only the values of the charge > cutoff are stored, together with the Miller indexes.
    Hence
    """

    nr1, nr2, nr3 = nr
    with h5py.File(filename, "r") as h5f:
        ngm_g = h5f.attrs.get('ngm_g')
        # Read the total charge
        aux = np.array(h5f['rhotot_g']).reshape([ngm_g,2])
        rhotot_g = np.array(list(map(lambda x: x.dot((1e0,1.j)), aux)))
        rho_temp = np.zeros([nr1,nr2,nr3],dtype=np.complex128)
        for el in zip( h5f['MillerIndices'],rhotot_g):
            (i,j,k), rho = el
            rho_temp[i,j,k]=rho
        rhotot_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3

        # Read the charge difference spin up - spin down if present (for magnetic calculations)
        try:
            aux = np.array(h5f['rhodiff_g']).reshape([ngm_g, 2])
            rhodiff_g = np.array(list(map(lambda x: x.dot((1e0, 1.j)), aux)))
            rho_temp = np.zeros([nr1, nr2, nr3], dtype=np.complex128)
            for el in zip(h5f['MillerIndices'], rhodiff_g):
                (i, j, k), rho = el
                rho_temp[i, j, k] = rho
            rhodiff_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3
        except:
            rhodiff_r = np.zeros([nr1, nr2, nr3], dtype=np.complex128)

    return rhotot_r.real, rhodiff_r.real


class Charge:
    """
    A class for charge density.
    """
    # TODO: include the Miller index as in HDF5 file?
    def __init__(self, *args, **kwargs):
        """Create charge object from """
        self.setvars(*args, **kwargs)

    def setvars(self, nr_temp, charge=None, charge_diff=None):
        nr = np.array(nr_temp)
        assert nr.shape[0] == 3
        self.nr = nr
        try:
            assert charge.shape == (nr[0], nr[1], nr[2])
            self.charge = charge
        except:
            pass
        try:
            assert charge_diff.shape == (nr[0], nr[1], nr[2])
            self.charge_diff = charge_diff
        except:
            pass

    def set_calculator(self, calculator):
        self.calculator = calculator

    def read(self, filename, nr=None):
        """
        Read the charge from a HDF5 file.

        :param filename: HDF5 file with charge data
        :param nr: a numpy array or list of length 3 containing the grid dimensions
        :return: nothing
        """
        if not nr:
            try:
                nr = self.nr
            except:
                raise AttributeError("nr not defined in this Charge object")
        charge, charge_diff = read_charge_file_hdf5(filename, np.array(nr))
        self.charge = charge
        self.charge_diff = charge_diff

    # TODO: to be implemented
    def write(self):
        #write_charge(plot_file + '_up', charge_up, header)
        #write_charge(plot_file + '_down', charge_down, header)
        pass

    def plot1D(self, x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, ifmagn='total'):
        """
        Plot a 1D section of the charge from x0 along e1 direction using Fourier interpolation.

        :param x0: 3D vector, origin of the line
        :param e1: 3D vector which determines the plotting line
        :param nx: number of points in the line
        :param ifmagn: for a magnetic calculation, 'total' plot the total charge, 'up' plot the charge with spin up, 'down' for spin down
        :return: a Matplotlib figure object
        """
        a = self.calculator.get_a_vectors()
        b = self.calculator.get_b_vectors()
        G = compute_G(b, self.nr)
        if not self.calculator.get_spin_polarized():  # non magnetic calculation
            fig = plot1D_FFTinterp(self.charge, G, a, x0, e1, nx)
            fig.show()
            return fig
        else:  # magnetic calculation, plot as ifmagn
            if ifmagn=='up':
                charge_up = (self.charge + self.charge_diff) / 2.0
                fig2 = plot1D_FFTinterp(charge_up, G, a, x0, e1, nx)
                fig2.show()
            elif ifmagn=='down':
                charge_down = (self.charge - self.charge_diff) / 2.0
                fig = plot1D_FFTinterp(charge_down, G, a, x0, e1, nx)
                fig.show()
            else:
                fig = plot1D_FFTinterp(self.charge, G, a, x0, e1, nx)
                fig.show()
            return fig
