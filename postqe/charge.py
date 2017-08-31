#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
import h5py
from .plot import plot1D_FFTinterp, plot2D_FFTinterp
from .compute_vs import compute_G, compute_v_bare, compute_v_h, compute_v_xc

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


def write_charge(filename, charge, header):
    """
    Write the charge or another quantity calculated by postqe into a text file *filename*.
    """

    fout = open(filename, "w")

    # The header contains some information on the system, the grid nr, etc.
    fout.write(header)
    nr = charge.shape
    count = 0
    # Loop with the order as in QE files
    for z in range(0, nr[2]):
        for y in range(0, nr[1]):
            for x in range(0, nr[0]):
                fout.write("  {:.9E}".format(charge[x, y, z]))
                count += 1
                if (count % 5 == 0):
                    fout.write("\n")

    fout.close()


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


    def write(self, filename):
        header='# Charge file\n'
        header+='# nr1= '+str(self.nr[0])+' nr2= '+str(self.nr[1])+' nr3= '+str(self.nr[2])+'\n'
        try:
            self.charge
            write_charge(filename, self.charge, header)
        except:
            pass
        if self.calculator.get_spin_polarized():  # non magnetic calculation
            charge_up = (self.charge + self.charge_diff) / 2.0
            charge_down = (self.charge - self.charge_diff) / 2.0
            write_charge(filename + '_up', charge_up, header)
            write_charge(filename + '_down', charge_down, header)


    def plot(self, x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, e2 = (1., 0., 0.), ny=50, dim=1, ifmagn='total'):
        """
        Plot a 1D or 2D section of the charge from x0 along e1 (e2) direction(s) using Fourier interpolation.

        :param x0: 3D vector, origin of the line
        :param e1, e2: 3D vectors which determines the plotting lines
        :param nx, ny: number of points along e1, e2
        :param dim: 1 for a 1D section, 2 for a 2D section
        :param ifmagn: for a magnetic calculation, 'total' plot the total charge, 'up' plot the charge with spin up, 'down' for spin down
        :return: a Matplotlib figure object
        """
        try:
            self.charge
        except:
            return
        a = self.calculator.get_a_vectors()
        b = self.calculator.get_b_vectors()
        G = compute_G(b, self.nr)

        if not self.calculator.get_spin_polarized():  # non magnetic calculation
            if dim == 1:  # 1D section
                fig = plot1D_FFTinterp(self.charge, G, a, x0, e1, nx)
            else:
                fig = plot2D_FFTinterp(self.charge, G, a, x0, e1, e2, nx, ny)
            fig.show()
            return fig
        else:  # magnetic calculation, plot as ifmagn
            if ifmagn == 'up':
                charge_up = (self.charge + self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot1D_FFTinterp(charge_up, G, a, x0, e1, nx)
                else:
                    fig = plot2D_FFTinterp(charge_up, G, a, x0, e1, e2, nx, ny)
                fig.show()
            elif ifmagn == 'down':
                charge_down = (self.charge - self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot1D_FFTinterp(charge_down, G, a, x0, e1, nx)
                else:
                    fig = plot2D_FFTinterp(charge_down, G, a, x0, e1, e2, nx, ny)
                fig.show()
            else:
                if dim == 1:  # 1D section
                    fig = plot1D_FFTinterp(self.charge, G, a, x0, e1, nx)
                else:
                    fig = plot2D_FFTinterp(self.charge, G, a, x0, e1, e2, nx, ny)
                fig.show()
            return fig


class Potential(Charge):
    """
    A class for a potential. This is derived from a Charge class and additionally contains the potential.
    """
    def __init__(self, *args, **kwargs):
        """Call the Charge constructor """
        self.setvars(*args, **kwargs)

    def write(self, filename):
        header='# Potential file '+self.pot_type+'\n'
        header+='# nr1= '+str(self.nr[0])+' nr2= '+str(self.nr[1])+' nr3= '+str(self.nr[2])+'\n'
        try:
            self.v
            write_charge(filename, self.v, header)
        except:
            pass

    def compute_potential(self, pot_type='v_tot'):
        try:
            self.charge
        except:
            return
        self.pot_type = pot_type
        alat = self.calculator.get_alat()
        ecutrho = self.calculator.get_ecutrho()
        functional = self.calculator.get_xc_functional()
        a = self.calculator.get_a_vectors()
        b = self.calculator.get_b_vectors()
        atomic_positions = self.calculator.get_atomic_positions()
        atomic_species = self.calculator.get_atomic_species()
        pseudodir = self.calculator.get_pseudodir()

        if (pot_type=='v_bare'):
            self.v = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], self.nr, atomic_positions, atomic_species, pseudodir)
        elif (pot_type=='v_h'):
            self.v = compute_v_h(self.charge, ecutrho, alat, b)
        elif (pot_type=='v_xc'):
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            self.v = compute_v_xc(self.charge, charge_core, str(functional))
        elif (pot_type=='v_tot'):
            v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], self.nr, atomic_positions, atomic_species, pseudodir)
            v_h =  compute_v_h(self.charge, ecutrho, alat, b)
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            v_xc = compute_v_xc(self.charge, charge_core, str(functional))
            self.v = v_bare + v_h + v_xc

    def plot(self, x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, e2 = (1., 0., 0.), ny=50, dim=1, ifmagn='total'):
        """
        Plot a 1D or 2D section of the charge from x0 along e1 (e2) direction(s) using Fourier interpolation.

        :param x0: 3D vector, origin of the line
        :param e1, e2: 3D vectors which determines the plotting lines
        :param nx, ny: number of points along e1, e2
        :param dim: 1 for a 1D section, 2 for a 2D section
        :param ifmagn: for a magnetic calculation, 'total' plot the total charge, 'up' plot the charge with spin up, 'down' for spin down
        :return: a Matplotlib figure object
        """
        try:
            self.v
        except:
            return
        a = self.calculator.get_a_vectors()
        b = self.calculator.get_b_vectors()
        G = compute_G(b, self.nr)

        if dim == 1:  # 1D section
            fig = plot1D_FFTinterp(self.v, G, a, x0, e1, nx)
        else:
            fig = plot2D_FFTinterp(self.v, G, a, x0, e1, e2, nx, ny)
        fig.show()
        return fig
