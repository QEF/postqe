#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
import h5py
from .plot import plot_1Dcharge, plot_2Dcharge, plot_3Dcharge
from .compute_vs import compute_G, compute_v_bare, compute_v_h, compute_v_xc


def read_charge_file_hdf5(filename, nr = None):
    """
    Reads a charge file written with QE in HDF5 format. *nr = [nr1,nr2,nr3]* (the dimensions of
    the charge k-points grid) are given as parameter (taken for the xml output file by the caller).

    Notes: In the new format, the values of the charge in the reciprocal space are stored.
    Besides, only the values of the charge > cutoff are stored, together with the Miller indexes.
    Hence
    """

    with h5py.File(filename, "r") as h5f:
        MI = h5f.get('MillerIndices')
        if nr is None:
            nr1 = 2*max(abs(MI[:,0]))+1
            nr2 = 2*max(abs(MI[:,1]))+1
            nr3 = 2*max(abs(MI[:,2]))+1
        else:
            nr1,nr2,nr3 = nr

        ngm_g = h5f.attrs.get('ngm_g')
        gamma_only = h5f.attrs.get('gamma_only')
        gamma_only = 'TRUE' in str(gamma_only).upper()
        # Read the total charge
        aux = np.array(h5f['rhotot_g']).reshape([ngm_g,2])
        rhotot_g = aux.dot([1.e0,1.e0j])
        rho_temp = np.zeros([nr1,nr2,nr3],dtype=np.complex128)
        for el in zip( MI,rhotot_g):
            (i,j,k), rho = el
            try:
                rho_temp[i,j,k]=rho
            except IndexError:
                pass
        if gamma_only:
            rhotot_g = aux.dot([1.e0,-1.e0j])
            for el in zip(MI, rhotot_g):
                (i,j,k), rho = el
                if i > 0:
                    try:
                        rho_temp[-i, j, k] = rho
                    except IndexError:
                        pass
                elif j > 0:
                    try:
                        rho_temp[0, -j, k] = rho
                    except IndexError:
                        pass
                elif k > 0:
                    try:
                        rho_temp[0, 0, -k] = rho
                    except IndexError:
                        pass

        rhotot_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3

        # Read the charge difference spin up - spin down if present (for magnetic calculations)
        if not h5f.get('rhodiff_g') is None:
            aux = np.array(h5f['rhodiff_g']).reshape([ngm_g, 2])
            rhodiff_g = aux.dot([1.0e0,1.0e0j])
            rho_temp = np.zeros([nr1, nr2, nr3], dtype=np.complex128)
            for el in zip(h5f['MillerIndices'], rhodiff_g):
                (i, j, k), rho = el
                try:
                    rho_temp[i, j, k] = rho
                except IndexError:
                    pass
            if gamma_only:
                rhodiff_g = aux.dot([1.e0, -1.e0j])
                for el in zip(MI, rhodiff_g):
                    (i, j, k), rho = el
                    if i > 0:
                        try:
                            rho_temp[-i, j, k] = rho
                        except IndexError:
                            pass
                    elif j > 0:
                        try:
                            rho_temp[0, -j, k] = rho
                        except IndexError:
                            pass
                    elif k > 0:
                        try:
                            rho_temp[0, 0, -k] = rho
                        except IndexError:
                            pass
            rhodiff_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3
            return rhotot_r.real, rhodiff_r.real
        else:
            return rhotot_r.real, None


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


    def plot(self, x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, e2 = (0., 1., 0.), ny=50, e3 = (0., 0., 1.), nz=50,
             radius=1, dim=1, ifmagn='total', plot_file='', method='FFT', format='gnuplot', show=True):
        """
        Plot a 1D, 2D or 3D section of the charge from x0 along e1 (e2, e3) direction(s) using Fourier interpolation
        or another method (see below). For 1D or 2D sections, the code produce a Matplotlib plot. For a 3D plot, the
        charge must be exported in 'plotfile' with a suitable format ('xsf' or 'cube') and can be visualized with
        the corresponding external codes.

        :param x0: 3D vector (a tuple), origin of the line
        :param e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
        :param nx, ny, nz: number of points along e1, e2, e3
        :param radius: radious of the sphere in the polar average method
        :param dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
        :param ifmagn: for a magnetic calculation, 'total' plot the total charge, 'up' plot the charge with spin up,
                       'down' for spin down
        :param plotfile: file where plot data are exported in the chosen format (Gnuplot, XSF, cube Gaussian, etc.)
        :param method: interpolation method. Available choices are:\n
                        'FFT' -> Fourier interpolation (default)
                        'polar' -> 2D polar plot on a sphere
                        'spherical' -> 1D plot of the spherical average
                        'splines' -> not implemented
        :param format: format of the (optional) exported file. Available choices are:\n
                        'gnuplot' -> plain text format for Gnuplot (default). Available for 1D and 2D sections.
                        'xsf' -> XSF format for the XCrySDen program. Available for 2D and 3D sections.
                        'cube' -> cube Gaussian format. Available for 3D sections.
                        'contour' -> format for the contour.x code of Quantum Espresso
                        'plotrho' -> format for the plotrho.x code of Quantum Espresso
        :param show: if True, show the Matplotlib plot (only for 1D and 2D sections)
        :return: a Matplotlib figure object for 1D and 2D sections, None for 3D sections
        """
        # TODO: implement a Matplotlib plot for polar 2D
        try:
            self.charge
        except:
            return None
        # Extract some structural info in a dictionary
        struct_info = {
            'a' : self.calculator.get_a_vectors(),
            'b' : self.calculator.get_b_vectors(),
            'alat' : self.calculator.get_alat(),
            'nat'  : len(self.calculator.get_atomic_positions()),
            'atomic_positions' : self.calculator.get_atomic_positions(),
            'atomic_species': self.calculator.get_atomic_species(),
        }
        G = compute_G(struct_info['b'], self.nr)

        if not self.calculator.get_spin_polarized():  # non magnetic calculation
            if dim == 1:    # 1D section ylab='charge', plot_file='', format='', method='FFT'
                fig = plot_1Dcharge(self.charge, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format)
            elif dim == 2:  # 2D section
                fig = plot_2Dcharge(self.charge, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format)
            else:           # 3D section
                fig = plot_3Dcharge(self.charge, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method, format)
        else:  # magnetic calculation, plot as ifmagn
            if ifmagn == 'up':
                charge_up = (self.charge + self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot_1Dcharge(charge_up, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format)
                elif dim == 2:  # 2D section
                    fig = plot_2Dcharge(charge_up, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format)
                else:  # 3D section
                    fig = plot_3Dcharge(charge_up, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method,
                                        format)
            elif ifmagn == 'down':
                charge_down = (self.charge - self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot_1Dcharge(charge_down, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format)
                elif dim == 2:  # 2D section
                    fig = plot_2Dcharge(charge_down, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format)
                else:  # 3D section
                    fig = plot_3Dcharge(charge_down, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method,
                                        format)
            else:
                if dim == 1:  # 1D section ylab='charge', plot_file='', format='', method='FFT'
                    fig = plot_1Dcharge(self.charge, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format)
                elif dim == 2:  # 2D section
                    fig = plot_2Dcharge(self.charge, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method,
                                        format)
                else:  # 3D section
                    fig = plot_3Dcharge(self.charge, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file,
                                        method, format)

        if dim < 3:
            if show == True:
                fig.show()
            return fig
        else:
            return None


class Potential(Charge):
    """
    A class for a potential. This is derived from a Charge class and additionally contains the potential.
    """
    def __init__(self, *args, pot_type='v_tot', **kwargs):
        """Call the Charge constructor """
        self.setvars(*args, **kwargs)
        try:
            self.pot_type = pot_type
        except:
            self.pot_type = 'v_tot'
        print(self.pot_type)

    def write(self, filename):
        """
        Write the potential in a text file. The potential must have been calculated before.
        :param filename: name of the output file
        :return:
        """
        header='# Potential file '+self.pot_type+'\n'
        header+='# nr1= '+str(self.nr[0])+' nr2= '+str(self.nr[1])+' nr3= '+str(self.nr[2])+'\n'
        try:
            self.v
            write_charge(filename, self.v, header)
        except:
            pass

    def compute_potential(self):
        """
        Compute the potential from the electronic charge. The type of potential is defined in self.pot_type when
        an instance of the class Potential is create (default 'v_tot')
        :return:
        """
        try:
            self.charge
        except:
            return
        alat = self.calculator.get_alat()
        ecutrho = self.calculator.get_ecutrho()
        functional = self.calculator.get_xc_functional()
        a = self.calculator.get_a_vectors()
        b = self.calculator.get_b_vectors()
        atomic_positions = self.calculator.get_atomic_positions()
        atomic_species = self.calculator.get_atomic_species()
        pseudodir = self.calculator.get_pseudodir()

        if self.pot_type=='v_bare':
            self.v = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], self.nr, atomic_positions, atomic_species, pseudodir)
        elif self.pot_type=='v_h':
            self.v = compute_v_h(self.charge, ecutrho, alat, b)
        elif self.pot_type=='v_xc':
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            self.v = compute_v_xc(self.charge, charge_core, str(functional))
        elif self.pot_type=='v_tot':
            v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], self.nr, atomic_positions, atomic_species, pseudodir)
            v_h =  compute_v_h(self.charge, ecutrho, alat, b)
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            v_xc = compute_v_xc(self.charge, charge_core, str(functional))
            self.v = v_bare + v_h + v_xc


    def plot(self, x0 = (0., 0., 0.), e1 = (1., 0., 0.), nx = 50, e2 = (0., 1., 0.), ny=50, e3 = (0., 0., 1.), nz=50,
             radius=1, dim=1, plot_file='', method='FFT', format='gnuplot', show=True):
        """
        Plot a 1D, 2D or 3D section of the potential from x0 along e1 (e2, e3) direction(s) using Fourier interpolation
        or another method (see below). For 1D or 2D sections, the code produce a Matplotlib plot. For a 3D plot, the
        charge must be exported in 'plotfile' with a suitable format ('xsf' or 'cube') and can be visualized with
        the corresponding external codes.

        :param x0: 3D vector (a tuple), origin of the line
        :param e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
        :param nx, ny, nz: number of points along e1, e2, e3
        :param radius: radious of the sphere in the polar average method
        :param dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
        :param plotfile: file where plot data are exported in the chosen format (Gnuplot, XSF, cube Gaussian, etc.)
        :param method: interpolation method. Available choices are:\n
                        'FFT' -> Fourier interpolation (default)
                        'polar' -> 2D polar plot on a sphere
                        'spherical' -> 1D plot of the spherical average
                        'splines' -> not implemented
        :param format: format of the (optional) exported file. Available choices are:\n
                        'gnuplot' -> plain text format for Gnuplot (default). Available for 1D and 2D sections.
                        'xsf' -> XSF format for the XCrySDen program. Available for 2D and 3D sections.
                        'cube' -> cube Gaussian format. Available for 3D sections.
                        'contour' -> format for the contour.x code of Quantum Espresso
                        'plotrho' -> format for the plotrho.x code of Quantum Espresso
        :param show: if True, show the Matplotlib plot (only for 1D and 2D sections)
        :return: a Matplotlib figure object for 1D and 2D sections, None for 3D sections
        """
        # TODO: implement a Matplotlib plot for polar 2D
        try:
            self.v
        except:
            self.compute_potential()

        # Extract some structural info in a dictionary
        struct_info = {
            'a' : self.calculator.get_a_vectors(),
            'b' : self.calculator.get_b_vectors(),
            'alat' : self.calculator.get_alat(),
            'nat'  : len(self.calculator.get_atomic_positions()),
            'atomic_positions' : self.calculator.get_atomic_positions(),
            'atomic_species': self.calculator.get_atomic_species(),
        }
        G = compute_G(struct_info['b'], self.nr)

        if dim == 1:    # 1D section ylab='charge', plot_file='', format='', method='FFT'
            fig = plot_1Dcharge(self.v, G, struct_info, x0, e1, nx, self.pot_type, plot_file, method, format)
        elif dim == 2:  # 2D section
            fig = plot_2Dcharge(self.v, G, struct_info, x0, e1, e2, nx, ny, radius, self.pot_type, plot_file, method, format)
        else:           # 3D section
            fig = plot_3Dcharge(self.v, G, struct_info, x0, e1, e2, e3, nx, ny, nz, self.pot_type, plot_file, method, format)

        if dim < 3:
            if show == True:
                fig.show()
            return fig
        else:
            return None
