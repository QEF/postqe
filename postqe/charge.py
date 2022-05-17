#
# Copyright (c), 2016-2019, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
import numpy as np
import h5py
from .plot import plot_1d_charge, plot_2d_charge, plot_3d_charge
from .compute_vs import compute_G, compute_v_bare, compute_v_h, compute_v_xc

####### TO BE MOVED TO QESCHEMA 1.2 #######

def read_charge_file_hdf5(filename):
    """
    Reads a PW charge file in HDF5 format.

    :param filename: 
    :return: a dictionary describing the content of file \
    keys=[nr, ngm_g, gamma_only, rhog_, MillerIndexes]
    """
    with h5py.File(filename, "r") as h5f:
        MI = h5f.get('MillerIndices')[:]
        nr1 = 2*max(abs(MI[:, 0]))+1
        nr2 = 2*max(abs(MI[:, 1]))+1
        nr3 = 2*max(abs(MI[:, 2]))+1
        nr = np.array([nr1, nr2, nr3])
        res = dict(h5f.attrs.items())
        res.update({'MillInd':MI,'nr_min': nr})
        rhog = h5f['rhotot_g'][:].reshape(res['ngm_g'],2).dot([1.e0,1.e0j])
        res.update({'rhotot_g':rhog})
        if 'rhodiff_g' in h5f.keys():
            rhog = h5f['rhodiff_g'][:].reshape(res['ngm_g'], 2).dot([1.e0,1.e0j])
            res.update({'rhodiff_g':rhog})
        return res


def get_minus_indexes(g1, g2, g3):
    """
    Used for getting the corresponding minus Miller Indexes. It is meant to be used 
    for converting Gamma Trick grids and is defined only for the for i >=0, in the i =0 plan
     is defined only for j >=0 and when i=0 j=0 k must be >=0. Out of this domain returns 
     None.

    :param g1: rank 1 array containing first Miller Index
    :param g2: rank 1 array containing second Miller Index
    :param g3: rank 1 array containing third Miller Index 
    :return: a rank 2 array with dimension (ngm/2,3) containing mirrored Miller indexes
    """

    def scalar_func(i, j, k):
        """
        scalar function to be vectorized
        :param i: 1st Miller Index
        :param j: 2nd 
        :param k: 3rd
        :return: the mirrored mirror indexes
        """
        if i > 0:
            return -i, j, k
        elif i == 0 and j > 0:
            return 0, -j, k
        elif i == 0 and j == 0 and k > 0:
            return 0, 0, -k
        else:
            return i, j, k

    vector_func = np.vectorize(scalar_func)

    res = np.array(vector_func(g1, g2, g3))
    return res.transpose()


def get_charge_r(filename, nr=None):
    """
    Reads a charge file written with QE in HDF5 format. *nr = [nr1,nr2,nr3]* (the dimensions of
    the charge k-points grid) are given as parameter (taken for the xml output file by the caller).

    Notes: In the new format, the values of the charge in the reciprocal space are stored.
    Besides, only the values of the charge > cutoff are stored, together with the Miller indexes.
    Hence
    """

    cdata = read_charge_file_hdf5(filename)
    if nr is None:
        nr1, nr2, nr3 = cdata['nrmin']
    else:
        nr1, nr2, nr3 = nr
    gamma_only = 'TRUE' in str(cdata['gamma_only']).upper()
    # Load the total charge
    rho_temp = np.zeros([nr1, nr2, nr3], dtype=np.complex128)
    for (i, j, k),rho in zip( cdata['MillInd'],cdata['rhotot_g']):
        try:
            rho_temp[i, j, k]=rho
        except IndexError:
            pass

    if gamma_only:
        rhotot_g = cdata['rhotot_g'].conjugate()
        MI = get_minus_indexes(cdata['MillInd'][:,0], cdata['MillInd'][:,1], cdata['MillInd'][:,2])
        print("MI", MI)
        for (i, j, k), rho  in zip(MI, rhotot_g):
            try:
                rho_temp[i, j, k] = rho
            except IndexError:
                pass

    rhotot_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3


    # Read the charge difference spin up - spin down if present (for magnetic calculations)
    if 'rhodiff_g' in cdata.keys():
        rho_temp = np.zeros([nr1, nr2, nr3], dtype=np.complex128)
        for (i,j,k), rho  in zip(cdata['MillInd'], cdata['rhodiff_g']):
            try:
                rho_temp[i, j, k] = rho
            except IndexError:
                pass
        if gamma_only:
            rhodiff_g = cdata['rhodiff_g'].conjugate()
            for (i, j, k),rho in zip(MI, rhodiff_g):
                try:
                    rho_temp[i, j, k] = rho
                except IndexError:
                    pass

        rhodiff_r = np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3
        return rhotot_r.real, rhodiff_r.real
    else:
        return rhotot_r.real, None


def charge_r_from_cdata(cdata, MI, gamma_only, nr ): 
    """
    Computes density in real space from data"
    :cdata: complex coefficients for 3D grid. correspondi G vectors are read from MI 
    :MI:    integer array of dim=3 with the corresponding indexs for data in the FFT 3D grid
    :gamma_only: boolean, if true data are real in real space, complex data for -G  are implicitly provided as conjg(rho(G)) 
    :nr:  integer array with the 3 dimensions on the FFT grid 
    :returns: 3D data in real space. 
    """
    nr1,nr2,nr3 = nr 
    rho_temp = np.zeros([nr1, nr2, nr3], dtype=np.complex128) 
    for (i,j,k),rho in zip ( MI, cdata):
        try:
            rho_temp [i,j,k] = rho
        except IndexError:
            pass 
    #
    if gamma_only:
        MI_minus = (get_minus_indexes(_) for _ in MI )
        rhog = ( _.conjugate() for _ in cdata) 
        for (i,j,k),rho in zip ( MI_minus, rhog): 
            try:
                rho_temp[i,j,k] = rhog 
            except IndexError:
                pass
    return np.fft.ifftn(rho_temp) * nr1 * nr2 * nr3 



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
                if count % 5 == 0:
                    fout.write("\n")

    fout.close()

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
####### TO BE MOVED TO QESCHEMA #######



class Charge:
    """
    A class for charge density (and data in real 3D grids) 
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
        """
        if not nr:
            try:
                nr = self.nr
            except:
                raise AttributeError("nr not defined in this Charge object")
        charge, charge_diff = get_charge_r(filename, np.array(nr))
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
        :param radius: radius of the sphere in the polar average method
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
            'a': self.calculator.get_a_vectors(),
            'b': self.calculator.get_b_vectors(),
            'alat': self.calculator.get_alat(),
            'nat': len(self.calculator.get_atomic_positions()),
            'atomic_positions': self.calculator.get_atomic_positions(),
            'atomic_species': self.calculator.get_atomic_species(),
        }
        G = compute_G(struct_info['b'], self.nr)

        if not self.calculator.get_spin_polarized():  # non magnetic calculation
            if dim == 1:    # 1D section ylab='charge', plot_file='', format='', method='FFT'
                fig = plot_1d_charge(self.charge, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format)
            elif dim == 2:  # 2D section
                fig = plot_2d_charge(
                    self.charge, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format
                )
            else:           # 3D section
                fig = plot_3d_charge(
                    self.charge, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method, format
                )
        else:  # magnetic calculation, plot as ifmagn
            if ifmagn == 'up':
                charge_up = (self.charge + self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot_1d_charge(
                        charge_up, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format
                    )
                elif dim == 2:  # 2D section
                    fig = plot_2d_charge(
                        charge_up, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format
                    )
                else:  # 3D section
                    fig = plot_3d_charge(
                        charge_up, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method, format
                    )
            elif ifmagn == 'down':
                charge_down = (self.charge - self.charge_diff) / 2.0
                if dim == 1:  # 1D section
                    fig = plot_1d_charge(
                        charge_down, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format
                    )
                elif dim == 2:  # 2D section
                    fig = plot_2d_charge(
                        charge_down, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format
                    )
                else:  # 3D section
                    fig = plot_3d_charge(
                        charge_down, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file, method, format
                    )
            else:
                if dim == 1:  # 1D section ylab='charge', plot_file='', format='', method='FFT'
                    fig = plot_1d_charge(
                        self.charge, G, struct_info, x0, e1, nx, 'charge', plot_file, method, format
                    )
                elif dim == 2:  # 2D section
                    fig = plot_2d_charge(
                        self.charge, G, struct_info, x0, e1, e2, nx, ny, radius, 'charge', plot_file, method, format
                    )
                else:  # 3D section
                    fig = plot_3d_charge(self.charge, G, struct_info, x0, e1, e2, e3, nx, ny, nz, 'charge', plot_file,
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
        """Call the Charge constructor, define self.pot_type """
        self.setvars(*args, **kwargs)
        try:
            self.pot_type = pot_type
        except:
            self.pot_type = 'v_tot'

    def write(self, filename):
        """
        Write the potential in a text file. The potential must have been calculated before.
        :param filename: name of the output file
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
        an instance of the class Potential is create (default 'v_tot').
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
            self.v = compute_v_bare(
                ecutrho, alat, a, self.nr, atomic_positions, atomic_species, pseudodir
            )
        elif self.pot_type=='v_h':
            self.v = compute_v_h(self.charge, ecutrho, alat, b)
        elif self.pot_type=='v_xc':
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            self.v = compute_v_xc(self.charge, charge_core, str(functional))
        elif self.pot_type=='v_tot':
            v_bare = compute_v_bare(
                ecutrho, alat, a[0], a[1], a[2], self.nr, atomic_positions, atomic_species, pseudodir
            )
            v_h =  compute_v_h(self.charge, ecutrho, alat, b)
            # TODO: core charge to be implemented
            charge_core = np.zeros(self.nr)
            v_xc = compute_v_xc(self.charge, charge_core, str(functional))
            self.v = v_bare + v_h + v_xc


    def plot(self, x0=(0., 0., 0.), e1=(1., 0., 0.), nx=50, e2=(0., 1., 0.), ny=50, e3=(0., 0., 1.), nz=50,
             radius=1, dim=1, plot_file='', method='FFT', format='gnuplot', show=True):
        """
        Plot a 1D, 2D or 3D section of the potential from x0 along e1 (e2, e3) direction(s) using Fourier interpolation
        or another method (see below). For 1D or 2D sections, the code produce a Matplotlib plot. For a 3D plot, the
        charge must be exported in 'plotfile' with a suitable format ('xsf' or 'cube') and can be visualized with
        the corresponding external codes.

        :param x0: 3D vector (a tuple), origin of the line
        :param e1, e2, e3: 3D vectors (tuples) which determines the plotting lines
        :param nx, ny, nz: number of points along e1, e2, e3
        :param radius: radius of the sphere in the polar average method
        :param dim: 1, 2, 3 for a 1D, 2D or 3D section respectively
        :param plot_file: file where plot data are exported in the chosen format (Gnuplot, XSF, cube Gaussian, etc.)
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
            'a': self.calculator.get_a_vectors(),
            'b': self.calculator.get_b_vectors(),
            'alat': self.calculator.get_alat(),
            'nat': len(self.calculator.get_atomic_positions()),
            'atomic_positions': self.calculator.get_atomic_positions(),
            'atomic_species': self.calculator.get_atomic_species(),
        }
        G = compute_G(struct_info['b'], self.nr)

        if dim == 1:    # 1D section ylab='charge', plot_file='', format='', method='FFT'
            fig = plot_1d_charge(self.v, G, struct_info, x0, e1, nx, self.pot_type, plot_file, method, format)
        elif dim == 2:  # 2D section
            fig = plot_2d_charge(
                self.v, G, struct_info, x0, e1, e2, nx, ny, radius, self.pot_type, plot_file, method, format
            )
        else:           # 3D section
            fig = plot_3d_charge(
                self.v, G, struct_info, x0, e1, e2, e3, nx, ny, nz, self.pot_type, plot_file, method, format
            )

        if dim < 3:
            if show is True:
                fig.show()
            return fig
        else:
            return None
