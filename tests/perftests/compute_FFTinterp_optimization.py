#!/usr/bin/env python3
#encoding: UTF-8
import timeit
import numpy as np
import xmlschema
try:
    from cythonplot import FFTinterp1D_Cython, FFTinterp2D_Cython
except ImportError:
    FFTinterp1D_Cython = FFTinterp2D_Cython = None

pi = np.pi


def compute_G(b, nr):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr. G are the vectors in the reciprocal lattice vector.
    b[0], b[1], b[2] are the reciprocal cell base vectors
    """
    G = np.zeros((nr[0], nr[1], nr[2], 3))
    for x in range(0, nr[0]):
        if (x >= nr[0] // 2):
            g0 = (x - nr[0])
        else:
            g0 = x
        for y in range(0, nr[1]):
            if (y >= nr[1] // 2):
                g1 = (y - nr[1])
            else:
                g1 = y
            for z in range(0, nr[2]):
                if (z >= nr[2] // 2):
                    g2 = (z - nr[2])
                else:
                    g2 = z

                G[x, y, z, :] = (g0 * b[0] + g1 * b[1] + g2 * b[2])  # compute the G vector

    return G


def read_charge_file_hdf5(filename, nr):
    """
    Reads a charge file written with QE in HDF5 format. *nr = [nr1,nr2,nr3]* (the dimensions of
    the charge k-points grid) are given as parameter (taken for the xml output file by the caller).

    Notes: In the new format, the values of the charge in the reciprocal space are stored.
    Besides, only the values of the charge > cutoff are stored, together with the Miller indexes.
    Hence
    """
    import h5py

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


def FFTinterp1D(charge, G, a, x0, e1, nx):

        # normalize e1
        m1 = np.linalg.norm(e1)
        if abs(m1) < 1.0E-6:  # if the module is less than 1.0E-6
            e1 = a[1]
            m1 = np.linalg.norm(e1)
        e1 = e1 / m1

        # Computes the FFT of the charge
        fft_charge = np.fft.fftn(charge)
        nr = charge.shape

        # Steps along the e1 direction...
        deltax = m1 / (nx - 1)
        X = np.zeros(nx)
        Y = np.zeros(nx, dtype=complex)

        for i in range(0, nx):
            xi = x0[0] + i * deltax * e1[0]
            yi = x0[1] + i * deltax * e1[1]
            zi = x0[2] + i * deltax * e1[2]

            # For each point, evaluate the charge by Fourier interpolation
            for x in range(0, nr[0]):
                for y in range(0, nr[1]):
                    for z in range(0, nr[2]):
                        arg = 2.0 * pi * (xi * G[x, y, z, 0] + yi * G[x, y, z, 1] + zi * G[x, y, z, 2])
                        Y[i] += fft_charge[x, y, z] * complex(np.cos(arg), np.sin(arg))

            X[i] = i * deltax
            Y[i] = Y[i] / (nr[0] * nr[1] * nr[2])
            #print(X[i], Y[i].real)

        return X, Y


def FFTinterp2D(charge, G, a, x0, e1, e2, nx, ny):
    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1) < 1.0E-6):  # if the module is less than 1.0E-6
        e1 = a[1]
        m1 = np.linalg.norm(e1)
    e1 = e1 / m1

    # normalize e2
    m2 = np.linalg.norm(e2)
    if abs(m2) < 1.0E-6:  # if the module is less than 1.0E-6
        e2 = a[2]
        m2 = np.linalg.norm(e2)
    e2 = e2 / m2

    # Computes the FFT of the charge
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    # Steps along the e1 and e2 directions...
    deltax = m1 / (nx - 1)
    deltay = m2 / (ny - 1)

    temp = np.zeros((nx, ny), dtype=complex)
    X = np.zeros((nx, ny))
    Y = np.zeros((nx, ny))
    Z = np.zeros((nx, ny))

    # loop(s) over the G points
    for x in range(0, nr[0]):
        for y in range(0, nr[1]):
            for z in range(0, nr[2]):

                # eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
                # compute these factors to save CPU time
                eigx = np.zeros(nx, dtype=complex)
                for i in range(0, nx):
                    eigx[i] = np.exp(2.0 * pi * complex(0.0, 1.0) * (i * deltax *
                                                                     (e1[0] * G[x, y, z, 0] + e1[1] * G[x, y, z, 1] +
                                                                      e1[2] * G[x, y, z, 2]) +
                                                                     (x0[0] * G[x, y, z, 0] + x0[1] * G[x, y, z, 1] +
                                                                      x0[2] * G[x, y, z, 2])))

                eigy = np.zeros(ny, dtype=complex)
                for j in range(0, ny):
                    eigy[j] = np.exp(2.0 * pi * complex(0.0, 1.0) * (j * deltax *
                                                                     (e2[0] * G[x, y, z, 0] + e2[1] * G[x, y, z, 1] +
                                                                      e2[2] * G[x, y, z, 2])))

                for i in range(0, nx):
                    for j in range(0, ny):
                        temp[i, j] += fft_charge[x, y, z] * eigx[i] * eigy[j]

    Z = temp.real / (nr[0] * nr[1] * nr[2])

    return X, Y, Z


if __name__ == "__main__":
    pw_output = xmlschema.to_dict("../../examples/example5/Ni.xml", "../../postqe/schemas/qes.xsd")['output']
    nr = np.array([
        pw_output["basis_set"]["fft_grid"]["@nr1"],
        pw_output["basis_set"]["fft_grid"]["@nr2"],
        pw_output["basis_set"]["fft_grid"]["@nr3"]
    ])
    a1 = np.array(pw_output["atomic_structure"]["cell"]["a1"])
    a2 = np.array(pw_output["atomic_structure"]["cell"]["a2"])
    a3 = np.array(pw_output["atomic_structure"]["cell"]["a3"])
    a = np.array([a1, a2, a3])
    b1 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b3"])
    b = np.array([b1, b2, b3])

    ecutrho = (pw_output["basis_set"]["ecutrho"])
    alat = (pw_output["atomic_structure"]["@alat"])

    charge, chargediff = read_charge_file_hdf5("../../examples/example5/tmp/Ni.save/charge-density.hdf5", nr)
    G = compute_G(b, charge.shape)

    x0 = (0,0,0)
    e1 = (1,0,0)
    e2 = (0,1,0)
    nx = 500
    ny = 500

    # Test conformance
#    X, Y = FFTinterp1D(charge, G, a, x0, e1, nx)
#    _X, _Y = FFTinterp1D_Cython(charge, G, a, x0, e1, nx)
#    if not np.array_equal(X, _X) or not np.array_equal(Y, _Y) :
#        raise ValueError("FFTinterp1D_Cython() result is different!")

    # Running performance tests
    print("##### 'FFTinterp1D' versions performance timing #####\n")

    setup = ("from __main__ import a, charge, G, x0, e1, nx, FFTinterp1D, FFTinterp1D_Cython,"
             "FFTinterp2D, FFTinterp2D_Cython")

    if FFTinterp1D_Cython is not None:
        print("FFTinterp1D_Cython:",
              timeit.repeat('FFTinterp1D_Cython(charge, G, a, x0, e1, nx)', setup=setup, number=1, repeat=3))
    print("FFTinterp1D:",
          timeit.repeat('FFTinterp1D(charge, G, a, x0, e1, nx)', setup=setup, number=1, repeat=3))

    print("##### 'FFTinterp2D' versions performance timing #####\n")

    if FFTinterp2D_Cython is not None:
        print("FFTinterp2D_Cython:",
              timeit.repeat('FFTinterp2D_Cython(charge, G, a, x0, e1, e2, nx, ny)', setup=setup, number=1, repeat=3))
    print("FFTinterp1D:",
          timeit.repeat('FFTinterp2D(charge, G, a, x0, e1, e2, nx, ny)', setup=setup, number=1, repeat=3))

