#!/usr/bin/env python3
#encoding: UTF-8
import timeit
import numpy as np
import xmlschema
pi = np.pi


def get_cell_data(xmlfile):
    """
    Gets some data about the unit cell from the xmlfile.

    :param xmlfile:
    :return ibrav, alat, a, b:
    """
    d = xmlschema.to_dict(xmlfile)
    dout = d["output"]
    alat = (dout["atomic_structure"]["@alat"])
    a1 = np.array(dout["atomic_structure"]["cell"]["a1"])
    a2 = np.array(dout["atomic_structure"]["cell"]["a2"])
    a3 = np.array(dout["atomic_structure"]["cell"]["a3"])
    ibrav = (dout["atomic_structure"]["@bravais_index"])
    b1 = np.array(dout["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(dout["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(dout["basis_set"]["reciprocal_lattice"]["b3"])
    a = np.array([a1, a2, a3])
    b = np.array([b1, b2, b3])
    a_p = (dout["atomic_structure"]["atomic_positions"]["atom"])
    a_s = (dout["atomic_species"]["species"])
    nat = (dout["atomic_structure"]["@nat"])
    ntyp = (dout["atomic_species"]["@ntyp"])

    # for subsequent loops it is important to have always lists for atomic_positions
    # and atomic_species. If this is not, convert
    if (type(a_s) == type([])):
        atomic_species = a_s
    else:
        atomic_species = [a_s]

    if (type(a_p) == type([])):
        atomic_positions = a_p
    else:
        atomic_positions = [a_p]

    return ibrav, alat, a, b, nat, ntyp, atomic_positions, atomic_species


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


def compute_G(b,nr):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr. G are the vectors in the reciprocal lattice vector.
    b[0], b[1], b[2] are the reciprocal cell base vectors
    """
    G = np.zeros((nr[0],nr[1],nr[2],3))
    for x in range(0,nr[0]):
        if (x>=nr[0]//2):
            g0 = (x-nr[0])
        else:
            g0 = x
        for y in range(0,nr[1]):
            if (y>=nr[1]//2):
                g1 = (y-nr[1])
            else:
                g1 = y
            for z in range(0,nr[2]):
                if (z>=nr[2]//2):
                    g2 = (z-nr[2])
                else:
                    g2 = z
 
                G[x,y,z,:] = (g0*b[0]+g1*b[1]+g2*b[2])    # compute the G vector
                
    return G


def compute_G_ver2(b, nr):
    """
    A version of compute_G with the slight optimization of calculating
    only once the "nr[x] // 2" thresholds.
    """
    hnr0 = nr[0] // 2
    hnr1 = nr[1] // 2
    hnr2 = nr[2] // 2

    G = np.empty((nr[0], nr[1], nr[2], 3))
    for x in range(nr[0]):
        if x >= hnr0:
            g0 = x - nr[0]
        else:
            g0 = x
        for y in range(nr[1]):
            if y >= hnr1:
                g1 = y - nr[1]
            else:
                g1 = y
            for z in range(nr[2]):
                if z >= hnr2:
                    g2 = z - nr[2]
                else:
                    g2 = z

                G[x, y, z, :] = g0 * b[0] + g1 * b[1] + g2 * b[2]  # compute the G vector

    return G


def compute_G_ver3(b, nr):
    """
    A version of compute_G that uses the np.dot API in the inner cycle.
    """
    hnr0 = nr[0] // 2
    hnr1 = nr[1] // 2
    hnr2 = nr[2] // 2

    G = np.empty((nr[0], nr[1], nr[2], 3))
    for x in range(nr[0]):
        for y in range(nr[1]):
            for z in range(nr[2]):
                g = (
                    x - nr[0] if x >= hnr0 else x,
                    y - nr[1] if y >= hnr1 else y,
                    z - nr[2] if z >= hnr2 else z
                )

                G[x, y, z, :] = np.dot(g, b) # compute the G vector

    return G


def compute_G_ver4(b, nr):
    """
    A version of compute_G that uses np.dot API on entire array.
    """
    hnr0 = nr[0] // 2
    hnr1 = nr[1] // 2
    hnr2 = nr[2] // 2

    g = np.empty((nr[0], nr[1], nr[2], 3))
    for x in range(nr[0]):
        for y in range(nr[1]):
            for z in range(nr[2]):
                g[x, y, z, :] = (
                    x - nr[0] if x >= hnr0 else x,
                    y - nr[1] if y >= hnr1 else y,
                    z - nr[2] if z >= hnr2 else z
                )

    G = np.dot(g, b) # compute the G vector
    return G


def compute_G_ver5(b, nr):
    """
    A version of compute_G that calculates coefficients on three 1-D arrays
    and after uses only vector operations.
    """
    nr = np.array(nr)
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b) # compute the G vector
    return G


if __name__ == "__main__":
    pw_output = xmlschema.to_dict("./example5/Ni.xml", "../postqe/schemas/qes.xsd")['output']
    nr = np.array([
        pw_output["basis_set"]["fft_grid"]["@nr1"],
        pw_output["basis_set"]["fft_grid"]["@nr2"],
        pw_output["basis_set"]["fft_grid"]["@nr3"]
    ])
    b1 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b3"])
    b = np.array([b1, b2, b3])

    charge, chargediff = read_charge_file_hdf5("./example5/tmp/Ni.save/charge-density.hdf5", nr)

    # Test conformance
    G = compute_G(b, charge.shape)

    if not np.array_equal(G, compute_G_ver2(b, charge.shape)):
        raise ValueError("compute_G_ver2() result is different!")
    if not np.array_equal(G, compute_G_ver3(b, charge.shape)):
        raise ValueError("compute_G_ver3() result is different!")
    if not np.array_equal(G, compute_G_ver4(b, charge.shape)):
        raise ValueError("compute_G_ver4() result is different!")
    if not np.array_equal(G, compute_G_ver5(b, charge.shape)):
        raise ValueError("compute_G_ver5() result is different!")

    print("##### 'compute_G' versions performance timing #####\n")

    setup = ("from __main__ import b, charge, compute_G, compute_G_ver2, compute_G_ver3, "
             "compute_G_ver4, compute_G_ver5")

    print("compute_G:",
          timeit.repeat('compute_G(b, charge.shape)', setup=setup, number=1, repeat=3))
    print("compute_G_ver2:",
          timeit.repeat('compute_G_ver2(b, charge.shape)', setup=setup, number=1, repeat=3))
    print("compute_G_ver3:",
          timeit.repeat('compute_G_ver3(b, charge.shape)', setup=setup, number=1, repeat=3))
    print("compute_G_ver4:",
          timeit.repeat('compute_G_ver4(b, charge.shape)', setup=setup, number=1, repeat=3))
    print("compute_G_ver5:",
          timeit.repeat('compute_G_ver5(b, charge.shape)', setup=setup, number=1, repeat=3))
