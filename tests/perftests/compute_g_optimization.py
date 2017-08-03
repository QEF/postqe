#!/usr/bin/env python3
#encoding: UTF-8
import timeit
import numpy as np
import xmlschema
try:
    from cythonfn import compute_G_Cython, compute_G_squared_Cython, compute_Gs_Cython
except ImportError:
    compute_G_Cython = compute_G_squared_Cython = compute_Gs_Cython = None

pi = np.pi


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
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector
    return G


def compute_G_squared(b, nr, ecutrho, alat):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing G^2 at each point, G is the
    corresponding reciprocal lattice vector. Also apply a proper cut-off
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16

    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    G2 = np.zeros((nr[0], nr[1], nr[2]))
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

                G = g0 * b[0] + g1 * b[1] + g2 * b[2]  # compute the G vector
                G2[x, y, z] = np.linalg.norm(G) ** 2  # compute the G^2
                if (G2[x, y, z] > ecutm) or (G2[x, y, z] == 0.0):
                    G2[x, y, z] = bignum  # dummy high value so that n(G)/G^2 is 0

    return G2


def compute_G_squared_ver2(b, nr, ecutrho, alat):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing G^2 at each point, G is the
    corresponding reciprocal lattice vector. Also apply a proper cut-off
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16
    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector

    def g_square(x, y, z):
        g2 = np.linalg.norm(G, axis=3) ** 2  # compute the G^2
        g2[(g2 == 0.0) | (g2 > ecutm)] = bignum  # dummy high value so that n(G)/G^2 is 0
        return g2

    G2 = np.fromfunction(function=g_square, shape=nr, dtype=int)
    return G2


def compute_Gs(b, nr, ecutrho, alat):
    """
    This function computes both a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr and the G^2 moduli. G are the vectors in the
    reciprocal space. Also apply a proper cut-off for G^2
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16

    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    G = np.zeros((nr[0], nr[1], nr[2], 3))
    G2 = np.zeros((nr[0], nr[1], nr[2]))
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
                G2[x, y, z] = np.linalg.norm(G[x, y, z, :]) ** 2  # compute the G^2
                if (G2[x, y, z] > ecutm) or (G2[x, y, z] == 0.0):
                    G2[x, y, z] = bignum  # dummy high value so that n(G)/G^2 is 0

    return G, G2


def compute_Gs_ver2(b, nr, ecutrho, alat):
    """
    This function computes both a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr and the G^2 moduli. G are the vectors in the
    reciprocal space. Also apply a proper cut-off for G^2
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16
    ecutm = 2.0 * ecutrho / ((2.0 * pi / alat) ** 2)  # spherical cut-off for g vectors
    gx = np.fromiter((x - nr[0] if x >= nr[0] // 2 else x for x in range(nr[0])), dtype=int)
    gy = np.fromiter((y - nr[1] if y >= nr[1] // 2 else y for y in range(nr[1])), dtype=int)
    gz = np.fromiter((z - nr[2] if z >= nr[2] // 2 else z for z in range(nr[2])), dtype=int)
    g = np.fromfunction(function=lambda x, y, z: np.array((gx[x], gy[y], gz[z])), shape=nr, dtype=int)
    G = np.dot(np.rollaxis(g, 0, 4), b)  # compute the G vector

    def g_square(x, y, z):
        g2 = np.linalg.norm(G, axis=3) ** 2  # compute the G^2
        g2[(g2 == 0.0) | (g2 > ecutm)] = bignum  # dummy high value so that n(G)/G^2 is 0
        return g2

    G2 = np.fromfunction(function=g_square, shape=nr, dtype=int)
    return G, G2



if __name__ == "__main__":
    pw_output = xmlschema.to_dict("../../examples/example5/Ni.xml", "../../postqe/schemas/qes.xsd")['output']
    nr = np.array([
        pw_output["basis_set"]["fft_grid"]["@nr1"],
        pw_output["basis_set"]["fft_grid"]["@nr2"],
        pw_output["basis_set"]["fft_grid"]["@nr3"]
    ])
    b1 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(pw_output["basis_set"]["reciprocal_lattice"]["b3"])
    b = np.array([b1, b2, b3])

    ecutrho = (pw_output["basis_set"]["ecutrho"])
    alat = (pw_output["atomic_structure"]["@alat"])

    charge, chargediff = read_charge_file_hdf5("../../examples/example5/tmp/Ni.save/charge-density.hdf5", nr)

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
    if compute_G_Cython is not None and not np.array_equal(G, compute_G_Cython(b, charge.shape)):
        raise ValueError("compute_G_Cython() result is different!")

    G2 = compute_G_squared(b, charge.shape, ecutrho, alat)
    if not np.array_equal(G2, compute_G_squared_ver2(b, charge.shape, ecutrho, alat)):
        raise ValueError("compute_G_squared_ver2() result is different!")

    G, G2 = compute_Gs(b, charge.shape, ecutrho, alat)
    _G, _G2 = compute_Gs_ver2(b, charge.shape, ecutrho, alat)
    if not np.array_equal(G, _G) or not np.array_equal(G2, _G2) :
        raise ValueError("compute_Gs_ver2() result is different!")

    # Running performance tests
    print("##### 'compute_G' versions performance timing #####\n")

    setup = ("from __main__ import b, charge, compute_G, compute_G_ver2, compute_G_ver3, "
             "compute_G_ver4, compute_G_ver5, compute_G_Cython")

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
    if compute_G_Cython is not None:
        print("compute_G_Cython:",
              timeit.repeat('compute_G_Cython(b, charge.shape)', setup=setup, number=1, repeat=3))


    print("\n##### 'compute_G_squared' versions performance timing #####\n")

    setup = ("from __main__ import b, charge, ecutrho, alat, compute_G_squared, "
             "compute_G_squared_ver2, compute_G_squared_Cython")

    print("compute_G_squared:",
          timeit.repeat(
              'compute_G_squared(b, charge.shape, ecutrho, alat)',
              setup=setup, number=1, repeat=3))
    print("compute_G_squared_ver2:",
          timeit.repeat(
              'compute_G_squared_ver2(b, charge.shape, ecutrho, alat)',
              setup=setup, number=1, repeat=3))
    if compute_G_squared_Cython is not None:
        print("compute_G_squared_Cython:",
              timeit.repeat(
                  'compute_G_squared_Cython(b, charge.shape, ecutrho, alat)',
                  setup=setup, number=1, repeat=3))

    print("\n##### 'compute_Gs' versions performance timing #####\n")

    setup = ("from __main__ import b, charge, ecutrho, alat, compute_Gs, compute_Gs_ver2, compute_Gs_Cython")

    print("compute_Gs:",
          timeit.repeat('compute_Gs(b, charge.shape, ecutrho, alat)', setup=setup, number=1, repeat=3))
    print("compute_Gs_ver2:",
          timeit.repeat('compute_Gs_ver2(b, charge.shape, ecutrho, alat)', setup=setup, number=1, repeat=3))
    if compute_G_squared_Cython is not None:
        print("compute_Gs_Cython:",
              timeit.repeat('compute_Gs_Cython(b, charge.shape, ecutrho, alat)', setup=setup, number=1, repeat=3))
