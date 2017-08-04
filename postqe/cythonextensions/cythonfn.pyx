import numpy as np
cimport numpy as np

cdef double pi = np.pi

# Some C functions and types for complex numbers
cdef extern from "complex.h":
    double cos(double x) nogil
    double sin(double x) nogil
    double complex I
    double complex cexp(double complex)

# return a complex number given rho and phi
cdef double complex complesso(double rho, double phi):
    cdef double a = cos(rho)
    cdef double b = sin(phi)
    return a + b * I

# return a complex number (x + y * j) given the doubles x and y
cdef double complex complesso2(double x, double y):
    return x + y * I


def compute_G_Cython(b,nr):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr. G are the vectors in the reciprocal lattice vector.
    b[0], b[1], b[2] are the reciprocal cell base vectors
    """

    cdef unsigned int x, y, z 
    cdef int g0, g1, g2
    cdef long nr0, nr1, nr2
    cdef double b00, b01, b02
    cdef double b10, b11, b12
    cdef double b20, b21, b22
    cdef double[:,:,:,:] G = np.zeros((nr[0],nr[1],nr[2],3))

    nr0, nr1, nr2 = nr[0], nr[1], nr[2]
    b00, b01, b02 = b[0][0], b[0][1], b[0][2]
    b10, b11, b12 = b[1][0], b[1][1], b[1][2]
    b20, b21, b22 = b[2][0], b[2][1], b[2][2]
    for x in range(0,nr0):
        if (x>=nr0//2):
            g0 = (x-nr0)
        else:
            g0 = x
        for y in range(0,nr1):
            if (y>=nr1//2):
                g1 = (y-nr1)
            else:
                g1 = y
            for z in range(0,nr2):
                if (z>=nr2//2):
                    g2 = (z-nr2)
                else:
                    g2 = z

                G[x,y,z,0] = (g0*b00+g1*b10+g2*b20)    # compute the G vector in x,y,z
                G[x,y,z,1] = (g0*b01+g1*b11+g2*b21)    
                G[x,y,z,2] = (g0*b02+g1*b12+g2*b22)  
                
    return G


def compute_G_squared_Cython(b,nr,ecutrho,alat):
    """
    This function computes a matrix nr[0]*nr[1]*nr[2] containing G^2 at each point, G is the 
    corresponding reciprocal lattice vector. Also apply a proper cut-off 
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """

    cdef double bignum = 1.0E16
    
    cdef unsigned int x, y, z 
    cdef int g0, g1, g2
    cdef unsigned int nr0, nr1, nr2
    cdef double b00, b01, b02
    cdef double b10, b11, b12
    cdef double b20, b21, b22
    cdef double GG0, GG1, GG2

    cdef double ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2)  # spherical cut-off for g vectors   
    cdef double ecutm2 =  ecutm * ecutm 
    cdef double[:,:,:] G2 = np.zeros((nr[0],nr[1],nr[2]))

    nr0, nr1, nr2 = nr[0], nr[1], nr[2]
    b00, b01, b02 = b[0][0], b[0][1], b[0][2]
    b10, b11, b12 = b[1][0], b[1][1], b[1][2]
    b20, b21, b22 = b[2][0], b[2][1], b[2][2]

    for x in range(0,nr0):
        if x>=nr0//2:
            g0 = x-nr0
        else:
            g0 = x
        for y in range(0,nr1):
            if y>=nr1//2:
                g1 = y-nr1
            else:
                g1 = y
            for z in range(0,nr2):
                if z>=nr2//2:
                    g2 = z-nr2
                else:
                    g2 = z
 
                GG0 = (g0*b00+g1*b10+g2*b20)    # compute the G vector in x,y,z
                GG1 = (g0*b01+g1*b11+g2*b21)    
                GG2 = (g0*b02+g1*b12+g2*b22)  
                G2[x,y,z] = GG0*GG0 + GG1*GG1 + GG2*GG2
                if (G2[x,y,z] > ecutm2) or (G2[x,y,z] ==0.0):
                    G2[x,y,z] = bignum     # dummy high value so that n(G)/G^2 is 0
 
    return G2


def compute_Gs_Cython(b,nr,ecutrho,alat):
    """
    This function computes both a matrix nr[0]*nr[1]*nr[2] containing the G vectors at each point
    of the mesh points defined by nr and the G^2 moduli. G are the vectors in the
    reciprocal space. Also apply a proper cut-off for G^2
    ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2). For G^2>ecutm G^2, G^2 should be 0.
    Here G^2 is set to a big number so that n(G)/G^2 is 0 in the inverse FFT.
    """
    bignum = 1.0E16
    
    cdef unsigned int x, y, z 
    cdef int g0, g1, g2
    cdef unsigned int nr0, nr1, nr2
    cdef double b00, b01, b02
    cdef double b10, b11, b12
    cdef double b20, b21, b22
    cdef double GG0, GG1, GG2

    cdef double ecutm = 2.0 * ecutrho / ((2.0*pi/alat)**2)  # spherical cut-off for g vectors  
    cdef double[:,:,:,:] G = np.zeros((nr[0],nr[1],nr[2],3)) 
    cdef double[:,:,:] G2 = np.zeros((nr[0],nr[1],nr[2]))

    nr0, nr1, nr2 = nr[0], nr[1], nr[2]
    b00, b01, b02 = b[0][0], b[0][1], b[0][2]
    b10, b11, b12 = b[1][0], b[1][1], b[1][2]
    b20, b21, b22 = b[2][0], b[2][1], b[2][2]

    for x in range(0,nr0):
        if x>=nr0//2:
            g0 = x-nr0
        else:
            g0 = x
        for y in range(0,nr1):
            if y>=nr1//2:
                g1 = y-nr1
            else:
                g1 = y
            for z in range(0,nr2):
                if z>=nr2//2:
                    g2 = z-nr2
                else:
                    g2 = z
 
                G[x,y,z,0] = (g0*b00+g1*b10+g2*b20)    # compute the G vector in x,y,z
                G[x,y,z,1] = (g0*b01+g1*b11+g2*b21)    
                G[x,y,z,2] = (g0*b02+g1*b12+g2*b22)  

                G2[x,y,z] = G[x,y,z,0]*G[x,y,z,0] + G[x,y,z,1]*G[x,y,z,1] + G[x,y,z,2]*G[x,y,z,2]
                if (G2[x,y,z] > ecutm) or (G2[x,y,z] ==0.0):
                    G2[x,y,z] = bignum     # dummy high value so that n(G)/G^2 is 0
 
    return G, G2


def compute_v_h_Cython(charge,ecutrho,alat,b):
    """
    This function computes the exchange-correlation potential from the charge and
    the type of functional given in input. The charge is a numpy matrix nr1*nr2*nr3.
    The functional is a string identifying the functional as in QE convention.
    
    """    
    # First compute the FFT of the charge          
    fft_charge = np.fft.fftn(charge)
    nr = charge.shape
        
    # Compute G^2 values for the mesh, applying the proper cutoff from ecutrho
    # and alat
    G2 = compute_G_squared_Cython(b,nr,ecutrho,alat)
    
    conv_fact = 2.0 / pi * alat**2
    v = np.fft.ifftn(fft_charge / G2) * conv_fact

    return v.real


def FFTinterp1D_Cython(np.ndarray[double, ndim=3] charge, np.ndarray[double, ndim=4] G, np.ndarray[double, ndim=2] a,
        x0, e1, nx):

    cdef unsigned int i, x, y, z
    cdef long nr0, nr1, nr2
    cdef double x00, x01, x02, e10, e11, e12, m1, deltax, xi, yi, zi, arg, ccos, ssin

    cdef double[:] X = np.zeros(nx)
    cdef np.ndarray[complex] Y = np.zeros(nx, dtype=complex)
    #cdef complex temp, Ytemp

    # normalize e1
    m1 = np.linalg.norm(e1)
    if abs(m1) < 1.0E-6:  # if the module is less than 1.0E-6
        e1 = a[1]
        m1 = np.linalg.norm(e1)
    e10 = e1[0] / m1
    e11 = e1[1] / m1
    e12 = e1[2] / m1

    # Computes the FFT of the charge
    cdef np.ndarray[complex, ndim = 3] fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    nr0, nr1, nr2 = nr[0], nr[1], nr[2]
    x00, x01, x02 = x0[0], x0[1], x0[2]

    # Steps along the e1 direction...
    deltax = m1 / (nx - 1)

    for i in range(0, nx):
        xi = x00 + i * deltax * e10
        yi = x01 + i * deltax * e11
        zi = x02 + i * deltax * e12

        # For each point, evaluate the charge by Fourier interpolation
        for x in range(0, nr0):
            for y in range(0, nr1):
                for z in range(0, nr2):
                    arg = 2.0 * pi * (xi * G[x, y, z, 0] + yi * G[x, y, z, 1] + zi * G[x, y, z, 2])
                    Y[i] = Y[i] + fft_charge[x, y, z] * complesso(arg, arg)

        X[i] = i * deltax
        Y[i] = Y[i] / (nr0 * nr1 * nr2)

    return np.asarray(X), np.asarray(Y)


def FFTinterp2D_Cython(np.ndarray[double, ndim=3] charge, np.ndarray[double, ndim=4] G, np.ndarray[double, ndim=2] a,
        x0, e1, e2, nx, ny):

    cdef unsigned int i, j, x, y, z
    cdef long nr0, nr1, nr2, nxx, nyy
    cdef double x00, x01, x02, e10, e11, e12, e20, e21, e22, m1, m2, deltax, deltay, xi, yi, zi
    cdef double complex arg

    cdef double[:,:] X = np.zeros((nx, ny))
    cdef double[:,:] Y = np.zeros((nx, ny))
    cdef np.ndarray[complex, ndim = 2] Z = np.zeros((nx, ny), dtype=complex)

    cdef np.ndarray[complex] eigx = np.zeros(nx, dtype=complex)
    cdef np.ndarray[complex] eigy = np.zeros(nx, dtype=complex)

    # normalize e1
    m1 = np.linalg.norm(e1)
    if (abs(m1) < 1.0E-6):  # if the module is less than 1.0E-6
        e1 = a[1]
        m1 = np.linalg.norm(e1)
    e10 = e1[0] / m1
    e11 = e1[1] / m1
    e12 = e1[2] / m1

    # normalize e2
    m2 = np.linalg.norm(e2)
    if abs(m2) < 1.0E-6:  # if the module is less than 1.0E-6
        e2 = a[2]
        m2 = np.linalg.norm(e2)
    e20 = e2[0] / m2
    e21 = e2[1] / m2
    e22 = e2[2] / m2

    # Computes the FFT of the charge
    cdef np.ndarray[complex, ndim = 3] fft_charge = np.fft.fftn(charge)
    nr = charge.shape

    nr0, nr1, nr2 = nr[0], nr[1], nr[2]
    x00, x01, x02 = x0[0], x0[1], x0[2]
    nxx, nyy = nx, ny

    # Steps along the e1 and e2 directions...
    deltax = m1 / (nx - 1)
    deltay = m2 / (ny - 1)

    for i in range(0, nxx):
        for j in range(0, nyy):
            X[i, j] = i * deltax
            Y[i, j] = j * deltay

    # loop(s) over the G points
    for x in range(0, nr0):
        for y in range(0, nr1):
            for z in range(0, nr2):

                # eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
                # compute these factors to save CPU time
                for i in range(0, nxx):
                    arg = (2.0 * pi * complesso2(0.0, 1.0) *
                            (i * deltax * (e10 * G[x, y, z, 0] + e11 * G[x, y, z, 1] +  e12 * G[x, y, z, 2]) +
                            (x00 * G[x, y, z, 0] + x01 * G[x, y, z, 1] + x02 * G[x, y, z, 2])))
                    eigx[i] = cexp(arg)     # much more efficient with cexp

                for j in range(0, nyy):
                    arg = (2.0 * pi * complesso2(0.0, 1.0) *
                            (j * deltax * (e20 * G[x, y, z, 0] + e21 * G[x, y, z, 1] + e22 * G[x, y, z, 2])))
                    eigy[j] = cexp(arg)     # much more efficient with cexp

                for i in range(0, nxx):
                    for j in range(0, nyy):
                        Z[i, j] = Z[i, j] + fft_charge[x, y, z] * eigx[i] * eigy[j]

    Z = Z / (nr[0] * nr[1] * nr[2])

    return np.asarray(X), np.asarray(Y), np.asarray(Z)
