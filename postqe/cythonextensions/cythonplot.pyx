import numpy as np
cimport numpy as np
#from libc.math cimport sin
#from libc.math cimport cos

cdef extern from "complex.h":
    double cos(double x) nogil
    double sin(double x) nogil
    float complex I

cdef double pi = np.pi

cdef float complex complesso(float x, float y):
    cdef float a = cos(x)
    cdef float b = sin(y)
    return a + b * I

def FFTinterp1D_Cython(np.ndarray[double, ndim=3] charge, np.ndarray[double, ndim=4] G, np.ndarray[double, ndim=2] a, x0, e1, nx):

    cdef unsigned int i, x, y, z
    cdef long nr0, nr1, nr2
    cdef double x00, x01, x02, e10, e11, e12, m1, deltax, xi, yi, zi, arg, ccos, ssin

    cdef double[:] X = np.zeros(nx)
    cdef np.ndarray[complex] Y = np.zeros(nx, dtype=complex)
    cdef complex temp, Ytemp

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
        Ytemp = complesso(0,0)
        for x in range(0, nr0):
            for y in range(0, nr1):
                for z in range(0, nr2):
                    arg = 2.0 * pi * (xi * G[x, y, z, 0] + yi * G[x, y, z, 1] + zi * G[x, y, z, 2])
                    temp = complesso(arg, arg)
                    Ytemp = Ytemp + fft_charge[x, y, z] * temp
                    #Y[i] = Y[i] + fft_charge[x, y, z] * temp

        X[i] = i * deltax
        Y[i] = Ytemp / (nr0 * nr1 * nr2)
        #Y[i] = Y[i] / (nr0 * nr1 * nr2)
        #print(X[i], Y[i].real)

    return X, Y


def FFTinterp2D_Cython(charge, G, a, x0, e1, e2, nx, ny):
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
                    eigx[i] = np.exp(2.0 * pi * complex(0.0, 1.0) * (i * deltax * \
                                                                     (e1[0] * G[x, y, z, 0] + e1[1] * G[x, y, z, 1] +
                                                                      e1[2] * G[x, y, z, 2]) + \
                                                                     (x0[0] * G[x, y, z, 0] + x0[1] * G[x, y, z, 1] +
                                                                      x0[2] * G[x, y, z, 2])))

                eigy = np.zeros(ny, dtype=complex)
                for j in range(0, ny):
                    eigy[j] = np.exp(2.0 * pi * complex(0.0, 1.0) * (j * deltax * \
                                                                     (e2[0] * G[x, y, z, 0] + e2[1] * G[x, y, z, 1] +
                                                                      e2[2] * G[x, y, z, 2])))

                for i in range(0, nx):
                    for j in range(0, ny):
                        temp[i, j] += fft_charge[x, y, z] * eigx[i] * eigy[j]

    Z = temp.real / (nr[0] * nr[1] * nr[2])

    return X, Y, Z