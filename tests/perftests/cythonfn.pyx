
import numpy as np
cimport numpy as np

cdef double pi = np.pi

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

                G2[x,y,z] = G[x,y,z,0]*G[x,y,z,0] + G[x,y,z,1]*G[x,y,z,1] + G[x,y,z,2]*G[x,y,z,2]
                if (G2[x,y,z] > ecutm) or (G2[x,y,z] ==0.0):
                    G2[x,y,z] = bignum     # dummy high value so that n(G)/G^2 is 0

    return G, G2

