import numpy as np
from readutils import read_pseudo_file, read_pseudo_file_old

# f2py module
from pyqe import pyqe_getcelldms, pyqe_recips, pyqe_latgen
from pyqe import pyqe_get_gg_list, pyqe_get_gl, pyqe_get_igtongl
from pyqe2 import pyqe_vloc_of_g, pyqe_struct_fact


def generate_glists(alat, at1, at2, at3, nr1, nr2, nr3, ecutrho):
    nrrr = nr1 * nr2 * nr3
    tpiba = 2 * np.pi / alat
    tpiba2 = tpiba**2
    bg1, bg2, bg3 = pyqe_recips(at1 / alat, at2 / alat, at3 / alat)
    g, gg, mill = pyqe_get_gg_list(nrrr, nr1, nr2, nr3, bg1, bg2, bg3)
    gzip = zip(gg, g, mill)
    gzip = (el for el in gzip if el[0] <= ecutrho / tpiba2)
    gzip = sorted(gzip, key=lambda el: el[0])
    mill_cut = [el[2] for el in gzip]
    g_cut = [el[1] for el in gzip]
    gg_cut =[el[0] for el in gzip]
    igtongl, ngl = pyqe_get_igtongl(gg_cut)
    gl = pyqe_get_gl(gg_cut, ngl)
    return g_cut, gg_cut, mill_cut, igtongl, gl


def vloc_of_g(rab, r, vloc_r, zp, alat, omega, gl):
    """
    :type omega: float
    :type rab: np.ndarray
    :type r: np.ndarray
    :type vloc_r: np.ndarray
    :type zp: float
    :type alat: float
    :type gl: np.ndarray
    """
    msh = len(rab)
    tpiba2 = (np.pi * 2.e0 / alat) ** 2
    return pyqe_vloc_of_g(msh, r=r, rab=rab, vloc_at=vloc_r, zp=zp, tpiba2=tpiba2, gl=gl, omega=omega)


def shift_and_transform(nr1, nr2, nr3, vlocs, strct_facs, mill, igtongl):
    aux = np.zeros([nr1, nr2, nr3], dtype="D")
    for nt in range(len(vlocs)):
        for ig in range(len(mill)):
            ijk = mill[ig]
            i = int(ijk[0])
            j = int(ijk[1])
            k = int(ijk[2])
            aux[i, j, k] += strct_facs[nt][ig ] * vlocs[nt][igtongl[ig]-1]
    return np.fft.ifftn(aux)*(nr1*nr2*nr3)


def compute_struct_fact(tau, alat, g):
    str_fact, checkgg, check_tau = pyqe_struct_fact(tau / alat, g)
    return str_fact


def wrap_setlocal(alat, at1, at2, at3, nr1, nr2, nr3, atomic_positions, species, ecutrho, pseudodir="./"):
    omega = abs(at1[0] * at2[1] * at3[2] + at1[1] * at2[2] * at3[0] + at1[2] * at2[0] * at3[1] -
                at3[0] * at2[1] * at1[2] - at3[1] * at2[2] * at1[0] - at3[2] * at2[0] * at3[1])

    g, gg, mill, igtongl, gl = generate_glists(alat, at1, at2, at3, nr1, nr2, nr3, ecutrho)
    #
    strct_facs = []

    for typ in species:
        tau_spec = []
        name = typ["@name"]
        for pos in atomic_positions:
            if pos["@name"] == name:
                coords = [float(x) for x in pos['$'] ]
                new_atomic_positions = np.array(coords  )*alat
                tau_spec.append(new_atomic_positions)
        tau_spec = np.array(tau_spec)
        str_fact = compute_struct_fact(tau_spec, alat, g)
        strct_facs.append(str_fact)
        #
    vlocs = []
    for typ in species:
        filename = typ["pseudo_file"]
        pseudo = read_pseudo_file(pseudodir+filename)
        vloc_r = pseudo["PP_LOCAL"]
        r = pseudo["PP_MESH"]["PP_R"]
        rab = pseudo["PP_MESH"]["PP_RAB"]
        #
        if ( "PP_HEADER" in pseudo.keys()):
            try:
                header = pseudo["PP_HEADER"].split('\n')
                my_line = [l for l in header if 'Z valence' in l][0]
                zp = float(my_line.split()[0])
            except AttributeError:
                zp = pseudo['PP_HEADER']['z_valence']
        else:
            with open (filename, 'r') as f:
                my_line = [ l for l in f.readlines() if 'z_valence=' in l][0]
                zp = float(my_line.split('"')[1])
        #
        vloc_g = vloc_of_g(rab, r, vloc_r, zp, alat, omega,  gl)

        vlocs.append(vloc_g)
        #
    vltot = shift_and_transform(nr1, nr2, nr3, vlocs, strct_facs, mill, igtongl)
    prova = np.real(vltot)  
    return prova
