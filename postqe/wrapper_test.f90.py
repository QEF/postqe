import numpy as np
from postqe import get_charge, get_cell_data, get_calculation_data, \
    compute_v_bare, compute_v_h, compute_v_xc, compute_G, plot1D_FFTinterp

from wtest import pyqe_v_xc

fin = "./Ni.xml"  # file xml produce by QE
ibrav, alat, a, b, nat, ntyp, \
atomic_positions, atomic_species = get_cell_data(fin)  # get some data on the unit cell
prefix, outdir, ecutwfc, ecutrho, functional, lsda, noncolin, pseudodir, nr, nr_smooth = \
    get_calculation_data(fin)  # get some data on the QE calculation

charge, chargediff = get_charge(fin)  # get the charge (and charge diff) from the HDF5 file

nr = charge.shape
print(nr[0])
print (charge[0,0,0],charge[0,0,1],charge[0,1,0],charge[1,0,0],charge[39,39,39])
ex1, ex2, v = pyqe_v_xc(charge, charge, charge, lsda, nr[0], nr[1], nr[2] )
print (v[0,0,0],v[0,0,1],v[0,1,0],v[1,0,0],)
exit()
# Compute the bare potential (v_bare), the Hartree potential (v_h) and the exhange-correlation potential (v_xc)
v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir)
v_h = compute_v_h(charge, ecutrho, alat, b)
charge_core = np.zeros(charge.shape)