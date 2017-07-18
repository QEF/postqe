#!/usr/bin/env python3
#encoding: UTF-8

#from pyqe import pyqe_xc as xc
#from pyqe import pyqe_getcelldms, pyqe_recips, pyqe_latgen, pyqe_struct_fact
#from pyqe import pyqe_get_gg_list, pyqe_get_gl, pyqe_get_igtongl
#from pyqe2 import pyqe_vloc_of_g

#exit()

import numpy as np
from pp import get_from_xml
from readutils import (
    read_line, read_n_real_numbers, read_charge_file_iotk, read_charge_file_hdf5,
    read_wavefunction_file_hdf5, write_charge, create_header
)

print ('Hello')

from compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from pyqe import pyqe_getcelldms

# get some needed values from the xml output
ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, \
nat, ntyp, nspin, noncolin, pseudodir, nr, nr_smooth = get_from_xml("../tests/Ni_pz_nc/Ni.xml")
celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

print (alat, a[0], a[1], a[2], ibrav)
print (celldms)
print (atomic_positions)
pseudodir = "../tests/Ni_pz_nc/"

charge_file = "../tests/Ni_pz_nc/charge-density.hdf5"
charge, chargediff = read_charge_file_hdf5(charge_file, nr)
header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

# TESTS

# plot_num = 0
write_charge('../tests/Ni_pz_nc/postqeout0', charge, header)

# plot_num = 6
write_charge('../tests/Ni_pz_nc/postqeout6', chargediff, header)

# plot_num = 2
v_bare = compute_v_bare(
    ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
)
write_charge('../tests/Ni_pz_nc/postqeout2', v_bare, header)


# plot_num = 11
v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,
                        atomic_species, pseudodir)
v_h = compute_v_h(charge, ecutrho, alat, b)
v_tot = v_bare + v_h
write_charge('../tests/Ni_pz_nc/postqeout11', v_tot, header)

# plot_num = 1
v_bare = compute_v_bare(
    ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
)
v_h = compute_v_h(charge, ecutrho, alat, b)
charge_core = np.zeros(charge.shape)
v_xc = compute_v_xc(charge, charge_core, str(functional))
v_tot = v_bare + v_h + v_xc
write_charge('../tests/Ni_pz_nc/postqeout1', v_tot, header)

