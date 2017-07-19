#!/usr/bin/env python3
#encoding: UTF-8


import numpy as np
from pp import get_from_xml
from readutils import ( read_charge_file_hdf5,
    read_wavefunction_file_hdf5, write_charge, create_header
)

print ('Hello')

from compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from pyqe import pyqe_getcelldms

# get some needed values from the xml output
ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, \
nat, ntyp, nspin, noncolin, pseudodir, nr, nr_smooth = get_from_xml("../tests/Si/Si.xml")
celldms = pyqe_getcelldms(alat, a[0], a[1], a[2], ibrav)

pseudodir = "../tests/Si/"
charge_file = "../tests/Si/charge-density.hdf5"

charge, chargediff = read_charge_file_hdf5(charge_file, nr)
header = create_header("Si", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)

# TESTS

# plot_num = 0
write_charge('../tests/Si/postqeout0', charge, header)

# plot_num = 6
write_charge('../tests/Si/postqeout6', chargediff, header)

# plot_num = 2
v_bare = compute_v_bare(
    ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
)
write_charge('../tests/Si/postqeout2', v_bare, header)


# plot_num = 11
v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,
                        atomic_species, pseudodir)
v_h = compute_v_h(charge, ecutrho, alat, b)
v_tot = v_bare + v_h
write_charge('../tests/Si/postqeout11', v_tot, header)

# plot_num = 1
v_bare = compute_v_bare(
    ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
)
v_h = compute_v_h(charge, ecutrho, alat, b)
charge_core = np.zeros(charge.shape)
v_xc = compute_v_xc(charge, charge_core, str(functional))
v_tot = v_bare + v_h + v_xc
write_charge('../tests/Si/postqeout1', v_tot, header)



#  Plotting tests, second step pp.x

from plot import plotcharge1D, plotcharge2D
from readutils import read_postqe_output_file, write_charge, create_header
from compute_vs import compute_G

charge = read_postqe_output_file('../tests/Si/postqeout0')

# Plot a 1D section
x0 = [0., 0., 0.]
e1 = [1., 0., 0.]
nx = 50
G = compute_G(b, charge.shape)
fig = plotcharge1D(charge, G, a, x0, e1, nx)
fig.show()

# Plot a 2D section
x0 = [0., 0., 0.]
e1 = [1., 0., 0.]
nx = 50
e2 = [0., 1., 0.]
ny = 50
G = compute_G(b, charge.shape)
fig = plotcharge2D(charge, G, a, x0, e1, e2, nx, ny)
fig.show()

import time
time.sleep(1000)