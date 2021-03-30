#!/usr/bin/env python3
#encoding: UTF-8
import os
import sys

# Adds the the package directory to sys.path, in order to make
# the development module loadable also without set PYTHONPATH.
pkg_search_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if sys.path[0] != pkg_search_path:
    sys.path.insert(0, pkg_search_path)

import numpy as np
from postqe.readutils import read_charge_file_hdf5, write_charge, create_header
from postqe.compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from postqe.pyqe import pyqe_getcelldms
from postqe.xmldata import PWData

print ('Hello')

data = PWData("./Ni_pbe_us/Ni.xml", '../postqe/schemas/qes.xsd')
celldms = pyqe_getcelldms(data.alat, data.a[0], data.a[1], data.a[2], data.ibrav)

pseudodir = "./Ni_pbe_us/"
charge_file = "./Ni_pbe_us/charge-density.hdf5"

charge, chargediff = read_charge_file_hdf5(charge_file, data.nr)
header = create_header("Ni", data.nr, data.nr_smooth, data.ibrav, celldms, data.nat,
                       data.ntyp, data.atomic_species, data.atomic_positions)

# TESTS

# plot_num = 0
write_charge('./Ni_pbe_us/postqeout0', charge, header)

# plot_num = 6
write_charge('./Ni_pbe_us/postqeout6', chargediff, header)

# plot_num = 2
v_bare = compute_v_bare(
    data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
    data.atomic_positions, data.atomic_species, pseudodir
)
write_charge('./Ni_pbe_us/postqeout2', v_bare, header)


# plot_num = 11
v_bare = compute_v_bare(data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
                        data.atomic_positions, data.atomic_species, pseudodir)
v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
v_tot = v_bare + v_h
write_charge('./Ni_pbe_us/postqeout11', v_tot, header)

# plot_num = 1, not working yet
v_bare = compute_v_bare(
    data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
    data.atomic_positions, data.atomic_species, pseudodir
)
v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
charge_core = np.zeros(charge.shape)
v_xc = compute_v_xc(charge, charge_core, str(data.functional))
v_tot = v_bare + v_h + v_xc
write_charge('./Ni_pbe_us/postqeout1', v_tot, header)


print("END")

#  Plotting tests, second step pp.x

from plot import plotcharge1D, plotcharge2D
from readutils import read_postqe_output_file, write_charge, create_header
from compute_vs import compute_G

charge = read_postqe_output_file('./Ni_pbe_us/postqeout0')

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