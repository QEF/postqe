#!/usr/bin/env python

import sys
sys.exit(0)

# TODO: refactor as test with latest API
import numpy as np
from postqe.pp import get_from_xml
from postqe.readutils import read_charge_file_hdf5, write_charge, create_header
from postqe.compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from postqe.pyqe import pyqe_getcelldms
from postqe.xmldata import PWData

print ('Hello')


# get some needed values from the xml output
data = PWData("./Ni_pz_nc/Ni.xml", '../postqe/schemas/qes.xsd')
celldms = pyqe_getcelldms(data.alat, data.a[0], data.a[1], data.a[2], data.ibrav)

pseudodir = "./Ni_pz_nc/"
charge_file = "./Ni_pz_nc/charge-density.hdf5"

charge, chargediff = read_charge_file_hdf5(charge_file, data.nr)
header = create_header("Ni", data.nr, data.nr_smooth, data.ibrav, celldms, data.nat,
                       data.ntyp, data.atomic_species, data.atomic_positions)

# TESTS

# plot_num = 0
write_charge('./Ni_pz_nc/postqeout0', charge, header)

# plot_num = 6
write_charge('./Ni_pz_nc/postqeout6', chargediff, header)

# plot_num = 2
v_bare = compute_v_bare(
    data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
    data.atomic_positions, data.atomic_species, pseudodir
)
write_charge('./Ni_pz_nc/postqeout2', v_bare, header)


# plot_num = 11
v_bare = compute_v_bare(data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
                        data.atomic_positions, data.atomic_species, pseudodir)
v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
v_tot = v_bare + v_h
write_charge('./Ni_pz_nc/postqeout11', v_tot, header)

# plot_num = 1
v_bare = compute_v_bare(
    data.ecutrho, data.alat, data.a[0], data.a[1], data.a[2], data.nr,
    data.atomic_positions, data.atomic_species, pseudodir
)
v_h = compute_v_h(charge, data.ecutrho, data.alat, data.b)
charge_core = np.zeros(charge.shape)
v_xc = compute_v_xc(charge, charge_core, str(data.functional))
v_tot = v_bare + v_h + v_xc
write_charge('./Ni_pz_nc/postqeout1', v_tot, header)



#  Plotting tests, second step pp.x
from postqe.plot import plotcharge1D, plotcharge2D
from postqe.readutils import read_postqe_output_file
from postqe.compute_vs import compute_G

charge = read_postqe_output_file('./Ni_pz_nc/postqeout0')

# Plot a 1D section
x0 = [0., 0., 0.]
e1 = [1., 0., 0.]
nx = 50
G = compute_G(data.b, charge.shape)
fig = plotcharge1D(charge, G, data.a, x0, e1, nx)
fig.show()

# Plot a 2D section
x0 = [0., 0., 0.]
e1 = [1., 0., 0.]
nx = 50
e2 = [0., 1., 0.]
ny = 50
G = compute_G(data.b, charge.shape)
fig = plotcharge2D(charge, G, data.a, x0, e1, e2, nx, ny)
fig.show()

import time
time.sleep(1000)