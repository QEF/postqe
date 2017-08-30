#!/usr/bin/env python
#encoding: UTF-8
#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
import time, sys
import numpy as np
from readutils import read_charge_file_iotk

from compute_vs import compute_G, compute_G_squared, compute_Gs, compute_v_h
from cythonfun import compute_G_Cython, compute_Gs_Cython, compute_G_squared_Cython, compute_v_h_Cython


if __name__ == "__main__":
    
    ecutrho=1000
    alat=2.052600000000000E+001
    b1 = np.array([-0.1000000E+01,      -0.1000000E+01,       0.1000000E+01])
    b2 = np.array([0.1000000E+01,       0.1000000E+01,       0.1000000E+01])
    b3 = np.array([-0.1000000E+01,       0.1000000E+01,      -0.1000000E+01])
    nr = np.array([75,75,75])
    b = np.array([b1,b2,b3])
    
    start_time = time.time()
    G1 = compute_G(b,nr)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_G... Elapsed time: "+str(elapsed_time)+" s.")
    
    start_time = time.time()
    G2 = compute_G_Cython(b,nr)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_G_Cython... Elapsed time: "+str(elapsed_time)+" s.")
    
    if (G1==G2).all():
        print ("Ok, got the same result!")
    
    start_time = time.time()
    G3 = compute_G_squared(b,nr,ecutrho,alat)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gsquared... Elapsed time: "+str(elapsed_time)+" s.")
    
    start_time = time.time()
    G4 = compute_G_squared_Cython(b,nr,ecutrho,alat)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gsquared_Cython... Elapsed time: "+str(elapsed_time)+" s.")
    
    if (np.sum(G3-G4)<1E-3):
        print ("Ok, got the same result!")
        
    
    start_time = time.time()
    G1, G1bis = compute_Gs(b,nr,ecutrho,alat)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gs... Elapsed time: "+str(elapsed_time)+" s.")
    
    start_time = time.time()
    G2, G2bis = compute_Gs_Cython(b,nr,ecutrho,alat)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gs_Cython... Elapsed time: "+str(elapsed_time)+" s.")
    
    if (G1==G2).all() and (np.sum(G3-G4)<1E-3):
        print ("Ok, got the same result!")
        
    exit()
    charge_file = "charge-density.dat"
    charge = read_charge_file_iotk(charge_file)
    
    start_time = time.time()
    v_h =  compute_v_h(charge,ecutrho,alat,b)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gsquared... Elapsed time: "+str(elapsed_time)+" s.")
    
    start_time = time.time()
    v_h =  compute_v_h_Cython(charge,ecutrho,alat,b)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Testing function compute_Gsquared_Cython... Elapsed time: "+str(elapsed_time)+" s.")