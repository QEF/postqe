#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
import os.path
from readutils import (
    read_line, read_n_real_numbers, read_charge_file_iotk, read_charge_file_hdf5,
    read_wavefunction_file_hdf5, write_charge, create_header
)
from compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from qeutils import py_getcelldms


def get_from_xml(filename):
    """
    Get some useful values from xml file
    """
    import xmlschema

    ##########################################################
    # TODO for whatever reason this is not working now
    #schemaLoc = xmlschema.fetch_schema(filename)
    #xs = xmlschema.XMLSchema(schemaLoc)
    #
    # temporary local solution
    xs = xmlschema.XMLSchema('schemas/qes.xsd')
    ##########################################################

    print ("Reading xml file: ", filename)
    d = xs.to_dict(filename)
    try:
        pseudodir = d["input"]["control_variables"]["pseudo_dir"]
    except KeyError:
        pseudodir = './'   # DB: ['./', os.path.dirname(filename)] ... alternatives?
    dout = d["output"]
    ecutwfc = (dout["basis_set"]["ecutwfc"])
    ecutrho = (dout["basis_set"]["ecutrho"])
    alat = (dout["atomic_structure"]["@alat"])
    a1 = np.array(dout["atomic_structure"]["cell"]["a1"])
    a2 = np.array(dout["atomic_structure"]["cell"]["a2"])
    a3 = np.array(dout["atomic_structure"]["cell"]["a3"])   
    ibrav = (dout["atomic_structure"]["@bravais_index"])
    b1 = np.array(dout["basis_set"]["reciprocal_lattice"]["b1"])
    b2 = np.array(dout["basis_set"]["reciprocal_lattice"]["b2"])
    b3 = np.array(dout["basis_set"]["reciprocal_lattice"]["b3"])
    functional = np.array(dout["dft"]["functional"])
    a_p = (dout["atomic_structure"]["atomic_positions"]["atom"])
    a_s = (dout["atomic_species"]["species"])
    nat = (dout["atomic_structure"]["@nat"])
    ntyp = (dout["atomic_species"]["@ntyp"])
    nspin = (dout["magnetization"]["do_magnetization"])
    noncolin = (dout["magnetization"]["noncolin"])
    nr = np.array([dout["basis_set"]["fft_grid"]["@nr1"],dout["basis_set"]["fft_grid"]["@nr2"],dout["basis_set"]["fft_grid"]["@nr3"]])
    nr_smooth = np.array([dout["basis_set"]["fft_smooth"]["@nr1"],dout["basis_set"]["fft_smooth"]["@nr2"],dout["basis_set"]["fft_smooth"]["@nr3"]])

    # for subsequent loops it is important to have always lists for atomic_positions
    # and atomic_species. If this is not, convert
    if (type(a_s)==type([])):
        atomic_species = a_s
    else:
        atomic_species = [a_s]

    if (type(a_p)==type([])):
        atomic_positions = a_p
    else:
        atomic_positions = [a_p]
    
    a = np.array([a1,a2,a3])
    b = np.array([b1,b2,b3])
    
    return (ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions,
            atomic_species, nat, ntyp, nspin, noncolin, pseudodir, nr, nr_smooth)


def run(pars):
    """Run pp with configuration parameters"""
    ### DB: creare un oggetto per i parametri??

    # get some needed values from the xml output
    ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, \
    nat, ntyp, nspin, noncolin, pseudodir, nr, nr_smooth = get_from_xml("Ni.xml")
    celldms = py_getcelldms(alat, a[0], a[1], a[2], ibrav)
      
    charge_file = pars.outdir + "/charge-density.hdf5"
    charge, chargediff  = read_charge_file_hdf5(charge_file)
    header = create_header("Ni", nr, nr_smooth, ibrav, celldms, nat, ntyp, atomic_species, atomic_positions)
        
    if (pars.plot_num==0):   # Read the charge and write it in filplot
        # TODO: handle different spin cases (easy)
        write_charge(pars.filplot, charge, header)
        
    elif (pars.plot_num==1):
        v_bare = compute_v_bare(
            ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
        )
        v_h =  compute_v_h(charge,ecutrho, alat, b)
        charge_core = np.zeros(charge.shape)    # only for now, later in input
        v_xc = compute_v_xc(charge, charge_core, str(functional))
        v_tot = v_bare + v_h + v_xc
        write_charge(pars.filplot,v_tot,header)     
            
    elif (pars.plot_num==2):
        v_bare = compute_v_bare(
            ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions, atomic_species, pseudodir
        )
        write_charge(pars.filplot,v_bare,header)

    if (pars.plot_num==6):   # Write the charge difference (spin up - spin down) for magnetic systems
        write_charge(pars.filplot, chargediff, header)

    elif (pars.plot_num==11):
        v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,
                                atomic_species, pseudodir)
        v_h =  compute_v_h(charge, ecutrho, alat, b)
        v_tot = v_bare + v_h
        write_charge(pars.filplot, v_tot, header)
        
    else:
        print ("Not implemented yet")
