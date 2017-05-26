#!/usr/bin/env python3
#encoding: UTF-8

import time, sys
import numpy as np
from readutils import read_line, read_n_real_numbers,\
read_charge_file_iotk, read_charge_file_hdf5, read_wavefunction_file_iotk,\
read_wavefunction_file_hdf5, write_charge, create_header
from compute_vs import compute_v_bare, compute_v_h, compute_v_xc
from pyQ import pyq_getcelldm as calcola_celldm
import settings

def get_from_xml(fname, schema = None):
    """
    Get some useful values from xml file
    """

    import xmlschema
    from xml.etree import ElementTree as ET
    try:
        from urllib.request import urlopen
        from urllib.parse import urlparse
    except ImportError:
        from urllib2 import urlopen
        from urllib2 import urlparse
        from urlparse import urlparse

    root = ET.parse(fname).getroot()
    schemaLoc = root.get('{http://www.w3.org/2001/XMLSchema-instance}schemaLocation').split()[1]
    schemaLoc = urlparse(schemaLoc)
    if schemaLoc.scheme =='':
        xd = xmlschema.XMLSchema(schemaLoc.path)
    else:
        lines = urlopen(schemaLoc.geturl() ).readlines()
        with open('temp.xsd','w') as tempSchemaFile:
            tempSchemaFile.writelines([l.decode("utf-8") for l in lines ])
        xd = xmlschema.XMLSchema('temp.xsd')
        import os
        os.remove('temp.xsd')
    print ("Reading xml file: ",fname)
    d = xd.to_dict(fname)
    try:
        settings.pseudodir = d["input"]["control_variables"]["pseudo_dir"]
    except KeyError:
        pass
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
    
    # for subsequent loops it is important to have always lists for atomic_positions
    # and atomic_species. If this is not, convert
    if (type(a_s)==type([])):
        atomic_species = a_s
    else:
        atomic_species = [a_s,]
    if (type(a_p)==type([])):
        atomic_positions = a_p
    else:
        atomic_positions = [a_p,]
    
    a = np.array([a1,a2,a3])
    b = np.array([b1,b2,b3])
    
    return ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species, nat, ntyp, nspin, noncolin



def get_input_parameters():
    """
    Get postprocessing input parameters using argparse.
    """
    import argparse 
    
    parser = argparse.ArgumentParser(description='QE post processing')
    
    parser.add_argument('-plot_num', type=int, nargs='?', default=0, choices=range(0, 21),
    help="""selects what to save in filplot:\n
    0  = electron (pseudo-)charge density\n
    1  = total potential V_bare + V_H + V_xc\n
    2  = local ionic potential V_bare\n
    3  = local density of states at e_fermi (number of states per volume, in bohr^3,\
         per energy unit, in Ry)\n\
    4  = local density of electronic entropy\n\
    5  = STM images\n\
        Tersoff and Hamann, PRB 31, 805 (1985)\n\
    6  = spin polarization (rho(up)-rho(down))\n
    
    7  = contribution of a selected wavefunction to the
        (pseudo-)charge density. For norm-conserving PPs,
        |psi|^2 (psi=selected wavefunction). Noncollinear case:
        contribution of the given state to the charge or
        to the magnetization along the direction indicated
        by spin_component (0 = charge, 1 = x, 2 = y, 3 = z )

    8  = electron localization function (ELF)

    9  = charge density minus superposition of atomic densities

    10 = integrated local density of states (ILDOS)
        from emin to emax (emin, emax in eV)
        if emax is not specified, emax=E_fermi

    11 = the V_bare + V_H potential

    12 = the sawtooth electric field potential (if present)

    13 = the noncollinear magnetization.

    17 = all-electron valence charge density
        can be performed for PAW calculations only
        requires a very dense real-space grid!

    18 = The exchange and correlation magnetic field in
        the noncollinear case

    19 = Reduced density gradient
        (J. Chem. Theory Comput. 7, 625 (2011))
        Set the isosurface between 0.3 and 0.6 to plot the
        non-covalent interactions (see also plot_num = 20)

    20 = Product of the electron density (charge) and the second
        eigenvalue of the electron-density Hessian matrix;
        used to colorize the RDG plot (plot_num = 19)

    21 = all-electron charge density (valence+core).
        For PAW calculations only; requires a very dense
        real-space grid.
    """

        )
        
    #default_prefix = "Si"
    default_prefix = "SiO2"
    #default_prefix = "SrTiO3"
    parser.add_argument('-prefix', type=str, nargs='?', default=default_prefix,
                    help='prefix of files saved by program pw.x')
    default_outdir = "../tests/"+default_prefix
    parser.add_argument('-outdir', type=str, nargs='?', default=default_outdir,
                    help='directory containing the input data, i.e. the same as in pw.x')
    parser.add_argument('-filplot', type=str, nargs='?', default="filplot",
                    help='file \"filplot\" contains the quantity selected by plot_num\
                        (can be saved for further processing)')

    parser.add_argument('-spin_component', type=int, nargs='?', default=0,
                    help='if plot_num==0: 0 = total charge (default value), 1 = spin up charge, 2 = spin down charge.\
                          if plot_num==1: 0 = spin averaged potential (default value), 1 = spin up potential,\
                          2 = spin down potential.')

    args = parser.parse_args()
    
    return args


if __name__ == "__main__":
    
    start_time = time.time()
    
    # get the input parameters
    pars = get_input_parameters()
    
    print (pars)
    
    # get some needed values from the xml output
    ecutwfc, ecutrho, ibrav, alat, a, b, functional, atomic_positions, atomic_species,\
    nat, ntyp, nspin, noncolin = get_from_xml(pars.outdir+"/"+pars.prefix+".xml",settings.schema)
    celldms = calcola_celldm(alat,a[0],a[1],a[2],ibrav)
      
    charge_file = pars.outdir+"/charge-density.dat"
    charge = read_charge_file_iotk(charge_file)
    nr = charge.shape
    header = create_header(pars.prefix,nr,ibrav,celldms,nat,ntyp,atomic_species,atomic_positions)
        
    if (pars.plot_num==0):   # Read the charge and write it in filplot
        # TO DO: handle different spin cases
        write_charge(pars.filplot,charge,header)
        
    elif (pars.plot_num==1):
        v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,\
        atomic_species, settings.pseudodir)      
        v_h =  compute_v_h(charge,ecutrho,alat,b)
        charge_core = np.zeros(charge.shape)    # only for now, later in input
        v_xc = compute_v_xc(charge,charge_core,str(functional))
        v_tot = v_bare + v_h + v_xc
        write_charge(pars.filplot,v_tot,header)     
            
    elif (pars.plot_num==2):
        v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,\
        atomic_species, pseudodir)
        write_charge(pars.filplot,v_bare,header)

    # TO DO: add plot_num == 6 case, should be easy


    elif (pars.plot_num==11):
        v_bare = compute_v_bare(ecutrho, alat, a[0], a[1], a[2], nr, atomic_positions,\
        atomic_species,settings.pseudodir)    
        v_h =  compute_v_h(charge,ecutrho,alat,b)
        v_tot = v_bare + v_h
        write_charge(pars.filplot,v_tot,header)        
        
    else:
        print ("Not implemented yet")

  
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")
