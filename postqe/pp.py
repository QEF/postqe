#!/usr/bin/env python3
#encoding: UTF-8

from .api import compute_eos, compute_band_structure, comput_dos, compute_charge

def run(pars):
    # just go through the possible commands and call the proper API function
    if (pars.commands == 'eos'):
        compute_eos(pars.prefix, pars.outdir, pars.eos_type, pars.fileout, pars.fileplot, pars.show)
    elif (pars.commands == 'bands'):
        compute_band_structure(pars.prefix, pars.outdir, pars.schema, pars.reference_energy, pars.emin, pars.emax,
                               pars.fileplot, pars.show)
    elif (pars.commands == 'dos'):
        if (pars.emin == None) or (pars.max == None):
            window = None
        else:
            window = (pars.emin, pars.emax)
        comput_dos(pars.prefix, pars.outdir, pars.schema, pars.width, window, pars.npts, pars.fileout, pars.fileplot,
                   pars.show)
    elif (pars.commands == 'charge'):
        compute_charge(pars.prefix, pars.outdir, pars.schema, pars.fileout, pars.x0, pars.e1, pars.nx,
                       pars.e2, pars.ny, pars.e3, pars.nz, pars.radius, pars.dim, pars.ifmagn, pars.plot_file,
                       pars.method, pars.format, pars.show)
    elif (pars.commands == 'potential'):
        raise NotImplementedError
    else:
        print('Command not implemented! Exiting...')