#!/usr/bin/env python3
#encoding: UTF-8

from .api import compute_eos, compute_band_structure, comput_dos

def run(pars):
    # just go through the possible commands and call the proper API function
    if (pars.commands == 'eos'):
        compute_eos(pars.prefix, pars.outdir, pars.eos_type, pars.fileout, pars.fileplot, pars.show)
    elif (pars.commands == 'bands'):
        compute_band_structure(pars.prefix, pars.outdir, pars.schema, pars.reference_energy,
                           pars.emin, pars.emax, pars.fileplot, pars.show)
    elif (pars.commands == 'dos'):
        comput_dos(pars.prefix, pars.outdir, pars.schema, pars.width, (pars.emin, pars.emax), pars.npts,
               pars.fileout,  pars.fileplot, pars.show)
    else:
        print('Command not implemented! Exiting...')