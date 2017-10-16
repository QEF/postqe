#!/usr/bin/env python3
#encoding: UTF-8

from .api import compute_eos, compute_band_structure

def run(pars):

    # TODO: implement submenus
    #compute_eos(pars.prefix, pars.outdir, pars.eos_type, pars.fileout, pars.fileplot, pars.show)

    compute_band_structure(pars.prefix, pars.outdir, pars.schema, pars.reference_energy,
                           pars.emin, pars.emax, pars.fileplot, pars.show)