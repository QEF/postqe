#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
import os.path

def run(pars):

    #if (pars. == 0):  # Read the charge and write it in filplot
    #else:
    #    print("Not implemented yet")
    from .api import compute_eos
    compute_eos(pars.prefix, pars.outdir, pars.eos_type, pars.fileout, pars.fileplot, pars.show)


