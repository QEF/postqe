#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the the band structure of silicon.
"""
    
if __name__ == "__main__":
    from postqe import get_band_structure

    bs = get_band_structure(prefix='Si', schema='../../schemas/qes.xsd', reference_energy=0)
    fig = bs.plot(emin=-20, emax=50, show=True, filename='Sibands.png')
