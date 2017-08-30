#encoding: UTF-8

"""
Functions for reading the charge in the old binary format. Tricky and not supported.

Warning: Reading binary files is tricky and platform dependent and may not work properly.
This is kept only for the old QE binary format, but we strongly recommend to use
the new hdf5 format.
"""

import numpy as np
import re
from struct import *


def read_line(f):
    """
    Read a line from a binary file. End character must be b'\n'
    """
    line = b''
    byte = f.read(1)
    while byte != b'\n':
        line += byte
        byte = f.read(1)
            
    return line


def read_n_real_numbers(f,nr1,nr2,kind=8):
    """
    Read n real numbers from a binary (fortran) file section. Extra characters
    are present in the file and discarded. Potentially troublesome...
    """
    a = np.zeros((nr1,nr2))
    byte = f.read(12)    # discard the first 12 bytes (extra info written by iotk/fortran)
    for j in range(0,nr2):
        for i in range(0,nr1):
            byte = f.read(kind)    
            x = unpack('d',byte)  
            a[i,j] = x[0]     # x is a list, take the first element only
        
    return a


def get_info(info_line):
    """
    Get nr1, nr2, nr3 from line <INFO...>
    """
    pat = re.compile(r'nr\d="\d+"')    
    info_line_nrx = pat.findall(info_line)
    
    nr = np.zeros(3,dtype=int)
    for i in range(0,len(info_line_nrx)):
        nr[i] = info_line_nrx[i].split("\"")[1]
    
    return nr[0], nr[1], nr[2]


def read_charge_file_iotk(filename):
    """
    Read a old binary charge file written with QE and iotk.
    Warning: the binary format is tricky and platform dependent. It will not be supported in the future.
    Use the hdf5 format instead.
    """
 
    tempcharge = []
    with open(filename, "rb") as f:
        line = b''
        while (line!=b'  <CHARGE-DENSITY>'):    # skip lines till the <CHARGE-DENSITY> tag is found
            line = read_line(f)
            
        discard = read_line(f)  # discard binary line here
        line = read_line(f)     # this is the line <INFO nr1=.. />
        info_line = line.decode('utf-8')
        nr1, nr2, nr3 = get_info(info_line)

        # data are grouped in nr1*nr2 sequential reals in the charge file, each
        # starting with <z.1, <z.2, etc. tag and ending with a corresponding
        # </z.1>, </z.2> containing extra information...
        charge = np.zeros((nr1, nr2, nr3))
        for i in range(0, nr3):
            discard = read_line(f)
            line = read_line(f)     # this is the line <z.1 type="real" size="2025" kind="8"> or similar
        
            temp = read_n_real_numbers(f, nr1, nr2)
            charge[:,:,i] = temp
            
            discard = read_line(f)
            line = read_line(f)     # this is end tag </z.1>
        
    return charge

