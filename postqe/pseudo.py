#!/usr/bin/env python3
#encoding: UTF-8

import numpy as np
from xml.etree import ElementTree as ET

# f2py module
#from .pyqe import pyqe_getcelldms, pyqe_recips, pyqe_latgen
#from .pyqe import pyqe_get_gg_list, pyqe_get_gl, pyqe_get_igtongl
from .setlocal import generate_glists

from .xmlfile import get_cell_data

def iter_upf_file(upffile):
    """
    Creates an iterator over the lines of an UPF file,
    inserting the root <UPF> tag when missing.
    """
    with open(upffile, 'r') as f:
        fake_root = None
        for line in f:
            if fake_root is not None:
                yield line.replace('&input', '&amp;input')
            else:
                line = line.strip()
                if line.startswith("<UPF") and line[4] in ('>', ' '):
                    yield line
                    fake_root = False
                elif line:
                    yield "<UPF>"
                    yield line
                    fake_root = True
    if fake_root is True:
        yield "</UPF>"


class Pseudo:
    """
    A class for a pseudopotential.
    """
    def __init__(self, *args, **kwargs):
        """Create charge object from """
        self.setvars(*args, **kwargs)

    def setvars(self):
        # TODO: check if it makes sense to implement a non-void constructor
        pass

    def set_calculator(self, calculator):
        self.calculator = calculator

    def read(self, filename):
        """
        This function reads a pseudopotential XML-like file in the QE UPF format (text) and stores
        the content of each tag in the class. The file is read in strings and completed with a root
         UPF tag when it lacks, to avoids an XML syntax error.

        :param filename: an UPF pseudopotential file
        """
        psroot = ET.fromstringlist(iter_upf_file(filename))

        # PP_INFO
        try:
            self.pp_info = psroot.find('PP_INFO').text
        except AttributeError:
            self.pp_info = ""
        try:
            self.pp_input = psroot.find('PP_INFO/PP_INPUTFILE').text
        except AttributeError:
            self.pp_input = ""

        # PP_HEADER
        self.pp_header = dict(psroot.find('PP_HEADER').items())

        # PP_MESH
        self.pp_mesh = dict(psroot.find('PP_MESH').items())
        self.pp_r = np.array([float(x) for x in psroot.find('PP_MESH/PP_R').text.split()])
        self.pp_rab = np.array([float(x) for x in psroot.find('PP_MESH/PP_RAB').text.split()])

        # PP_LOCAL
        node = psroot.find('PP_LOCAL')
        if not node is None:
            self.pp_local = np.array([x for x in map(float, node.text.split())])
        else:
            self.pp_local = None

        # PP_RHOATOM
        node = psroot.find('PP_RHOATOM')
        if not node is None:
            self.pp_rhoatom = np.array([v for v in map(float, node.text.split())])
        else:
            self.pp_rhoatom = None

        # PP_NONLOCAL
        node = psroot.find('PP_NONLOCAL')
        if not node is None:
            betas = list()
            dij = None
            pp_aug = None
            pp_q = None
            for el in node:
                if 'PP_BETA' in el.tag:
                    beta = dict(el.items())
                    val = np.array([x for x in map(float, el.text.split())])
                    beta.update(dict(beta=val))
                    betas.append(beta)
                elif 'PP_DIJ' in el.tag:
                    text = '\n'.join(el.text.strip().split('\n')[1:])
                    dij = np.array([x for x in map(float, text.split())])
                elif 'PP_AUGMENTATION' in el.tag:
                    pp_aug = dict(el.items () )
                    pp_qijl = list()
                    pp_qij  = list()
                    for q in el:
                        if 'PP_QIJL' in q.tag:
                            qijl = dict( q.items() )
                            val = np.array( [ x for x in map(float, q.text.split())])
                            qijl.update(dict(qijl = val))
                            pp_qijl.append(qijl)
                        elif 'PP_QIJ' in q.tag:
                            qij = dict(q.items() )
                            val = np.array( [x for x in map(float,q.text.split())])
                            qij.update(dict(qij = val))
                            pp_qij.append(qij)
                        elif q.tag =='PP_Q':
                            pp_q = np.array( [x for x in map(float, q.text.split() )])
                    pp_aug.update(dict(PP_QIJL=pp_qijl, PP_QIJ = pp_qij, PP_Q = pp_q) )
            # TODO: check this is the best way to store these data
            self.pp_nonlocal = dict(PP_BETA = betas, PP_DIJ = dij, PP_AUGMENTATION = pp_aug )
        else:
            self.pp_nonlocal = None

    def get_pseudo_charge(self, G, tau, spin):
        pass

