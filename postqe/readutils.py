#encoding: UTF-8

import time, sys
import numpy as np
from struct import *

################################################################################
# Readers of binary files.
#
# Warning: all the following functions for reading binary files are platform 
# dependent and may not work properly. Note that the long term solution is to 
# use hdf5 format.
def read_line(f):
    """
    Read a line from a binary file. End character must be b'\n'
    """
    line = b''
    byte = f.read(1)
    while (byte != b'\n'):
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
    import re
    
    pat = re.compile(r'nr\d="\d+"')    
    info_line_nrx = pat.findall(info_line)
    
    nr = np.zeros(3,dtype=int)
    for i in range(0,len(info_line_nrx)):
        nr[i] = info_line_nrx[i].split("\"")[1]
    
    return nr[0], nr[1], nr[2]



def read_charge_file_iotk(fname):
    """
    Read a binary charge file written with QE and iotk. Warning: platform dependent,
    use hdf5 format when available
    """
 
    tempcharge = []
    
    with open(fname, "rb") as f:
        line = b''
        while (line!=b'  <CHARGE-DENSITY>'):    # skip lines till the <CHARGE-DENSITY> tag is found
            line = read_line(f)
            
        discard = read_line(f)  # discard binary line here
        line = read_line(f)     # this is the line <INFO nr1=.. />
        info_line = line.decode('utf-8')
        nr1, nr2, nr3 = get_info(info_line)
        #print ("Reading from charge file: nr1= ",nr1,"nr2= ",nr2,"nr3= ",nr3)

        # data are grouped in nr1*nr2 sequential reals in the charge file, each
        # starting with <z.1, <z.2, etc. tag and ending with a corresponding
        # </z.1>, </z.2> containing extra information...
        charge = np.zeros((nr1,nr2,nr3))
        for i in range(0,nr3):   
            discard = read_line(f)
            line = read_line(f)     # this is the line <z.1 type="real" size="2025" kind="8"> or similar
        
            temp = read_n_real_numbers(f,nr1,nr2)
            charge[:,:,i] = temp
            
            discard = read_line(f)
            line = read_line(f)     # this is end tag </z.1>
        
    return charge


def read_wavefunction_file_iotk(fname):
    """
    To be implemented.
    """

    return None


def read_charge_file_hdf5(fname):
    """
    Reads an hdf5 charge file written with QE. nr1, nr2, nr3 (the dimensions of
    the charge k-points grid) are read from the charge file.
    """
 
    import h5py
    
    f = h5py.File(fname, "r")
    nr1 = f.attrs.get("nr1")
    nr2 = f.attrs.get("nr2")
    nr3 = f.attrs.get("nr3")
    charge = np.zeros((nr1,nr2,nr3))
    for i in range(0,nr3):   
        dset_label = "K"+str(i+1)   # numbered from 1 to nr3 in file hdf5
        tempcharge = np.array(f[dset_label])
        charge[:,:,i] = np.reshape(np.array(tempcharge),(nr1,nr2))
    
    return charge


def read_wavefunction_file_hdf5(fname):
    """
    Reads an hdf5 wavefunction file written with QE. Returns a dictionary with
    the data structure in the hdf5 file. 
    """
 
    import h5py
    
    f = h5py.File(fname, "r")    
    nkpoints = len(f["KPOINT1"].attrs.values())
    #print ("nkpoints = ",nkpoints)

    wavefunctions = {}

    for i in range(0,nkpoints):
        temp = {}
        kpoint_label = "KPOINT"+str(i+1)
        # read the attributes at each kpoint
        attrs_to_read = ["gamma_only", "igwx", "ik", "ispin","ngw","nk","nbnd","nspin","scale_factor"]
        for attr in attrs_to_read:
            temp[attr] = f[kpoint_label].attrs.get(attr)
        for iband in range(0,temp["nbnd"]):
            band_label = "BAND"+str(iband+1)
            temp[band_label] = np.array(f[kpoint_label][band_label])
        
        wavefunctions[kpoint_label] = temp
        
    return wavefunctions
    
    
################################################################################
# Reader of pseudopotential files and auxiliary functions.
################################################################################
# Warning: limited functionalities for now
#
def get_tag(line):
    """
    Check if the input string is an xml tag, i.e. it contains between "<" and ">" the "/" and any alphanumeric character.
    If so, returns the string between  "<" and ">". If "/" is present, it is in [0] position in the returned string.
    If not a tag, returns None
    """

    import re

    line = line.strip(" \n\t")         # first strip spaces, endline chars, tabs
    pat = re.compile(r'<\/?[\w\ =\+\-\.\"]+\/?>')    # a tag may contain between "<" and ">" the "/" and any alphanumeric character
    if pat.match(line):
        return line.strip(" <>")     # returns the tag without "<" and ">"

################################################################################
# Read tags in list_tags and returns a dictionary
#
def read_tags(lines,list_tags):
    """
    This function loops over the lines in input and returns a dictionary 
    with the content of each tag, but only if they are in list_tags
    """

    dic = {}
    content = ""
    i = 0
    add_content = False
    while (i<len(lines)):
        tag = get_tag(lines[i])
        if ( tag!=None): tag = tag.split()[0]
        if ((tag!=None) and (tag[0]!="/") and (tag  in list_tags)):
            #print ("Reading tag... ",tag)
            add_content = True
        elif ((tag!=None) and (tag[0]=="/") and (tag.strip("/") in list_tags)):
            tag = tag.strip("/")
            #print ("Closing tag... ",tag)
            dic[tag] = content
            content = ""
            add_content = False
        elif ( add_content):
            content += lines[i]
        i += 1
    return dic

def read_upf2 (doc):
    """
    extracts pseudo data from upf v.2 files
    :param doc: doc is the ElementTree object representing the pseudo file in xml format
    :return: a dictionary 
    """
    pseudo ={}
    psroot = doc.getroot()
    #PP_INFO
    pp_info = psroot.find('./PP_INFO').text
    pp_input = psroot.find('./PP_INFO/PP_INPUTFILE').text
    pseudo.update(PP_INFO=dict(INFO=pp_info,PP_INPUT = pp_input))
    #
    #PP_HEADER
    pp_header = dict(psroot.find('./PP_HEADER').items())
    pseudo.update(PP_HEADER=pp_header)
    #PP_MESH
    pp_mesh = dict(psroot.find('./PP_MESH').items() )
    pp_r   = np.array([float(x) for x in psroot.find('./PP_MESH/PP_R').text.split() ])
    pp_rab = np.array([float(x) for x in psroot.find('./PP_MESH/PP_RAB').text.split()])
    pp_mesh.update(PP_R=pp_r, PP_RAB = pp_rab)
    pseudo.update(PP_MESH = pp_mesh)
    #PP_LOCAL
    pp_local = None
    node = psroot.find('./PP_LOCAL')
    if not node is None:
        pp_local = np.array([x for x in map(float, node.text.split() )])
    pseudo.update(PP_LOCAL = pp_local )
    #
    node = psroot.find('./PP_RHOATOM')
    if not node is None :
        pp_rhoatom = np.array([v for v in map(float,node.text.split()) ])
    else:
        pp_rhoatom = None
    pseudo.update(PP_RHOATOM=pp_rhoatom)
    #
    node = psroot.find('./PP_NONLOCAL')
    if not node is None:
        betas = list()
        dij = None
        pp_aug = None
        elements = node.getchildren()
        for el in elements:
            if 'PP_BETA' in el.tag:
                beta = dict( el.items() )
                val  = np.array([x for x in map(float,el.text.split() )])
                beta.update(beta = val)
                betas.append(beta)
            elif 'PP_DIJ' in el.tag:
                dij = np.array([ x for x in map(float,el.text.split() )])
            elif 'PP_AUGMENTATION':
                pp_aug = dict(el.items () )
                aug_elements = el.getchildren()
                pp_qijl = list()
                pp_qij  = list()
                for q in aug_elements:
                    if 'PP_QIJL' in q.tag:
                        qijl = dict( q.items() )
                        val = np.array( [ x for x in map(float, q.text.split() )])
                        qijl.update(qijl = val)
                        pp_qijl.append(qijl)
                    elif 'PP_QIJ' in q.tag:
                        qij = dict(q.items() )
                        val = np.array( [x for x in map(float,q.text.split() )])
                        qij.update(qij = val)
                        pp_qij.append(qij)
                    elif q.tag =='PP_Q':
                        pp_q = np.array( [x for x in map(float, q.text.split() )])
                pp_aug.update(PP_QIJL=pp_qijl, PP_QIJ = pp_qij, PP_Q=pp_q)
        pp_nonlocal = dict(PP_BETA = betas, PP_DIJ = dij, PP_AUGMENTATION = pp_aug )
    else:
        pp_nonlocal = None
    pseudo.update(PP_NONLOCAL = pp_nonlocal)
    return pseudo

def read_pseudo_file(fname):
    from xml.etree import ElementTree as ET
    from os  import remove
    """
    This function reads a pseudopotential file in the QE UPF format (text), returning
        the content of each tag in a dictionary. 
    WARNING: does not handle multiple tags with the same name yet and has limited 
    functionalities for now. It is meant to be used only for postqe, not as a general
    reader of pseudopotentials files.
    """

    list_tags = ["PP_INFO","PP_HEADER","PP_MESH","PP_NLCC", "PP_LOCAL","PP_NONLOCAL","PP_PSWFC","PP_RHOATOM"]
    with open (fname, 'r') as psfile, open('temp.pseudo','w') as temp:
        for line in psfile:
            temp.write(line.replace('&input','&amp;input'))
    try:
        psdoc = ET.parse('temp.pseudo')
        remove('temp.pseudo')
        return read_upf2(psdoc)
    except ET.ParseError:
        remove('temp.pseudo')
    lines = []
    with open(fname, "r") as temp:
        for line in temp:
            lines.append(line)
    pseudo = read_tags(lines,list_tags)
    ### convert strings to float numpy arrays
    rloc = np.array( [ float(x) for x in pseudo["PP_LOCAL"].split()] )
    pseudo["PP_LOCAL"] = rloc
    #### Read subfields
    list_tags_PP_MESH = ["PP_R","PP_RAB"]
    pseudo["PP_MESH"] = read_tags(pseudo["PP_MESH"].splitlines(),list_tags_PP_MESH)
    rr =  np.array( [float(x) for x in pseudo["PP_MESH"]["PP_R"].split() ] )
    rrab= np.array ([ float (x) for x in pseudo["PP_MESH"]["PP_RAB"].split()  ] )
    tempdict = dict(PP_R=rr, PP_RAB=rrab )
    pseudo["PP_MESH"] = tempdict

    rhoat = np.array([ float(x) for x in pseudo["PP_RHOATOM"].split() ] )
    pseudo["PP_RHOATOM"] = rhoat
    
    list_tags_PP_NONLOCAL = ["PP_BETA","PP_DIJ"] 
    pseudo["PP_NONLOCAL"] = read_tags(pseudo["PP_NONLOCAL"].splitlines(),list_tags_PP_NONLOCAL)
        
    return pseudo


################################################################################
# Other readers, writers, auxiliary functions.
################################################################################
def write_charge(fname,charge,header):
    """
    Write the charge or another quantity calculated by postqe into a file fname.    
    """
    
    fout = open(fname, "w")
    
    # The header contains some information on the system, the grid nr, etc.
    fout.write(header)
    
    nr = charge.shape
    count = 0
    # Loop with the order as in QE files
    for z in range(0,nr[2]):
        for y in range(0,nr[1]):
            for x in range(0,nr[0]):
                fout.write("  {:.9E}".format(charge[x,y,z]))
                count += 1
                if (count%5==0):
                    fout.write("\n")
                        
    fout.close()
    
    
def create_header(prefix,nr,ibrav,celldms,nat,ntyp,atomic_species,atomic_positions):
    """
    Creates the header lines for the output file. A few fields are different from QE.
    """

    text = "# "+prefix+"\n"
    text += "# {:8d} {:8d} {:8d} {:8d} {:8d} {:8d} {:8d} {:8d}\n".format(nr[0],nr[1],nr[2],nr[0],nr[1],nr[2],nat,ntyp)
    text += "# {:6d}    {:8E}  {:8E}  {:8E}  {:8E}  {:8E}  {:8E}\n".format(ibrav,*celldms)
    text += "#      "+4*"XXXX   "+"\n"
    
    ityp = 1
    for typ in atomic_species:
        text += "# {:4d} ".format(ityp)+typ["@name"]+"\n"
        ityp += 1

    ipos = 1
    for pos in atomic_positions:
        text += "# {:4d}  ".format(ipos)
        coords = [float(x) for x in pos['#text'] ]
        text += " {:9E} {:9E} {:9E}\n".format(*coords)
        text += pos["@name"]+"\n"
        ipos += 1
    
    return text
    
 
def read_postqe_output_file(fname):
    """
    This function reads the output charge (or other quantity) as the output 
    format of postqe. 
    """
    
    tempcharge = []
    count = 0
    nr = np.zeros(3,dtype=int)
    with open(fname, "r") as lines:
        for line in lines:
            linesplit=line.split()
            if count==1:
                nr[0] = int(linesplit[1])
                nr[1] = int(linesplit[2])
                nr[2] = int(linesplit[3])
            if linesplit[0]!='#':                       # skip the first lines beginning with #
                for j in range(0,len(linesplit)):       # len(linesplit)=5 except maybe for the last line
                    tempcharge.append(float(linesplit[j]))
            count += 1
   
    
    charge = np.zeros((nr[0],nr[1],nr[2]))
    count = 0
    # Loops according to the order it is done in QE
    for z in range(0,nr[2]):
        for y in range(0,nr[1]):
            for x in range(0,nr[0]):
                charge[x,y,z] = tempcharge[count]
                count += 1

    return charge
    
    
    
###########################################
#
# This is only for testing the functions in this module
if __name__ == "__main__":
    prefix = "../tests/"
    from readutils import read_pseudo_file
    pseudo = read_pseudo_file(prefix+"Al.pz-vbc.UPF")
 
    print ("PP_INFO\n")
    print (pseudo["PP_INFO"])
    print ("PP_MESH\n")
    print (pseudo["PP_MESH"])
    
    print (pseudo["PP_MESH"]["PP_R"])
    print (pseudo["PP_MESH"]["PP_RAB"])
    