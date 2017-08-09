#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################

from ase.atoms import Atoms, Atom

import numpy as np

def split_atomic_symbol(x):
    import re

    regex = r"([a-zA-Z]{1,2})(\d?)"

    match = re.match(regex, x)
    if match:
        return match.groups()
    else:
        return False


def read_espresso_xml(filename):

    from ase import units
    import xmlschema

    ##########################################################
    # TODO for whatever reason this is not working now
    # schemaLoc = xmlschema.fetch_schema(filename)
    # xs = xmlschema.XMLSchema(schemaLoc)
    #
    # temporary local solution
    xs = xmlschema.XMLSchema('schemas/qes.xsd')
    ##########################################################

    print ("Reading xml file: ", filename)
    d = xs.to_dict(filename)
    dout = d["output"]
    a1 = np.array(dout["atomic_structure"]["cell"]["a1"])
    a2 = np.array(dout["atomic_structure"]["cell"]["a2"])
    a3 = np.array(dout["atomic_structure"]["cell"]["a3"])
    a_p = (dout["atomic_structure"]["atomic_positions"]["atom"])

    atoms = Atoms()

    # First define the unit cell from a1, a2, a3 and alat
    cell = np.zeros((3, 3))
    cell[0] = a1
    cell[1] = a2
    cell[2] = a3
    atoms.set_cell(cell)

    # Now the atoms in the unit cell
    for atomx in a_p:
        # TODO: extent to all possible cases the symbol splitting (for now, only numbering up to 9 work). Not a very common case...
        symbol = split_atomic_symbol(atomx['@name'])[0]
        x = float(atomx['$'][0])
        y = float(atomx['$'][1])
        z = float(atomx['$'][2])
        atoms.append(Atom(symbol, (x, y, z)))

    return atoms


if __name__ == "__main__":
    from ase.build import bulk
    from ase.visualize import view


    FeO = read_espresso_xml('feo_af.xml')
    print (FeO.get_atomic_numbers())
    print (FeO.get_cell(True))
    print (FeO.get_positions())
    view(FeO)


    print (100*'#')

    Ni = bulk('Ni', 'fcc', a=6.65)
    print (Ni.get_atomic_numbers())
    print (Ni.get_cell(True))
    print (Ni.get_positions())
    view(Ni)

    Ni2 = read_espresso_xml('Ni.xml')
    print (Ni2.get_atomic_numbers())
    print (Ni2.get_cell(True))
    print (Ni2.get_positions())
    view(Ni2)