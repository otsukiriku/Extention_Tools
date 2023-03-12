
import numpy as np
import sys
from MolCop import mmpystream as mmps
from MolCop.analysis import topology


def create_molnum(atoms:mmps.Stream().sdat):
    connect_list = atoms.connect_list
    m_list = topology.create_molecule(atoms)
    if not 'molnum' in atoms.particles:
        atoms.add_particles_property("molnum", _dtype=int, dim=1)
    else:
        pass

    for l in m_list:
        num_mol = len(l)
        for ind in l:
            atoms.particles["molnum"][ind]=num_mol

#    return atom
