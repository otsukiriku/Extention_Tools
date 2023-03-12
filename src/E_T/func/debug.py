
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
#from evalueate_structure import get_center as gc
from E_T.io import read_center as rc
import numpy as np
import sys
import copy
import math



def visualize_by_add_atom(atoms:mmps.Stream().sdat, tar_pos:np.array, tar_velo = False, tar_type=False, tar_mask = False):
    if atoms.particles['pos'].ndim != tar_pos.ndim:
        print("dimension should be same")

    atom_num=atoms.total_particle

    add_idx = [i+atom_num+1 for i in range(len(tar_pos))]
    atoms.particles['id']=np.append(atoms.particles['id'],np.array(add_idx))

    atoms.particles['pos']=np.append(atoms.particles['pos'],tar_pos,axis=0)

    if 'velo' in atoms.particles:
        if not tar_velo:
            tar_velo = np.zeros_like(tar_pos)
        atoms.particles['velo']=np.append(atoms.particles['velo'],tar_velo,axis=0)
    if 'type' in atoms.particles:
        if type(tar_type) is not np.ndarray:
            type_list = [10 for _ in range(len(tar_pos))]
        else:
            type_list = tar_type
        atoms.particles['type']=np.append(atoms.particles['type'],type_list)

    if 'mask' in atoms.particles:
        if not tar_mask:
            mask_list = [0 for _ in range(len(tar_pos))]
        else:
            mask_list = tar_mask
        atoms.particles['mask']=np.append(atoms.particles['mask'],np.array(mask_list))

    atoms.total_particle += len(tar_pos)


def create_cellmargin(atoms:mmps.Stream().sdat, margin):
    #topology.unwrap_molecule(atoms)
    pos_max = np.max(atoms.particles["pos"],axis=0)
    pos_min = np.min(atoms.particles["pos"],axis=0)
    pos_max += margin
    pos_min -= margin
    atoms.cell = [pos_max[i] for i in range(0,3)]
    atoms.dcell = [pos_min[i] for i in range(0,3)]
    for i in range(0,3):
        atoms.newcell[i] = atoms.cell[i] - atoms.dcell[i]
    atoms.shift_particles()
    atoms.wrap_particles()
    return atoms

def create_xycellmargin(atoms:mmps.Stream().sdat, margin):
    #topology.unwrap_molecule(atoms)
    pos_max = np.max(atoms.particles["pos"],axis=0)
    pos_min = np.min(atoms.particles["pos"],axis=0)
    pos_max += margin
    pos_min -= margin
    for i in range(0,2):
        atoms.cell[i] = pos_max[i]
        atoms.dcell[i] = pos_min[i]
        atoms.newcell[i] = atoms.cell[i] - atoms.dcell[i]
    atoms.shift_particles()
    atoms.wrap_particles()
    return atoms
