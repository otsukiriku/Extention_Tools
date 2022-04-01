
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
#from evalueate_structure import get_center as gc
from E_T.io import read_center as rc
import numpy as np
import sys
import copy
import math



def add_probe(atoms:mmps.Stream().sdat, tar_pos:np.array, tar_velo = False):
    if atoms.particles['pos'].ndim != tar_pos.ndim:
        print("dimension should be same")
    atom_num=atoms.total_particle
    add_idx = [i+atom_num+1 for i in range(len(tar_pos))]
    atoms.particles['id']=np.append(atoms.particles['id'],np.array(add_idx))
    atoms.particles['pos']=np.append(atoms.particles['pos'],tar_pos,axis=0)
    if not tar_velo:
        tar_velo = np.zeros_like(tar_pos)
    atoms.particles['pos']=np.append(atoms.particles['pos'],tar_velo,axis=0)
    atoms.particles['type']=np.append(atoms.particles['type'],10)
    if 'mask' in atoms.particles:
        atoms.particles['mask']=np.append(atoms.particles['mask'],0)
    atoms.total_particle += len(tar_pos)
