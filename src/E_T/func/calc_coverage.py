from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import read_center as rc
from E_T.func import get_center as gc
import numpy as np
import sys
import copy
import math

def distance_flag(atoms: mmps.Stream().sdat, center=[],radius=100):
    distflag = np.zeros_like(atoms.particles['id'], dtype=np.bool)
    #print(distflag)
    for cent in center:
        dist = (atoms.particles['pos'] - cent)**2
        dist = np.sum(dist, axis=1)
        #print(*dist)
        distarr = dist < radius**2
        distflag += distarr
    #flag反転し，中心から一定以上離れたものを残す．
    distflag = np.logical_not(distflag)
    return distflag

def get_coordination(atoms: mmps.Stream().sdat):
    coordination_num=[]
    coordination_num.extend([len(i) for i in atoms.connect_list])
    coordination_num = np.array(coordination_num)
    return coordination_num

def calc_coverage(tar:mmps.Stream(), atoms:mmps.Stream().sdat, mask_conditions=[],type_condition=None,cutoff=9):
    count=[]
    cover=[0 for _ in mask_conditions]
    for idx,condtion in enumerate(mask_conditions):
        if type_condition is None:
            pt_pos = atoms.particles['pos'][atoms.particles['mask'] == condtion]
        else:
            flag = (atoms.particles['mask'] == condtion) & (atoms.particles["type"] == type_condition)
            pt_pos = atoms.particles['pos'][flag]
        distflag = np.zeros_like(tar.particles['id'], dtype=np.bool)
        for p in pt_pos:
            dist = (tar.particles['pos'] - p)**2
            dist = np.sum(dist, axis=1)
            #[print(math.sqrt(d)) for d in dist if d < cutoff**2]
            flag = dist < (cutoff**2)
            sumflag = np.sum(flag)
            coverflag=0
            if sumflag != 0:
                coverflag = 1
            cover[idx]+=coverflag
            distflag += flag
        temp=np.sum(distflag)
        temp.tolist()
        count.append(temp)
    return count,cover

