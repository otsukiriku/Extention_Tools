import copy
import sys
import os
import numpy as np
from MolCop import mmpystream as mmps
sys.path.append("/nfshome12/rotsuki/molcop/src")
import evaluate_structure.get_center as g_c
import modeling.unwrap_particles as u_p

def get_Pt_G(atoms:mmps.Stream().sdat, Pt_num:int, each_cluster_atoms:int):
    Pt_flag = atoms.particles['type'] == 4
    Pt = copy.deepcopy(atoms)
    Pt.trimming_particles(Pt_flag)
    #print(Pt.particles['id'])

    Pt.add_particles_property('flag', _dtype=int, dim=1)

    Pt_m_list = []
    flag = []
    for j in range(0,Pt_num):
        Pt_m_list.append(j)
        for i in range(0,each_cluster_atoms):
            flag.append(j)
    Pt.particles['flag']=np.array(flag)

    Pt_G=[]
    for i in Pt_m_list:
        flag = Pt.particles['flag'] == i
        Pt_g = g_c.g_c(Pt.particles, flag)    
        Pt_g.tolist()
        Pt_G.append(Pt_g[0])
    Pt_G = np.array(Pt_G)
    return Pt_G 

def get_CB_G(atoms:mmps.Stream().sdat, CB_num:int, ):
    CB_list = [i for i in range(1,CB_num+1)]
    CB = copy.deepcopy(atoms)
    CB_flag = CB.particles['type']==1
    CB.flagconnect = True
    CB.trimming_particles(CB_flag) 
    u_p.unwrap_p(CB)
    CB_G = np.empty((0,3))
    for mask in CB_list:
        #for wrap particles
        flag = CB["mask"] == mask
        CB_g = [CB["pos"][flag].mean(axis=0)]
        CB_G = np.append(CB_G, CB_g, axis=0)
    return CB_G

def correct_center_bymask(atoms : mmps.Stream().sdat, center:np.array, mask):
    half_cell=np.array(atoms.cell)/2
    for idx,cent in enumerate(center):
        flag=atoms.particles['mask']==mask
        over_flag = (atoms.particles['pos'][flag]-cent) > half_cell
        under_flag = (atoms.particles['pos'][flag]-cent) < -half_cell
        atoms.particles['pos'][flag] -= over_flag * atoms.newcell
        atoms.particles['pos'][flag] += under_flag * atoms.newcell
        center = g_c.g_c(atoms.particles)
        return center

