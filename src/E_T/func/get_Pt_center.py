import copy
import sys
import os
import numpy as np
from MolCop import mmpystream as mmps
sys.path.append("/nfshome12/rotsuki/molcop/src")
import evaluate_structure.get_center as g_c

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
