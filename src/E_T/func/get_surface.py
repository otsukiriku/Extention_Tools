
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import numpy as np
import sys
import copy
import math

def surface(atoms:mmps.Stream().sdat, coordination_cutoff, _moloption = False)
    atoms.create_connect_list(0.3)
    connect_list = atoms.connect_list

    target = copy.deepcopy(atoms)

    #配位数でtargetの表面を管理することにした．
    coordination_num=[]
    coordination_num.extend([len(i) for i in target.connect_list])
    coordination_num = np.array(coordination_num)

    flag = (coordination_num < coordination_cutoff)
    target.flagconnect = True
    target.trimming_particles(flag, reindex=True)
    
    if _moloption == True:
        m_list = topology.create_molecule(target)
    #trim only target particle
    return target
