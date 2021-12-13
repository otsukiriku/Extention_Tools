from MolCop import mmpystream as mmps
import numpy as np
import sys
import copy

def r_c(read_file:str):
    ifn = "./"+ read_file 
    center=open(ifn, mode='r')
    gravity = []
    line = center.readline()
    while True:
        line = center.readline().split()
        if not line:
            break
        if len(line) == 4:
            mask=line.pop(0)
            gravity.append(line)
            #print(line)
    Grav=[]
    for line in gravity:
        Grav.append(list(map(float, line)))
    Grav=np.array(Grav)
    #print(CB_G)
    center.close()
    return Grav
