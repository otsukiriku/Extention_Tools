from MolCop import mmpystream as mmps
from MolCop.analysis import topology
#from evalueate_structure import get_center as gc
import numpy as np
import sys
import copy
import math

def g_c_by_mask(atoms : mmps.Stream().sdat.particles, mask):
    flag = atoms["mask"] == mask
    center = atoms["pos"][flag].mean(axis=0)
    return center

#atom1 should be larger one
def get_between_dis(atom1, atom2):
#    for i in atom2
    dis = (atom1 - atom2)**2
    disarr = np.sum(dis,axis = 1)
    return disarr

def Pt_location(atoms:mmps.Stream().sdat, nonzero_option=True, prove_range=60,CB_num=4, Pt_num=56, Pt_cluster_pnum=309, Pt_mask_start=7):
    # get center
    dcell=atoms.cell[2]
    pos_shift = [0, 0, dcell]
    CB_list=[i+1 for i in range(CB_num)]
    CB_G = np.empty((0,3))
    for mask in CB_list:
        #for wrap particles
        CB = atoms.particles
        flag = ((CB['mask']==1) & (CB['pos'][:,2]>dcell*2/3))
        CB['pos'][flag]-=pos_shift
        flag = ((CB['mask']==CB_num) & (CB['pos'][:,2]<dcell*1/3))
        CB['pos'][flag]+=pos_shift

        flag = CB["mask"] == mask
        #print(CB["mask"]) 
        #print(flag) 
        CB_g = [CB["pos"][flag].mean(axis=0)]
        CB_G = np.append(CB_G, CB_g, axis=0)
    #print("CB_G")
    #print(CB_G)
    bottom = CB_G[CB_num-1,] - np.array(pos_shift)
    CB_G = np.append(CB_G, bottom)
    top = CB_G[0,] + np.array(pos_shift)
    CB_G = np.append(CB_G, top)
    CB_G = np.reshape(CB_G,(CB_num+2,3))

    #get each Ptcluster gravity center
    Pt=atoms
    Pt_list=[i for i in range(Pt_mask_start,Pt_mask_start+Pt_num)]
    Ptmask=Pt.particles['mask'][Pt.particles['type'] == 4]
    Ptmask_2d=np.reshape(Ptmask, (Pt_num, Pt_cluster_pnum))

    for idx, msk in enumerate(Ptmask_2d):
        #一列ずつPtmaskから切り取る．
        msk[:] = Pt_list[idx]
    #print(Ptmask)

    atoms.particles['mask'][Pt.particles['type'] == 4] = Ptmask

    Pt_list=np.array(Pt_list)
    Pt_G = np.array([g_c_by_mask(Pt.particles, l) for l in Pt_list])

    #get closest & 2nd closest dis from center
    # original
    #m_l1 = [i+1 for i in range(0,CB_num)]
    m_l1 = [i+1 for i in range(0,CB_num+2)]
    mask_l = [m_l1 for j in range(Pt_num)]
    Pt_CB_dis=np.array([get_between_dis(j, CB_G) for j in Pt_G])
    #arguments of mask of closest CB from Pt
    sort = np.sort(Pt_CB_dis, axis=1)
    #print(sort)
    s_mask=np.empty((0,1), dtype=int)
    s_dis=np.empty((0,1))
    l_mask=np.empty((0,1),dtype=int)
    l_dis=np.empty((0,1))
    c_dis=np.empty((0,1))
    for idx, arg in enumerate(sort):
        #各Ptから最小のCBのmaskを取り出す．
        #print( np.where(Pt_CB_dis[idx:idx+1][0]==arg[0])[0][0])
        s_mask_idx = np.where(Pt_CB_dis[idx,:]==arg[0])[0][0]
        s_dis = np.append(s_dis,Pt_CB_dis[idx,s_mask_idx])
        s_mask = np.append(s_mask, s_mask_idx)
        #二番目に距離が短いmaskとかを取り出す．
        l_mask_idx = np.where(Pt_CB_dis[idx,:]==arg[1])[0][0]
        l_dis = np.append(l_dis,Pt_CB_dis[idx,l_mask_idx])
        l_mask = np.append(l_mask, l_mask_idx)
        #近い二つCBの中心間距離c_dis
        #Ptの重心から最近接，第二近接CBの中間の座標までの距離を出す
        CB_between_pos = ( CB_G[s_mask_idx] + CB_G[l_mask_idx] )/2
        Pt_to_CBbet = CB_between_pos - Pt_G
        
    dis = ( CB_G[s_mask,:] - CB_G[l_mask,:])**2
    c_dis = np.sum(dis,axis = 1)
    #c_dis = np.append(c_dis, disarr)

    s_dis=np.sqrt(s_dis)
    l_dis=np.sqrt(l_dis)
    c_dis=np.sqrt(c_dis)

    #prove radius
    CB_to_Pt_vec=((Pt_G-CB_G[s_mask,:]))
    #print(CB_to_Pt_vec)
    pr = s_dis*c_dis**2/(-l_dis**2+s_dis**2+c_dis**2) - s_dis

    #get prove center(for test)
    scalor=(s_dis+abs(pr))/s_dis
    scalor=np.reshape(scalor, (Pt_num,1))
    correct = np.sum(s_dis)/Pt_num-100
    #print(correct)
    pr = pr + correct
    pr_pos = CB_G[s_mask,:] + CB_to_Pt_vec*(scalor)
    count = []
    [count.append(i) for i in range(0,100,1)]
    freq = []
    mask = [i for i in range(Pt_mask_start, Pt_mask_start+Pt_num)]
    pr_list = []
    for i in pr:
        if (0 < i) & ( i < prove_range ):
            pr_list.append(i)
        else:
            if nonzero_option == False:
                pr_list.append(0.0)
            elif nonzero_option == True:
                pass

    return mask, pr_list
