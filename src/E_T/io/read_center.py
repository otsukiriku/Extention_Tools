import numpy as np

def r_c:
    center=open("./center.txt", mode='r')
    line = center.readline()
    while True:
        line = center.readline().split()
        if not line:
            break
        elif len(line)==4:
            CB_mask.append(line.pop(0))
            CB_g.append(line)
            #print(line)
    CB_mask=list(map(int, CB_mask))
    #print(CB_mask)
    CB_G=[]
    for line in CB_g:
        CB_G.append(list(map(float, line)))
    CB_G=np.array(CB_G)
    return CB_G
    #print(CB_G)
