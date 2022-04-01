import numpy as np

def import_ryudo(filename = None):
    ifn = "./"+filename
    file=open(file=ifn ,mode='r')

    pos = []
    flag=0
    while True:
        line = file.readline().split()
        if not line:
            break
        elif line[0] == "atoms":
            total = line[1]
        elif line[0] == "axis":
            laxis=[line[1],line[2],line[3]]
            laxis = [float(i) for i in laxis]
        elif line[0] == "part":
            flag = 1
        elif flag == 1:
            pos.append([float(line[5]),float(line[6]),float(line[7])])

    return pos, laxis
    #print(CB_G)
