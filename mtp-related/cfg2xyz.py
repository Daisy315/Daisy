# -*- coding: utf-8 -*-
"""
CONVERT cfg flie of MTP input format to xyz file of sGDml input format
@author: Daisy
"""

import sys
import os

if len(sys.argv) <=1:
     raise SystemExit("wrong need filename")
print("usage:python cfg2xyz.py cfg_filemane atom_type,for single element only")
atom_type=str(sys.argv[-1])
#atom_type="Si"
with open (sys.argv[1]) as reader:
#with open ("new.cfg") as reader:
    PATH=os.getcwd()
    PATH=os.path.join(PATH,"new.xyz")
    string = reader.readlines()
    new_to_xyz=[]
    for index,i in enumerate(string):
        if 'BEGIN_CFG' in i:
            num = int(string[index+2])
            to_xyz=[]
            for j in range(1,num+1):
                to_xyz.append(string[index+7+j])
            energy=float(string[index+num+9])
            new_to_xyz.append(num)
            new_to_xyz.append(energy)
            new_to_xyz.append(to_xyz)
    file = open(PATH, 'w')
    for index,line in enumerate(new_to_xyz):
        if (index+1) % 3 == 0:
            for i in range(63):#change here!!!value stand for atom_nums in each conf.
                file.write("%s  %s  %s  %s  %s  %s  %s\n" %(atom_type,line[i].split()[2],line[i].split()[3],line[i].split()[4],line[i].split()[5],line[i].split()[6],line[i].split()[7]))
        if (index+2) % 3 == 0:
            file.write("%6.f\n" %(line))
        if (index+3) % 3 == 0:
            file.write("%d\n" %(line))
        else:
            pass
