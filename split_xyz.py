#!/usr/bin/env python
#to split total qm7b xyz into single xyz file

import os

PATH=os.getcwd()
with open("QM7B_polar_dipole_CCSD.xyz","r") as reader:
    string = reader.readlines()
    count=0
    for index,i in enumerate(string):
        if 'Lattice=' in i:
            count+=1
            coordinate_list=""
            atom_num=int(string[index-1])
            for j in range(1,atom_num+1):
                coordinate_list+=string[index+j]
            tmp = os.path.join(PATH,"%04d" %count)
            os.mkdir(tmp)
            tmp = os.path.join(tmp,'frame.xyz')
            with open(tmp,"w") as writer:
                writer.write("%d\n\n" %atom_num)
                writer.write("%s" %coordinate_list)
            
# from ase.io import read
# import numpy as np

# frames = read("bak.xyz",index=":")
# for frame in frames:
    # print(frame.info)
