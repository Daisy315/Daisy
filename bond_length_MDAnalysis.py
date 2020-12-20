#!/usr/bin/env python

import numpy as np
import sys
import math
import re
import MDAnalysis as mda 
from MDAnalysis.analysis.distances import dist,distance_array
from subprocess import Popen,PIPE

#https://userguide.mdanalysis.org/stable/formats/reference/pdb.html#
#cp2k pdb is not a good pdb 
_end = [-1] #for skipping end
time = []
with open("cp2k.pdb",'r') as reader, open("cp2k_format.pdb",'w') as writer:
    for index, line in enumerate(reader):
        if line.startswith("END")or line.startswith("ENDMDL"):
            _end.append(index)
        elif "time =" in line:
            time.append(line)
    reader.seek(0)  #将END和ENDMDL所在行存为结束符列表
    for i in range(1,len(_end)):
        writer.write("MODEL         %d\n" %(i))  #label the start of a frame
        for j in range(_end[i-1]+1,_end[i]):
            writer.write(reader.readline())
        reader.readline() #end
        writer.write("TER\nENDMDL\n")            #label the end of a frame
dt = 1.0
if time != []:  
    now = lambda x:re.findall(r'time\s\=(.*?)\,',x)
    dt = float(now(time[-1])[0])-float(now(time[-2])[0])
u = mda.Universe('cp2k_format.pdb',dt=dt)
Fe_group = u.select_atoms('name Fe') # 确保找到所有帧的Fe
S_group = u.select_atoms('name S') # 确保找到所有帧的S
cutoff_1 = 2.96  #CSD共价半径的1.15倍，具体见http://sobereva.com/414和Dalton Trans., 2008, 2832-2838
cutoff_2 = 3.496
cutoff_3 = 2.415

count = 0
total_bond_num_1 = 0
total_length_1 = 0
total_bond_num_2 = 0
total_length_2 = 0
total_bond_num_3 = 0
total_length_3 = 0
for ts in u.trajectory[:]:
    s1=distance_array(Fe_group.positions, S_group.positions, box=ts.dimensions) #计算原子A在每帧中与其他B原子的距离，返回10271个*32个ndarray，每个ndarray有64个元素。len(s1)=32，len(s1[1])=64，循环有10271次
    s2=distance_array(Fe_group.positions, Fe_group.positions, box=ts.dimensions)
    s3=distance_array(S_group.positions, S_group.positions, box=ts.dimensions)
    for i in range(len(s1)):
        acceptable_donors_1 = np.where((s1[i] < cutoff_1) & (s1[i] > 0))
        for j in acceptable_donors_1:
            acc_distance_to_1 = s1[i][j]
            length_1 = np.sum(acc_distance_to_1)
        bond_num_1 = len(acc_distance_to_1)
        total_length_1 += length_1 
        total_bond_num_1 += bond_num_1
    for m in range(len(s2)):
        acceptable_donors_2 = np.where((s2[m] < cutoff_2) & (s2[m] > 0))
        for n in acceptable_donors_2:
            acc_distance_to_2 = s2[m][n]
            length_2 = np.sum(acc_distance_to_2)
        bond_num_2 = len(acc_distance_to_2)
        total_length_2 += length_2 
        total_bond_num_2 += bond_num_2
    for j in range(len(s3)):
        acceptable_donors_3 = np.where((s3[j] < cutoff_3) & (s3[j] > 0))
        for k in acceptable_donors_3:
            acc_distance_to_3 = s3[j][k]
            length_3 = np.sum(acc_distance_to_3)
        bond_num_3 = len(acc_distance_to_3)
        total_length_3 += length_3 
        total_bond_num_3 += bond_num_3
    count += 1
    print(count)

bond_length_1 = total_length_1/total_bond_num_1  #Fe-Fe
bond_length_2 = total_length_2/total_bond_num_2  #Fe-S
bond_length_3 = total_length_3/total_bond_num_3  #S-S
def print_function(group,bond_num,length):
    if bond_num == 0:
        print("there is no bond in %s"  %(group))
    else:
        print("%s内键的平均键长是: %f  A"  %(group,length))
print_function("Fe_to_S_group",total_bond_num_1,bond_length_1)
print_function("Fe_to_Fe_group",total_bond_num_2,bond_length_2)
print_function("S_to_S_group",total_bond_num_3,bond_length_3)
