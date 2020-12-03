# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 11:46:22 2020

@author:Daisy
"""


import re
import os
import numpy as np 

#导出 gamma > 2的帧数
s = np.loadtxt("a.txt")
m = s.tolist()
n = ""
for i in range(len(m)):
    if m[i] > 2:
        n += str(i)+" "
    else:
        pass  
o = n.split(" ")

#根据索引在out.cfg中导出对应的结构并保存
with open ("out.cfg") as reader:
    string = reader.read()
    pattern = r'BEGIN_CFG.*?END_CFG'
    tmp = re.findall(pattern,string,re.S)
    for i in range(len(o)-1):
        string2write = tmp[int(o[i])]
        dirname = "%d" %i
        os.mkdir(dirname)
        path = os.path.join(dirname,"single.cfg")
        with open(path,mode = "w") as writer:
            writer.write(string2write +"\n")
