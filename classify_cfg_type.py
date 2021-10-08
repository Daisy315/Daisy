# -*- coding: utf-8 -*-
import re
import os 
import sys
from collections import defaultdict

print("To use it : python classify_cfg_type.py cfgname")
if (len(sys.argv) <2):
    raise IOError("cfg filename needed")
with open(sys.argv[1],'r') as reader:
# with open("shuffled.cfg",'r') as reader:
    counts = []
    all = reader.read()
    separtor = '\n' 
    if '\r\n' in all:
        separtor = '\r\n'
    pattern = re.compile("BEGIN_CFG.*?END_CFG",re.S)
    ret = re.findall(pattern, all)
    for i in ret:
        alllines = i.split(separtor)
        for index, j in enumerate(alllines):
            if "Size" in j:
                atomnums = int(alllines[index+1])
            elif "AtomData" in j:
                atom_start = index+1
                break
            
        indices = []
        for j in range(atom_start,atom_start+atomnums):
            atomj = alllines[j]
            indices.append(int(atomj.split()[1]))#获取单帧的所有原子序号
        atomtypes = sorted(list(set(indices)))#获取单帧的最大原子类型
        tmp = []
        for i in atomtypes:
            tmp.extend([i,indices.count(i)])#count单帧的每类原子有多少个
        counts.append(tuple(tmp))#将单帧的结果存到counts
 
index_4_same_config = defaultdict(list)
for k, va in [(v,i) for i, v in enumerate(counts)]:
    index_4_same_config[k].append(va)
for key, value in index_4_same_config.items():
    a = "".join(str(key).strip())
    b = a.split() 
    dirname = "".join(b)
    dirname = dirname.replace("(","")
    dirname = dirname.replace(",","_")
    dirname = dirname.replace(")","")
    dir_name = "data%s" %dirname
    # if os.path.exists(path):
    #     for root, dirs, files in os.walk(path, topdown=False):
    #         for name in files:
    #             os.remove(os.path.join(root, name)) #删除文件
    #         for name in dirs:
    #             os.rmdir(os.path.join(root, name)) #删除文件夹
    #     os.rmdir(path) #删除mydata文件夹
    path = dir_name
    os.mkdir(path,mode=0o755)
    content = ""
    for m in value:
        to_written = ret[m] 
        content += "%s  \r\n\n" %(to_written)
    with open(os.path.join(path,"shuffled.cfg"),'w') as writer:
        writer.write(content)
    writer.close()
