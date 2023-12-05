#读取XDATCAR文件并每隔目标帧数存一次POSCAR

# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 18:04:11 2020

@author: nxu Daisy

Plan A:以"Direct configuration"字样或其他明显的重复单元标志区别帧数
读取XDATCAR输入文件
读取XDATCAR前7行获取原子类型及个数，基矢等信息，并存为字符串
s = 读取到的帧数
if  所在行有"Direct configuration":
        读取到的帧数+1
        if  读取到的帧数是想要的帧数：
            往下读取“原子总数”行并将坐标写入创建的POSCAR文件
        else：
            读取“原子总数”行的内容
            
Plan B:以所在行数区分帧数,适用于每帧都只记录原子坐标的输入文件，无明显的重复单元标志
读取XDATCAR输入文件
读取XDATCAR前7行获取原子类型及个数，基矢等信息，并存为字符串
将剩下文件全部读入为列表
遍历列表，if 索引=关注帧数的起始行数（需要计算一下表达式）：
    将index in range(index+1,index+atom_nums+1)对应的Line存到一个新的字符串中
    将XDATCAR前7行和新的字符串写入POSCAR中



"""
import numpy as np
import os
import sys
import re

if len(sys.argv) == 1:
    Interval = 15
else:
    try:
        Interval = int(sys.argv[1])
    except:
        raise SyntaxError("Interval is needed.")  #获取帧数间隔
        
file = "XDATCAR"
filesize = os.path.getsize(file)
with open(file,'r') as reader:
    title = reader.readline()
    scale = float(reader.readline())
    base_vector = np.zeros([3,3],dtype = np.float64)
    for i in range(3):
        base_vector[i] = tuple(map(float,reader.readline().split())) 
    element = reader.readline().split() 
    num_str = reader.readline().split()
    num_str = [i.strip() for i in num_str] 
    num=list(map((int),num_str)) 
    atom_nums = sum(num)
    tmp = '  '.join(element) 
    tmp_num = '  '.join(num_str) #读取XDATCAR前7行
    
    content = "%s\n" %("POSCAR by daisy")
    content += "%d\n " %(scale)
    for i in range(3):
        content += "{0:.8f} {1:.8f} {2:.8f} \n".format(*(base_vector[i]))
    content += "%s\n" %(tmp)
    content += "%s\n" %(tmp_num)
    content += "%s\n" %("Direct")  #将XDATCAR前7行存入字符串content
    #PLAN A
    index = 0
    while (reader.tell() < filesize):  #当未读取到文件最后一行时一直往下执行if判断语句
        tmp = reader.readline()
        if "Direct configuration" in tmp:
            index += 1  #记录所读取的帧数
            if index % Interval == 0:  #取出需要保存为POSCAR的帧数
                file_name = re.sub(r"[^A-Za-z0-9]+","",tmp.split("=")[-1])  #将读取的帧数保存为待创建的文件名
                coordinates = np.zeros([atom_nums,3],dtype = np.float)
                path = os.path.join("POSCARs","%s") %(file_name)
                os.makedirs(path,mode=0o777,exist_ok=1)
                with open(os.path.join(path,"POSCAR"),'w') as writer:  #创建POSCARs/所提取帧数/POSCAR文件
                    writer.write(content) #首先将content（POSCAR的前七行固定格式）写入POSCAR中
                    for i in range(atom_nums):
                        coordinates[i] = tuple(map(float,reader.readline().split()))
                        writer.write("{0:.8f} {1:.8f} {2:.8f} \n".format(*(coordinates[i]))) #读取下面atom_nums行的内容并写入POSCAR中
                writer.close()
            else:
                for i in range(atom_nums): #若读取到不需要的帧数，读取原子坐标行数即可（Direct configuration一行在if外读取）
                    reader.readline()
    
    #PLAN B
    rest_content = reader.readlines()
    file_name = 0
    for index,line in enumerate(rest_content):
        coordinates_content = ""
        if (index-((Interval-1)*(atom_nums+1))) % ((atom_nums+1)*Interval) == 0: #索引是关注帧数的起始行数
            file_name += Interval  #保存关注的帧数用作文件名
            path = os.path.join("POSCARs","%d") %(file_name)
            os.makedirs(path,mode=0o777,exist_ok=1)
            with open(os.path.join(path,"POSCAR"),'w') as writer:
                writer.write(content)
                for i in range(index+1,index+atom_nums+1):  #索引在关注的重复单元内循环
                    tmp = rest_content[i] #将关注帧数的原子坐标存在tmp中
                    coordinates_content += "%s" %(tmp)
                writer.write(coordinates_content)
            writer.close()
