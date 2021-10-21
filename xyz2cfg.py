# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:38:09 2021

@author: Daisy

To use it:python xyz2cfg.py xyz.filename

Atention: no viral forces written by this script
"""

import re
import sys
import os

if len(sys.argv) <=1:
    raise SystemExit("wrong need filename")
with open (sys.argv[1]) as reader:
    os.mkdir("cfgs")
    PATH=os.getcwd()
    PATH=os.path.join(PATH,"cfgs")
    string = reader.readlines()
    for index,i in enumerate(string):
        if 'config_type=amorph' in i:
            atom_list = string[index-1]
            lattice_list= string[index]
            lattice = re.findall('Lattice="(.*?)"',lattice_list)
            vector = lattice[0].split( )
            energy = re.findall('dft_energy=(.*?) ',lattice_list)
            coordinate_list = []
            tmp = os.path.join(PATH,"%d" %index)
            os.mkdir(tmp)
            tmp = os.path.join(tmp,'single.cfg')
            file = open(tmp, 'w')
            for j in range(1,int(atom_list)+1):
                coordinate_list.append(string[index+j])
                coordinate2write = []
                for m,n in enumerate(coordinate_list):
                    b = n.split( )
                    coordinate = ' '.join([b[1],b[2],b[3],b[-6],b[-5],b[-4]])
                    coordinate_2 = "             "+str(m+1)+"    1       "+coordinate
                    coordinate2write.append(coordinate_2)
            file.write("BEGIN_CFG\n Size\n")
            file.write("    %s" %(atom_list))
            file.write(" Supercell\n")
            file.write('         %s %s %s\n' %(vector[0],vector[1],vector[2]))
            file.write('         %s %s %s\n' %(vector[3],vector[4],vector[5]))
            file.write('         %s %s %s\n' %(vector[6],vector[7],vector[8]))
            file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
            for line in coordinate2write:
                file.write(line+'\n')
            file.write(" Energy\n")
            file.write("       %s\n" %(energy[0]))
            file.write(" END_CFG\n \n")
