# -*- coding: utf-8 -*-
#data2poscar
import sys

print("Function: python data2poscar.py atom_type_0 atom_type_1 ……")
with open("init.data","r") as reader,open("POSCAR","w") as writer:
    # element_0 = sys.argv[0]
    # element_1 = sys.argv[1]
    element_0 = "Li"
    element_1 = "Si"
    lines = reader.readlines()
    for i,j in enumerate(lines):
        line = j.split( )
        if "atoms" in line:
            total_num = int(line[0])
        elif "xlo" in line:
            x_coordinate = eval(line[1])
        elif "ylo" in line:
            y_coordinate = eval(line[1])
        elif "zlo" in line:
            z_coordinate = eval(line[1])
        elif "Atoms" in line:
            coordinates = []
            count_0 = 0
            count_1 = 0
            for k in range(total_num):
                coordinates.append(lines[i+k+2])
                atom_type = lines[i+k+2].split()[1]
                if int(atom_type) == 1:
                    count_0 += 1
                elif int(atom_type) == 2:
                    count_1 += 1
    writer.write("POSCAR by daisy\n")
    writer.write("1.0\n")
    writer.write('%.3f  0.0000  0.0000\n' %x_coordinate)
    writer.write("0.0000  %.3f  0.0000\n" %y_coordinate)
    writer.write("0.0000  0.0000  %.3f\n" %z_coordinate)
    writer.write("%s %s\n" %(element_0,element_1))
    writer.write("%d %d\n" %(count_0,count_1))
    writer.write("Cart\n")
    for i in range(len(coordinates)):
        x=eval(coordinates[i].split()[2])
        y=eval(coordinates[i].split()[3])
        z=eval(coordinates[i].split()[4])
        writer.write('%.3f %.3f %.3f\n' %(x,y,z))
