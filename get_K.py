#根据POSCAR的基矢得到K点密度并生成KPOINTS

#By nxu + Daisy

import numpy as np
import sys

tmp = sys.argv
if len(tmp)==1:
    _density = 0.03
else:
    _density = float(tmp[1])

lattice = np.zeros((3,3))
with open("POSCAR") as reader:
    reader.readline()
    scale = float(reader.readline())
    for i in range(3):
        lattice[i] = list(map(float,reader.readline().split()))
    #suppose orthogonal
    lattice *= scale
    a1 = lattice[0][0]
    a2 = lattice[1][1]
    a3 = lattice[2][2]
    volume = a1 * a2 * a3 
    # reciprocal lattice
    b1 = 2.0 * np.math.pi * (a2 * a3) / volume
    b2 = 2.0 * np.math.pi * (a1 * a3) / volume
    b3 = 2.0 * np.math.pi * (a2 * a1) / volume
    density = _density #units 2pi/A
    density *= 2.0 * np.math.pi
    K_1, K_2, K_3 = (int(round(b1/density)),int(round(b2/density)),int(round(b3/density)))
    K_1, K_2, K_3 = (max(1,K_1),max(1,K_2),max(1,K_3))
    file = open("KPOINTS", "w")
    content = "KPT-Resolved Value to Generate K-Mesh:%.2f\n" %(_density)
    content += "%s\n" %("0")
    content += "%s\n" %("Gamma")
    content += "%d %d %d\n" %(K_1,K_2,K_3)
    content += "%s\n " %("0 0 0")
    file.write(content)
    file.close()
