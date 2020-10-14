#By nxu
from pymatgen import Structure
from monty.serialization import loadfn
import os
data = loadfn('test.json')
train_structures = [d['structure'] for d in data]
train_energies = [d['outputs']['energy'] for d in data]
train_forces = [d['outputs']['forces'] for d in data]
os.mkdir("POSCARs")
PATH=os.getcwd()
PATH=os.path.join(PATH,"POSCARs")
index=1
for i in train_structures:
    tmp=os.path.join(PATH,"%d" %index)
    os.mkdir(tmp)
    tmp=os.path.join(tmp,'POSCAR')
    i.to("poscar",tmp)
    index+=1
