# -*- coding: utf-8 -*-
import re
import os 
import sys
import numpy as np

if (len(sys.argv) <2):
    raise IOError("cfg filename needed")
with open(sys.argv[1],'r') as reader:
    everything_one_cfg = []
    all = reader.read()
    separtor = '\n' 
    if '\r\n' in all:
        separtor = '\r\n'
    pattern = re.compile("BEGIN_CFG.*?END_CFG",re.S)
    ret = re.findall(pattern, all)
    for index_cfg, i in enumerate(ret):
        stress = None
        alllines = i.split(separtor)
        for index, j in enumerate(alllines):
            if "Size" in j:
                atomnums = int(alllines[index+1])
            elif "Supercell" in j:
                cell_start = index+1                
            elif "AtomData" in j:
                atom_start = index+1
            elif "Energy" in j:
                Energy_start = index+1
            elif "Stress" in j:
                stress = index+1                
        atoms = []
        atoms_type = []
        cells = []
        for j in range(atom_start,atom_start+atomnums):
            atom_j = alllines[j].split()
            atoms_type.append(int(atom_j[1]))
            atoms.append(list(map(float,atom_j[2:])))
        #atomtypes_nums = sorted(list(set(atoms_type)))

        for j in range(cell_start,cell_start+3):
            cell_j = alllines[j].split()
            cells.extend(list(map(float,cell_j)))
        energy = float(alllines[Energy_start])

        if stress:
            stress = list(map(float,alllines[stress].split()))
        
        everything_one_cfg.append({"atom_types":atoms_type,"atom_xyz_f":atoms,"cell":cells,"energy":energy,"stress":stress})
del alllines  
all_atomtypes = [i["atom_types"] for i in everything_one_cfg]

max_all_atomtypes = max([max(i) for i in all_atomtypes])   

all_component = ['_'.join(["%d" %(i.count(j)) for j in range(max_all_atomtypes+1)]) for i in all_atomtypes]     
component_index = {i:[] for i in list(set(all_component))}
for index,i in enumerate(all_component):
    component_index[i].append(index)

all_atoms_xyz_f = [i["atom_xyz_f"] for i in everything_one_cfg]
all_cell = np.array([i["cell"] for i in everything_one_cfg])
all_energy = np.array([i["energy"] for i in everything_one_cfg])
if stress:
    all_stress = np.array([i["stress"] for i in everything_one_cfg])

if sys.version[0] == "3":
    input = input
else:
    input = raw_input

element_type = []

for i in range(max_all_atomtypes+1):
    element_type.append(input("input element name for type %d-->" %i))  

for key,value in component_index.items():
    current_atomtypes = np.array([all_atomtypes[i] for i in value]) 
    current_energy = all_energy[np.ix_(value)]
    current_cell = all_cell[np.ix_(value)]
    tmp = np.array([all_atoms_xyz_f[i] for i in value])

    current_xyz = tmp[:,:,:3]
    current_f =  tmp[:,:,3:]

    #check element sequence
    nframes = current_atomtypes.shape[0]
    for i in range(nframes-1):
        for j in range(i+1,nframes):
            if not (current_atomtypes[i]==current_atomtypes[j]).all():
                raise SystemExit("Element sequence of component %s %d and %d is not the same" %(key,i,j))
         
    # folder = key
    # os.makedirs(folder, exist_ok = True)
    folder = "set.000"
    os.makedirs(folder, exist_ok = True)
    
    np.savetxt(os.path.join(folder, 'type.raw'),    current_atomtypes[0], fmt = '%d')
    np.savetxt(os.path.join(folder, 'type_map.raw'),    np.array(element_type), fmt = '%s')
    np.savetxt(os.path.join(folder, 'box.raw'),  current_cell)
    np.savetxt(os.path.join(folder, 'coord.raw'),   np.reshape(current_xyz,   [nframes, -1]))
    np.savetxt(os.path.join(folder, 'energy.raw'),  current_energy )
    np.savetxt(os.path.join(folder, 'force.raw'),   np.reshape(current_f,   [nframes, -1]))
    if stress:
        all_stress = all_stress[np.ix_(value)]
        stress_writer = np.array([[i[0],i[-1],i[-2],i[-1],i[1],i[-3],i[-2],i[-3],i[2]] for i in all_stress])
        np.savetxt(os.path.join(folder, 'virial.raw'), stress_writer)
