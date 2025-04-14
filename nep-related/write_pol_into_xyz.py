import ase.io
from ase.io import read
import shutil

with open("single-polar-raman-1_0.data") as reader:
    all_lines = reader.readlines()
    tensor = []
    for index, line in enumerate(all_lines):
        if "POLARIZABILITY TENSOR (atomic units):" in line:
            break
    xx, yy, zz = list(map(float, all_lines[index + 1].split()[1:]))
    xy, xz, yz = list(map(float, all_lines[index + 2].split()[1:]))
    yx, zx, zy = list(map(float, all_lines[index + 3].split()[1:]))

with open("cp2k.out") as reader:
    content = reader.readlines()
    for line in content:
        if "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:" in line:
            energy = float(line.split()[-1]) * 27.2114

# 复制原始文件到新文件
shutil.copyfile("mol.xyz", "new.xyz")

# 读取新文件
atoms = ase.io.read('new.xyz')

# 替换 "energy" 标签后的值
atoms.info["energy"] = f"{energy}"

# 在注释行处增加 "pol" 标签
atoms.info["pol"] = f"{xx} {xy} {xz} {yx} {yy} {yz} {zx} {zy} {zz}"

# 重新写入文件
ase.io.write("new.xyz", atoms)


# Read polarization data from "out" file
def read_alpha_from_out(num):
    with open("out") as reader:
        content = reader.readlines()
        start_index = None
        for i, line in enumerate(content):
            if "covCN         q      C6AA      α(0)" in line:
                start_index = i + 1
                break
        if start_index is not None:
            alpha_data = []
            for j in range(start_index, start_index + num):
                alpha_values = content[j].split()
                alpha_a0 = float(alpha_values[-1])
                alpha_data.append(alpha_a0)
            return alpha_data
        else:
            return None

# Read and modify XYZ file content
def modify_xyz_with_alpha(xyz_filename, alpha_data):
    modified_lines = []
    with open(xyz_filename, 'r') as xyz_file:
        num_atoms = int(xyz_file.readline().strip())
        modified_lines.append(str(num_atoms)+ f"\n")  # Write the num of atoms
        modified_lines.append(xyz_file.readline().strip()+ f"\n")  # Write the comment line
        for i, line in enumerate(xyz_file):
            if i < num_atoms:
                line = line.strip() + f"   {alpha_data[i]}    {alpha_data[i]}    {alpha_data[i]}    {0}    {0}    {0}\n"  # Append alpha(0) value to each line
            modified_lines.append(line)
    return modified_lines

# Read the number of atoms from XYZ file
def get_num_atoms_from_xyz(xyz_filename):
    with open(xyz_filename, 'r') as xyz_file:
        num_atoms = int(xyz_file.readline().strip())
        return num_atoms

# Read polarization data from "out" file
num_atoms = get_num_atoms_from_xyz("new.xyz")
alpha_data = read_alpha_from_out(num_atoms)

# Modify XYZ file content and write to a new file
modified_content = modify_xyz_with_alpha("new.xyz", alpha_data)
with open("atomic_pol.xyz", 'w') as outfile:
    for line in modified_content:
        outfile.write(line)
