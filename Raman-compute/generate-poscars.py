#generate poscars from phonons for FDM to calculate Raman activity

import re
import math

def parse_outeigenvectors(outcar_file):
    modes = []
    current_mode = None
    found_section = False

    with open(outcar_file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "Eigenvectors after division by SQRT(mass)" in line:
            found_section = True
            continue
        
        if found_section:
            mode_match = re.match(r'^\s*(\d+)\s+f\s*=', line)
            if mode_match:
                # 开始新振动模式
                if current_mode is not None:
                    modes.append(current_mode)
                
                mode_number = int(mode_match.group(1))
                parts = line.split()
                cm_index = parts.index('cm-1') - 1
                cm = float(parts[cm_index])
                
                current_mode = {
                    'number': mode_number,
                    'cm': cm,
                    'atoms': []
                }
                
                # 寻找坐标行
                for j in range(idx+1, len(lines)):
                    if 'X         Y         Z           dx          dy          dz' in lines[j]:
                        start_idx = j + 1
                        break
                else:
                    start_idx = idx + 1
                
                # 收集原子数据
                atoms = []
                for k in range(start_idx, len(lines)):
                    atom_line = lines[k].strip()
                    if not atom_line:
                        continue
                    if re.match(r'^\s*\d+\s+f\s*=', atom_line):
                        break
                    
                    parts = atom_line.split()
                    if len(parts) < 6:
                        break
                    
                    try:
                        x, y, z = map(float, parts[0:3])
                        dx, dy, dz = map(float, parts[3:6])
                        atoms.append((x, y, z, dx, dy, dz))
                    except:
                        break
                
                current_mode['atoms'] = atoms
    
    if current_mode is not None:
        modes.append(current_mode)
    
    return modes

def write_poscars(modes, step_size=0.01, disps=[-1, 1]):
    poscar_header = """POSCAR file written by OVITO Basic 3.9.1
   1.00000000000000     
     3.9794829184134670    0.0000000000000000    0.0000000000000000
     0.0000000000000000    4.5055757901396056    0.0000000000000000
     0.0000000000000000    0.0000000000000000   12.7651373897190048
   C    H 
    10    20
Cartesian
"""

    for mode in modes:
        atoms = mode['atoms']
        if not atoms:
            continue
        
        # 计算归一化因子
        sum_sq = sum(dx**2 + dy**2 + dz**2 for (_, _, _, dx, dy, dz) in atoms)
        norm = math.sqrt(sum_sq)
        print(norm)
        if norm == 0:
            continue
        
        for disp in disps:
            scaled_disp = step_size * disp / norm
            new_coords = []
            for (x, y, z, dx, dy, dz) in atoms:
                new_x = x + dx * scaled_disp
                new_y = y + dy * scaled_disp
                new_z = z + dz * scaled_disp
                new_coords.append((new_x, new_y, new_z))
            
            # 生成文件名
            filename = f"POSCAR.{mode['number']:04d}.{'+1' if disp > 0 else '-1'}"
            
            # 写入文件
            with open(filename, 'w') as f:
                f.write(poscar_header)
                for coord in new_coords:
                    f.write(f"  {coord[0]:16.8f}  {coord[1]:16.8f}  {coord[2]:16.8f}\n")

if __name__ == "__main__":
    modes = parse_outeigenvectors("OUTCAR")
    write_poscars(modes)
