#根据起始结构（init.data）判定碳链，计算MD轨迹的setting angle
import numpy as np
import sys

def read_lammps(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # 找到原子类型信息
    masses_start_index = lines.index("Masses\n") + 1
    carbon_type = hydrogen_type = None
    
    # 读取质量信息
    for i in range(masses_start_index+1, masses_start_index + 3):  
        parts = lines[i].split()
        atom_type = int(parts[0])
        mass = float(parts[1])
        
        if mass == 12.011:
            carbon_type = atom_type
        elif mass == 1.00794:
            hydrogen_type = atom_type

    # 找到原子数据开始的位置
    atoms_start_index = lines.index("Atoms # full\n") + 2
    
    positions = []
    types = []
    atom_indices = []  # 添加这个列表来存储原子的序号
    
    for i in range(atoms_start_index, len(lines)):
        if i >= 6340:  # 防止超出读取限制
            break
        parts = lines[i].split()
        if len(parts) < 7:  # 确保数据完整
            continue
        atom_index = int(parts[0])  # 原子序号为第一个数字
        atom_type = int(parts[2])  # 原子类型
        x, y, z = map(float, parts[4:7])  # 读取坐标
        types.append(atom_type)
        positions.append([x, y, z])
        atom_indices.append(atom_index)  # 存储原子序号
    
    return np.array(types), np.array(positions), np.array(atom_indices), carbon_type

def find_chains(types, positions, carbon_type, threshold=1.6):
    carbon_indices = [i for i, t in enumerate(types) if t == carbon_type]
    chains = []
    visited = set()
    
    for c_index in carbon_indices:
        if c_index in visited:
            continue

        chain = []
        to_visit = [c_index]
        
        while to_visit:
            current = to_visit.pop()
            if current in visited:
                continue
            
            visited.add(current)
            chain.append(current)

            # 检查邻近的碳原子
            for other in carbon_indices:
                if other not in visited and np.linalg.norm(positions[current] - positions[other]) <= threshold:
                    to_visit.append(other)

        if len(chain) >= 2:  # 至少需要2个碳原子
            chains.append(chain)
    
    return chains

def compute_setting_angles(chain, positions, tilt_threshold=18):
    setting_angles = []
    for i in range(len(chain) - 2):
        c1, c2, c3 = chain[i], chain[i + 1], chain[i + 2]
        
        # 计算C1和C3的中点O
        midpoint_O = (positions[c1] + positions[c3]) / 2
        
        # 从O到C2的向量
        vec_O_C2 = positions[c2] - midpoint_O
        
        # 计算与x轴的夹角
        angle_with_x = np.degrees(np.arctan2(vec_O_C2[1], vec_O_C2[0]))

        # 归一化角度
        angle_with_x = angle_with_x % 360
        if angle_with_x > 180:
            angle_with_x = 360 - angle_with_x
            
        # 计算 C1-C3 的向量，更改为与 z 轴的夹角计算
        vec_C1_C3 = positions[c3] - positions[c1]
        vec_C1_C3_normalized = vec_C1_C3 / np.linalg.norm(vec_C1_C3)  # 归一化以计算夹角

        # 计算 C1-C3 向量与 z 轴的夹角
        z_axis = np.array([0, 0, 1])  # z 轴上的单位向量
        cos_angle = np.dot(vec_C1_C3_normalized, z_axis)
        angle_C1_C3_z = np.degrees(np.arccos(cos_angle))

        # 进行验证
        if angle_C1_C3_z <= tilt_threshold:
            setting_angles.append(angle_with_x)  # 仅当夹角小于等于阈值时保存角度

    return setting_angles

def read_and_extract_atoms(file_path, indices):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # 找到原子数据开始的位置
    atoms_start_index = lines.index("Atoms # full\n") + 2
    
    atom_coordinates = {}
    
    for i in range(atoms_start_index, len(lines)):
        if i >= 6340:  # 防止超出读取限制
            break
        parts = lines[i].split()
        if len(parts) < 7:  # 确保数据完整
            continue
        atom_index = int(parts[0])  # 原子序号为第一个数字
        if atom_index in indices:
            x, y, z = map(float, parts[4:7])  # 读取坐标
            atom_coordinates[atom_index] = (x, y, z)
    
    return atom_coordinates

def write_angles_to_file(file_path, angles):
    with open(file_path, 'w') as f:
        for angle in angles:
            f.write(f"{angle}\n")

def main(file_path, external_file):
    types, positions, atom_indices, carbon_type = read_lammps(file_path)
    chains = find_chains(types, positions, carbon_type)

    valid_chains = []
    setting_angles_init = []
    for i, chain in enumerate(chains):
        atom_chain_indices = [atom_indices[c] for c in chain]  # 获取链的原子序号
        if len(chain) == 21:  # 检查链是否确切包含21个碳原子
            angles = compute_setting_angles(chain, positions)
            valid_chains.append((i + 1, atom_chain_indices, angles))  # 保存链的序号、原子序号和角度
            setting_angles_init.extend(angles)  # 收集所有初始的设定角
            print(f"Chain {i + 1} Atom Indices: {atom_chain_indices}, Setting Angles: {angles}")

    print(f"Total Carbon Chains Identified: {len(valid_chains)}")
    if len(valid_chains) == 96:
        print("Successfully identified the expected number of chains (96).")
    else:
        print(f"Warning: Expected 96 chains, but found {len(valid_chains)}.")
        
    # 将设定角写入到文件
    write_angles_to_file("setting_angles_init.data", setting_angles_init)

    setting_angles_external = []
    # 从外部文件提取坐标并重新计算设定角
    for chain_info in valid_chains:
        chain_index, atom_indices, _ = chain_info  # 获取链的序号和原子序号
        atom_coords = read_and_extract_atoms(external_file, atom_indices)

        # 将提取的坐标转换为NumPy数组，便于计算
        positions_external = np.array([atom_coords[index] for index in atom_indices if index in atom_coords])
        
        # 计算设定角
        angles_external = compute_setting_angles(range(len(positions_external)), positions_external)
        setting_angles_external.extend(angles_external)  # 收集所有外部计算的设定角

        print(f"Chain {chain_index} Setting Angles (from external): {angles_external}")
        
    # 将外部设定角写入到文件
    write_angles_to_file("setting_angles_traj.data", setting_angles_external)

# 示例用法：
file_path = 'init.data'  # 输入文件路径
external_file = sys.argv[1]  # 其他结构文件路径
main(file_path, external_file)
