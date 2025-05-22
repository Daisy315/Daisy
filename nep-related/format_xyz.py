import re
import math

# ！需要修改，原子编号与符号对应关系
atom_dict = {
    '1': 'C',
    '2': 'H'
}

# 读取数据文件
with open('dumped.xyz', 'r') as f:
    data = f.readlines()

# 获取原子数
num_atoms = 0
for i in range(len(data)):
    if data[i].startswith("ITEM: NUMBER OF ATOMS"):
        num_atoms = int(data[i+1])
        break

output = []
i = 0
while i < len(data):
    if data[i].startswith("ITEM: TIMESTEP"):
        # 找到时间步
        timestep_line = i
        i += 1
        timestep = data[i].strip()
        i += 1

        # 读取原子数
        if data[i].startswith("ITEM: NUMBER OF ATOMS"):
            i += 1  # 跳过数字行
            i += 1  # ITEM: BOX BOUNDS
        else:
            raise ValueError("格式错误，找不到 'ITEM: NUMBER OF ATOMS'")

        # 提取盒子信息
        if data[i].startswith("ITEM: BOX BOUNDS xy xz yz"):
            i += 1
            xlo_bound, xhi_bound, xy = map(float, data[i].strip().split())
            i += 1
            ylo_bound, yhi_bound, xz = map(float, data[i].strip().split())
            i += 1
            zlo_bound, zhi_bound, yz = map(float, data[i].strip().split())
            i += 1
        else:
            raise ValueError("格式错误，找不到 'ITEM: BOX BOUNDS xy xz yz'")

        # 转换为真实盒子信息
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound

        # 盒子尺寸
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo

        # 晶格计算
        b2 = ly**2 + xy**2
        c2 = lz**2 + xz**2 + yz**2
        b = math.sqrt(b2)
        c = math.sqrt(c2)

        cos_alpha = (xy * xz + ly * yz) / (b * c)
        cos_beta = xz / c
        cos_gamma = xy / b

        # 基矢分量
        a_x = lx
        a_y = a_z = 0.0
        b_x = b * cos_gamma
        b_y = math.sqrt(b2 - b_x**2)
        b_z = 0.0
        c_x = c * cos_beta
        c_y = (xy * xz + ly * yz - b_x * c_x) / b_y
        c_z = math.sqrt(c**2 - c_x**2 - c_y**2)

        lattice_info = f'Lattice="{a_x} {a_y} {a_z} {b_x} {b_y} {b_z} {c_x} {c_y} {c_z}" Properties=species:S:1:pos:R:3 pol="0 0 0 0 0 0 0 0 0" energy=0 pbc="T T T"\n'

        # 读取原子信息头部
        if not data[i].startswith("ITEM: ATOMS"):
            raise ValueError("格式错误，找不到 'ITEM: ATOMS'")
        i += 1

        # 写入每帧头部信息
        output.append(f"{num_atoms}\n")
        output.append(lattice_info)

        # 替换原子编号为符号并写入
        for atom_index in range(num_atoms):
            line = data[i + atom_index].strip().split()
            if len(line) >= 4 and line[0] in atom_dict:
                line[0] = atom_dict[line[0]]
            output.append(' '.join(line) + '\n')

        i += num_atoms  # 跳过本帧的原子数据

    else:
        i += 1  # 跳过其他无关项

# 写入新文件
with open('new.xyz', 'w') as f:
    f.writelines(output)
