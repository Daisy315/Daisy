from ase.io import read
import numpy as np

atoms = read("cluster.data", format="lammps-data")

# 动态获取所有存在的原子类型及其数量
atomic_numbers = atoms.get_atomic_numbers()
num_counts = {}
for num in atomic_numbers:
    num_counts[num] = num_counts.get(num, 0) + 1

# 确定H原子类型（数量最大的原子类型）
H_type = max(num_counts.items(), key=lambda x: (x[1], x[0]))[0]

# 获取所有非H原子的索引
C_index = np.where(atomic_numbers != H_type)[0] + 1

with open("fix_file", 'w') as writer:
    writer.write("group C id %s\nfix c_fix C setforce 0.0 0.0 0.0\n" % (' '.join(map(str, C_index))))
