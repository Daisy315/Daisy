###这个脚本的作用是读取lammps的dynamical_matrix，并且计算出振动频率和简正坐标，参考 http://bbs.keinsci.com/thread-39620-1-1.html 和 https://ajz34.readthedocs.io/zh-cn/latest/QC_Notes/Freq_Series/freq_1.html#id3， 没有做空间点群的判定。后续的计算拉曼强度的脚本可参考 https://zhuanlan.zhihu.com/p/11609263666， 对坐标进行扰动后用三点差分法计算出alpha和beta2。

###注意，这里的原子质量和原子坐标是直接给出的，可以从final.data中读取，这里没有加这一步，并且需要进行单位转化（lammps的real转为高斯的单位制）

import numpy as np
from functools import partial
import scipy

np.set_printoptions(5, linewidth=150, suppress=True)
np.einsum = partial(np.einsum, optimize=["greedy", 1024 ** 3 * 2 / 8])


# In[8]:


# https://docs.scipy.org/doc/scipy/reference/constants.html
from scipy.constants import physical_constants

E_h = physical_constants["Hartree energy"][0]
a_0 = physical_constants["Bohr radius"][0]
N_A = physical_constants["Avogadro constant"][0]
c_0 = physical_constants["speed of light in vacuum"][0]
e_c = physical_constants["elementary charge"][0]
e_0 = physical_constants["electric constant"][0]
mu_0 = physical_constants["mag. constant"][0]

# 读取 dynmat.dat
with open("dynmat.dat", "r") as f:
    lines = [line.strip() for line in f if line.strip()]

# 获取原子数 Natoms
num_lines = len(lines)
Natoms = int(np.sqrt(num_lines / 3))  # 总行数 = 3*Natoms^2 → Natoms = sqrt(lines/3)
matrix_size = 3 * Natoms
D = np.zeros((matrix_size, matrix_size))

# 按顺序填充矩阵
line_idx = 0
for i in range(Natoms):          # LAMMPS 原子索引 i (1-based → 代码中 0-based)
    for alpha in range(3):  # α = x(0), y(1), z(2)
        row = 3*i + alpha   # 当前行索引
        for j in range(Natoms):  # LAMMPS 原子索引 j (1-based → 代码中 0-based)
            # 读取三个元素 (对应 β=x,y,z)
            values = list(map(float, lines[line_idx].split()))
            line_idx += 1
            for beta in range(3):
                col = 3*j + beta  # 当前列索引
                D[row, col] = values[beta]
np.save("dynmat.npy", D)
# 对称性修正：确保 D[iα,jβ] = D[jβ,iα] （消除数值误差）
D = 0.5 * ((D*4.46e-4) + ((D*4.46e-4).T))
'''
# 输出文件名
output_file = "hessian.dat"

# 写入文件
with open(output_file, "w") as f:
    # 可选：添加注释头（矩阵维度、原子数等）
    f.write(f"# Hessian Matrix (Dynamical Matrix)\n")
    f.write(f"# Dimensions: {D.shape[0]} x {D.shape[1]}\n")
    f.write(f"# Units: kcal/(mol*Å²) (assuming 'real' units in LAMMPS)\n")
    f.write(f"# Format: Each row corresponds to a 3N-dimensional row of the matrix.\n\n")
    
    # 逐行写入矩阵数据
    for row in D:
        # 将行数据转换为字符串（科学计数法，空格分隔）
        line = " ".join([f"{val:.8e}" for val in row])
        f.write(line + "\n")

print(f"Hessian matrix saved to {output_file}")
'''
mol_hess = D.reshape((Natoms, 3, Natoms, 3))
#np.save("hessian.npy", mol_hess)


data_str = "12.01115,12.01115,12.01115,12.01115,12.01115,12.01115,12.01115,12.01115,12.01115,12.01115,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797,1.00797"
mol_weight = np.array([float(x) for x in data_str.split(",")])

#data_str2 = "3.69127,3.4622,3.61841,5.25166,2.0764,3.61841,1.91556,2.36543,3.61841,3.69127,3.4622,8.44293,3.69127,3.4622,13.26744,5.25166,2.0764,8.44293,1.91556,2.36543,8.44293,5.25166,2.0764,13.26744,1.91556,2.36542,13.26744,3.69127,3.4622,18.09195,3.69127,3.4622,22.91646,5.25166,2.0764,18.09195,1.91556,2.36542,18.09195,5.25166,2.0764,22.91646,1.91556,2.36542,22.91646,3.82914,5.05204,1.20616,3.82914,5.05204,6.03067,5.60427,6.14972,1.20616,2.26802,6.43704,1.20616,5.60427,6.14972,6.03067,2.26802,6.43704,6.03067,3.82914,5.05204,10.85518,3.82914,5.05204,15.67969,5.60427,6.14972,10.85518,2.26802,6.43704,10.85518,5.60427,6.14972,15.67969,2.26802,6.43704,15.67969,3.82914,5.05204,20.5042,5.60427,6.14972,20.5042,2.26802,6.43704,20.5042"
#mol_coord = np.array([float(x) for x in data_str2.split(",")]).reshape((Natoms, 3))

#theta = np.einsum("AtBs, A, B -> AtBs", mol_hess, 1 / np.sqrt(mol_weight), 1 / np.sqrt(mol_weight)).reshape(3 * Natoms, 3 * Natoms)


# In[35]:


#print("shape of theta is :", theta.shape)


# 随后，我们对其进行对角化，可以立即得到本征值 `λ` 与简正坐标 `q`，且维度分别是 $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom})$ 与 $(3 n_\mathrm{Atom},)$。注意到 `λ` 的单位是 E<sub>h</sub> Bohr<sup>-2</sup> amu<sup>-1</sup>，而 `q` 现在是无量纲量。

# In[73]:

theta = D
λ, q = np.linalg.eigh(theta)
#print("Eigenvalue is :", λ)


# 现在获得的`λ`事实上是力常数除以质量的结果，或者按照 Levine (7ed) p63, eq (4.23) 的表达，为 $k/m$。因此，化为以波数表示的频率 `freq_cm_1` 的公式是
# 
# $$
# \tilde \nu = \frac{1}{2 \pi c_0} \sqrt{\lambda}
# $$
# 
# 其中，$c_0$ 表示真空光速。在实行具体计算前，需要将单位转换为国际单位制。最终会将频率转成 cm<sup>-1</sup> 单位。
# 对于λ为负值的情况，代入上式时应使用其绝对值，然后对频率取负值（以此代表虚频）。

# In[113]:


freq_cm_1 = np.sqrt(np.abs(λ * E_h * 1000 * N_A / a_0**2)) / (2 * np.pi * c_0 * 100) * ((λ > 0) * 2 - 1)


# In[111]:



# 需要留意，复数的频率实际上是虚数频率，或者说是现实中不存在的频率；使用复数表示这些频率仅仅是为了程序方便，以及约定俗称的原因。
# 先将3N个频率按照绝对值从小到大排序。由于该分子是非线性分子，因此其中有 6 个频率不应当归属于振动频率中，舍去。

# In[114]:


freq_cm_1=freq_cm_1.tolist()
print(freq_cm_1)
#freq_cm_1_sorted=np.array(sorted(freq_cm_1,key=lambda x:abs(x))[6:])
freq_cm_1_sorted=np.array(sorted(freq_cm_1,key=lambda x:abs(x)))



# 将剩下的3N-6个模式（对应分子内振动）按照频率数值从小到大排序，这样虚频对应的负值频率就会出现在输出的最前头，便于考察。

# In[115]:


freq_cm_1_sorted.sort()
print(freq_cm_1_sorted)

# 将每个原子质量重复三次（匹配 x, y, z 坐标）
mol_weight_expanded = np.repeat(mol_weight, 3)  #单位是amu^(1/2)

# 质量加权归一化
q_normalized = q / np.sqrt(mol_weight_expanded)[:, np.newaxis]

# 保存为 q.npy
np.save("q.npy", q_normalized)
import numpy as np

# --------------------------
# 读取原子坐标并排序 (按原子序号)
# --------------------------
atoms = []
# 假设原子坐标数据存储在"atom_coords.txt"中，格式如下：
# 3 1 1.95333 1.83212 1.91478
# 15 2 2.77906 1.09878 1.91478
# ...
with open("atom_coords.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        atom_num = int(parts[0])
        atom_type = int(parts[1])
        x, y, z = map(float, parts[2:5])
        atoms.append( (atom_num, atom_type, x, y, z) )

# 按原子序号排序 (1,2,3...30)
sorted_atoms = sorted(atoms, key=lambda x: x[0])
original_numbers = [atom[0] for atom in atoms]  # 原始顺序的原子序号

# --------------------------
# 加载计算数据
# --------------------------
q_normalized = np.load("q.npy")                    # 质量归一化特征向量 (3N×3N)

# --------------------------
# 生成VASP格式输出
# --------------------------
with open("vasp_freq_output.txt", "w") as f:
    # 写入文件头
    f.write("Eigenvectors after division by SQRT(mass)\n\n")
    f.write(" Eigenvectors and eigenvalues of the dynamical matrix\n")
    f.write(" ----------------------------------------------------\n\n")

    # 遍历每个振动模式
    for mode_idx, freq in enumerate(freq_cm_1_sorted, 1):
        # 频率行
        f.write(f"   {mode_idx} f  =   {freq:.6f} cm-1\n")
        f.write("             X         Y         Z           dx          dy          dz\n")

        # 获取当前模式的特征向量
        mode_vector = q_normalized[:, mode_idx-1]  # 假设每列对应一个模式

        # 重新排序位移向量 (按原子序号1-30)
        sorted_displacements = np.zeros(3*30)
        for i in range(30):
            orig_start = 3 * i
            target_start = 3 * (original_numbers[i] - 1)
            sorted_displacements[target_start:target_start+3] = mode_vector[orig_start:orig_start+3]

        # 写入原子坐标和位移
        for atom in sorted_atoms:
            atom_num, _, x, y, z = atom
            idx = atom_num - 1
            dx = sorted_displacements[3*idx]
            dy = sorted_displacements[3*idx+1]
            dz = sorted_displacements[3*idx+2]
            
            # 格式化输出行
            line = f"    {x:.5f}   {y:.5f}   {z:.5f}    {dx:10.6f}    {dy:10.6f}    {dz:10.6f}\n"
            f.write(line)
        
        f.write("\n")  # 模式间空行

'''
# 大多数情况下，舍去绝对值最小的六个频率即可。
# Gaussian输出的分子频率为[-1968.6121,-480.7998,183.0363,346.1818,422.7588,467.5898,607.3214,669.8050,768.4299,918.8506,974.9894,1061.4280,1268.9356,1482.8253,1651.0721]。**直接舍去前6个振动频率后计算的分子谐振频率与 Gaussian 给出的结果有少许的不同。**
# 
# 简正坐标在这里我们暂时不进行更多说明；在叙述去除平动、转动的频率后，我们再讨论简正坐标的导出。  

# ## 去除平动、转动的频率

# 去除平动、转动对频率的贡献，其过程大致是预先将平动、转动的模式求取，随后将力常数张量投影到平动、转动模式的补空间 ($3 n_\mathrm{Atom} - 6$ 维度空间)，得到新的力常数张量。
# 
# 其中的大部分内容应当在 Wilson et al.(Wilson, E. B.; Decius, J. C.; Cross, P. C. *Molecular Vibrations*; Dover Pub. Inc., 1980) 的 Chapter 2 可以找到。

# ### 质心坐标

# `center_coord` $C_t$ 表示质心坐标，维度 $(3,)$，单位 Bohr。
# 
# $$
# C_{t} = \frac{\sum_{A} A_{t} w_A}{\sum_A w_A}
# $$

# In[116]:


center_coord = (mol_coord * mol_weight[:, None]).sum(axis=0) / mol_weight.sum()



# `centered_coord` $A^\mathrm{C}_t$ 是将质心平移至原点后的原子坐标，维度 $(n_\mathrm{Atom}, 3)$，单位 Bohr。
# 
# $$
# A^\mathrm{C}_t = A_t - C_t
# $$

# In[117]:


centered_coord = mol_coord - center_coord


# ### 转动惯量本征向量

# `rot_tmp` $I_{ts}$ 是转动惯量相关的矩阵，在初始化时维度为 $(n_\mathrm{Atom}, 3, 3)$，最终结果通过求和得到 $(3, 3)$ 的矩阵，单位 Bohr<sup>2</sup> amu。
# 
# $$
# \begin{split}
# I_{ts} =
# \begin{cases}
#     \sum_{A} w_A \left( - (A_t^\mathrm{C})^2 + \sum_r (A_r^\mathrm{C})^2 \right) \,, & t = s \\
#     \sum_{A} w_A \left( - A_t^\mathrm{C} A_s^\mathrm{C} \right) \,, & t \neq s
# \end{cases}
# \end{split}
# $$

# In[118]:


rot_tmp = np.zeros((Natoms, 3, 3))
rot_tmp[:, 0, 0] = centered_coord[:, 1]**2 + centered_coord[:, 2]**2
rot_tmp[:, 1, 1] = centered_coord[:, 2]**2 + centered_coord[:, 0]**2
rot_tmp[:, 2, 2] = centered_coord[:, 0]**2 + centered_coord[:, 1]**2
rot_tmp[:, 0, 1] = rot_tmp[:, 1, 0] = - centered_coord[:, 0] * centered_coord[:, 1]
rot_tmp[:, 1, 2] = rot_tmp[:, 2, 1] = - centered_coord[:, 1] * centered_coord[:, 2]
rot_tmp[:, 2, 0] = rot_tmp[:, 0, 2] = - centered_coord[:, 2] * centered_coord[:, 0]
rot_tmp = (rot_tmp * mol_weight[:, None, None]).sum(axis=0)


# `rot_eig` $R_{ts}$ 是转动惯量相关的对称矩阵 $I_{ts}$ 所求得的本征向量，维度 $(3, 3)$，无量纲。

# In[119]:


_, rot_eig = np.linalg.eigh(rot_tmp)



# ### 平动、转动投影矩阵

# `proj_scr` $P_{A_t q}$ 是平动、转动的 $(3 n_\mathrm{Atom}, 6)$ 维度投影矩阵，其目的是将 $\Theta^{A_t B_s}$ 中不应对分子振动产生贡献的部分投影消去，剩余的 $3 n_\mathrm{Atom} - 6$ 子空间用于求取实际的分子振动频率。但在初始化 `proj_scr` $P_{A_t q}$ 时，先使用 $(n_\mathrm{Atom}, 3, 6)$ 维度的张量。
# 
# 在计算投影矩阵前，我们先生成 `rot_coord` $\mathscr{R}_{Asrw}$ 转动投影相关量，维度 $(n_\mathrm{Atom}, 3, 3, 3)$：
# 
# $$
# \mathscr{R}_{Asrw} = \sum_{t} A^\mathrm{C}_t R_{ts} R_{rw}
# $$

# In[120]:


rot_coord = np.einsum("At, ts, rw -> Asrw", centered_coord, rot_eig, rot_eig)



# 随后我们给出 `proj_scr` 的计算表达式。`proj_scr` 的前三列表示平动投影，当 $q \in (x, y, z) = (0, 1, 2)$ 时，
# 
# $$
# P_{A_t q} = \sqrt{w_A} \delta_{tq}
# $$

# 而当 $q \in (x, y, z) = (3, 4, 5)$ 时，
# 
# $$
# \begin{split}
# P_{A_t q} = \sqrt{w_A} \times
# \begin{cases}
#     \mathscr{R}_{Aytz} - \mathscr{R}_{Azty} \,, & q = x \\
#     \mathscr{R}_{Aztx} - \mathscr{R}_{Axtz} \,, & q = y \\
#     \mathscr{R}_{Axty} - \mathscr{R}_{Aytx} \,, & q = z
# \end{cases}
# \end{split}
# $$

# 最终，我们会将 $P_{A_t q}$ 中关于 $A_t$ 的维度进行归一化，因此最终获得的 $P_{A_t q}$ 是无量纲的。

# In[121]:


proj_scr = np.zeros((Natoms, 3, 6))
proj_scr[:, (0, 1, 2), (0, 1, 2)] = 1
proj_scr[:, :, 3] = (rot_coord[:, 1, :, 2] - rot_coord[:, 2, :, 1])
proj_scr[:, :, 4] = (rot_coord[:, 2, :, 0] - rot_coord[:, 0, :, 2])
proj_scr[:, :, 5] = (rot_coord[:, 0, :, 1] - rot_coord[:, 1, :, 0])
proj_scr *= np.sqrt(mol_weight)[:, None, None]
proj_scr.shape = (-1, 6)
proj_scr /= np.linalg.norm(proj_scr, axis=0)



# 最后我们声明，在经过上述投影后的力常数矩阵几乎表现为零：
# 
# $$
# \mathbf{P}^\dagger \mathbf{\Theta} \mathbf{P} \simeq \mathbf{0}
# $$

# In[122]:


proj_scr.T @ theta @ proj_scr


# 对上述矩阵进行对角化所给出的平动、转动频率如下：

# In[123]:


e_tr, _ = np.linalg.eigh(proj_scr.T @ theta @ proj_scr)
np.sqrt(np.abs(e_tr * E_h * 1000 * N_A / a_0**2)) / (2 * np.pi * c_0 * 100) * ((e_tr > 0) * 2 - 1)


# ### 平动、转动投影矩阵的补空间

# 既然我们已经得到了平动、转动的投影，那么根据矩阵的原理，相应地我们也能获得其补空间的投影。我们令 `proj_inv` $Q_{A_t q}$ 为 $P_{A_t q}$ 的补空间投影。获得补空间的大致方式是预先定义一个仅有一个分量为 $1$ 的 $(3 n_\mathrm{Atom}, )$ 维度向量，随后通过 Schmit 正交的方式给出已有投影空间的补空间向量。组合这些 Schmit 正交的向量便获得了 $Q_{A_t q}$。
# 
# $Q_{A_t q}$ 的维度本应当是 $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom} - 6)$ 维。但为了程序编写方便，我们先规定 `proj_inv` 是 $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom})$ 维度，并且其中的前 6 列填入 $P_{A_t q}$；在进行 Schmit 正交化后，再将前 6 列剔除。

# In[124]:


proj_inv = np.zeros((Natoms * 3, Natoms * 3))
proj_inv[:, :6] = proj_scr
cur = 6
for i in range(0, Natoms * 3):
    vec_i = np.einsum("Ai, i -> A", proj_inv[:, :cur], proj_inv[i, :cur])
    vec_i[i] -= 1
    if np.linalg.norm(vec_i) > 1e-8:
        proj_inv[:, cur] = vec_i / np.linalg.norm(vec_i)
        cur += 1
    if cur >= Natoms * 3:
        break
proj_inv = proj_inv[:, 6:]


# 我们最后获得的 $Q_{A_t q}$ 是列正交切归一的矩阵，且形式大致是下三角矩阵。但需要留意，对于当前的分子，最后一列只有 6 个非零值，与倒数第二列非零值的数量相差 2 个。

# In[125]:


proj_inv[:, :8]


# In[126]:


proj_inv[:, 8:]


# ### 去除平动、转动部分的频率

# 我们将对矩阵 $\mathbf{Q}^\dagger \mathbf{\Theta} \mathbf{Q}$ 进行对角化；且获得的第 $q$ 个简正坐标的频率相关量 `e` $K_q = k_q / m_q$ 与原始简正坐标 `q` $\mathbf{q}^\mathrm{orig}$ 表示如下：
# 
# $$
# \mathbf{Q}^\dagger \mathbf{\Theta} \mathbf{Q} \mathbf{q}^\mathrm{orig} = \mathbf{q}^\mathrm{orig} \mathrm{diag} (\boldsymbol{K})
# $$

# In[127]:


e, q = np.linalg.eigh(proj_inv.T @ theta @ proj_inv)


# 由此，我们就可以立即获得去除平动、转动部分的，以 cm<sup>-1</sup> 为单位的，总数为 $3 n_\mathrm{Atom} - 6$ 的分子频率 `freq_cm_1`：

# In[128]:


freq_cm_1 = np.sqrt(np.abs(e * E_h * 1000 * N_A / a_0**2)) / (2 * np.pi * c_0 * 100) * ((e > 0) * 2 - 1)

print(freq_cm_1.shape)

# 下面的计算结果显示与Gaussian输出的分子频率之间的误差非常小。

# In[130]:

print("Calculated frequency here:",freq_cm_1)

# ### 归一化的简正坐标

# 方才通过对角化，我们获得的原始简正坐标 `q` 的维度是 $3 n_\mathrm{Atom} - 6$。我们需要通过 `q` $\mathbf{q}^\mathrm{orig}$ 重塑回正常的简正坐标的维度 $q_{A_t q}$ $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom} - 6)$。
# 
# 我们首先给出未经过归一化的简正坐标，命名为 `q_unnormed` $q_{A_t q}^\mathrm{unnorm}$，其单位是 amu<sup>-1/2</sup>。该量将会用于后续的红外强度计算上。其计算过程大致是
# 
# $$
# \mathbf{q}^\mathrm{unnorm} = \mathbf{Q} \mathbf{q}^\mathrm{orig} / \sqrt{\mathbf{w}}
# $$

# In[131]:


q_unnormed = np.einsum("AtQ, A -> AtQ", (proj_inv @ q).reshape(Natoms, 3, (proj_inv @ q).shape[-1]), 1 / np.sqrt(mol_weight))
q_unnormed = q_unnormed.reshape(-1, q_unnormed.shape[-1])


# 而将每一个简正坐标的振动强度归一化的矩阵称为 `q_normed` $q_{A_t q}$；它是我们的目标的简正坐标。

# In[132]:


q_normed = q_unnormed / np.linalg.norm(q_unnormed, axis=0)


# 我们可以以下述代码核对前三个简正坐标。这些坐标应当与 Gaussian 所输出的坐标几乎相同，或刚好相差正负号。

# In[133]:


q_normed.reshape(Natoms, 3, 3 * Natoms - 6)[:, :, :3].transpose((2, 0, 1))

np.save('q_unnormed.npy', q_unnormed)
np.save('q_normed.npy', q_normed)
# ## 修订记录

# - 2021-06-21：重新写了 Schmidt 正交化代码。我不太理解当时的代码到底为什么是对的 (>.<)

# In[ ]:


'''
