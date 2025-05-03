#!/usr/bin/env python
# coding: utf-8

# # 频率分析得到分子频率与简正模式



# 分子对应的输入文件`C2O4H.gjf`、输出文件`C2O4H.log` 与 fchk 文件`C2O4H.fchk` 在`files`文件夹中。这份文档的目标将是重复输出文件中的分子频率 `Frequencies` (单位 cm<sup>-1</sup>) 与简正坐标部分；以下是其中一部分频率分析的输出：

# In[6]:


# with open("files/C2O4H.log", "r") as f:
    # while "and normal coordinates" not in f.readline(): continue
    # for _ in range(14): print(f.readline()[:-1])


# >频率分析 (1) 文档的目的与卢天 (Sobereva) 的 `Hess2freq` [程序](http://sobereva.com/328) 的程序基本相同，文档的编写过程也受到不少启发。

# ## 环境准备

# 下述引入的包中，
# 
# * `FormchkInterface` 可以用来读取 fchk 文件的信息；文件出自 [pyxdh](https://github.com/ajz34/Py_xDH/tree/master) 项目。
# 
# * 文档中我们可能会使用众多物理常数。这些由 SciPy 提供，数据来源是 CODATA 2014。

# In[7]:


from files.formchk_interface import FormchkInterface
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


# 现在我们准备分子的数据：
# 
# * `natm` 原子数量 $n_\mathrm{Atom}$
# 
# * `mol_weight` 原子质量 $w_A$，向量长度 $(n_\mathrm{Atom},)$，单位 amu
# 
# * `mol_coord` 原子坐标 $A_t$，矩阵大小 $(n_\mathrm{Atom}, 3)$，单位 Bohr
# 
# * `mol_hess` 坐标二阶梯度 (Hessian 矩阵) $E_\mathrm{tot}^{A_t B_s}$，张量大小 $(n_\mathrm{Atom}, 3, n_\mathrm{Atom}, 3)$，单位 E<sub>h</sub> Bohr<sup>-2</sup>
# 
# 本文档中，$A, B$ 指代原子，$t, s$ 指代坐标分量 $x, y$ 或 $z$。

# In[10]:


fchk = FormchkInterface("files/C2O4H.fchk")


# In[13]:


mol_weight = fchk.key_to_value("Real atomic weights")
natm = mol_weight.size
mol_coord = fchk.key_to_value("Current cartesian coordinates").reshape((natm, 3))
mol_hess = fchk.hessian()  #下三角矩阵元素（储存需要）转为对称矩阵
mol_hess = (mol_hess + mol_hess.T) / 2 #满足严格对称性
print(mol_hess)
mol_hess = mol_hess.reshape((natm, 3, natm, 3)) #四维数组，转为NMs个3*natoms的数组
print(mol_hess)

# In[110]:


mol_hess.shape


# ## 包含平动、转动的频率

# 这里按照 `Hess2freq` 程序的思路进行叙述。我们首先生成带原子质量权重的力常数张量 `theta`
# 
# $$
# \Theta^{A_t B_s} = E_\mathrm{tot}^{A_t B_s} / \sqrt{w_A w_B}
# $$
# 
# 但为了程序便利，我们重定义 `theta` 的维度信息为 $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom})$；单位是 E<sub>h</sub> Bohr<sup>-2</sup> amu<sup>-1</sup>。

# In[19]:


theta = np.einsum("AtBs, A, B -> AtBs", mol_hess, 1 / np.sqrt(mol_weight), 1 / np.sqrt(mol_weight)).reshape(3 * natm, 3 * natm)


# In[35]:


theta.shape


# 随后，我们对其进行对角化，可以立即得到本征值 `λ` 与简正坐标 `q`，且维度分别是 $(3 n_\mathrm{Atom}, 3 n_\mathrm{Atom})$ 与 $(3 n_\mathrm{Atom},)$。注意到 `λ` 的单位是 E<sub>h</sub> Bohr<sup>-2</sup> amu<sup>-1</sup>，而 `q` 现在是无量纲量。

# In[73]:


λ, q = np.linalg.eigh(theta)
λ


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
freq_cm_1


# 接下来获取分子点群，判断是不是线性分子,决定转动和平动自由度 '''added'''

# In[111]:


from files.point_group_detection import *
element_dict = dict(zip(NUC.values(), NUC.keys()))
mol_coord = fchk.key_to_value("Current cartesian coordinates").reshape((natm, 3))
atomic_list = [element_dict[i] for i in fchk.key_to_value("Atomic numbers")]
gpname,_,_ =detect_symm(list(zip(atomic_list,mol_coord.tolist()))) 
if gpname=="Dooh" or gpname=="Coov":
    print("转动和平动自由度是5.")
else:
    print("转动和平动自由度是6.")


# 需要留意，复数的频率实际上是虚数频率，或者说是现实中不存在的频率；使用复数表示这些频率仅仅是为了程序方便，以及约定俗称的原因。
# 先将3N个频率按照绝对值从小到大排序。由于该分子是非线性分子，因此其中有 6 个频率不应当归属于振动频率中，舍去。

# In[114]:


freq_cm_1=freq_cm_1.tolist()
freq_cm_1_sorted=np.array(sorted(freq_cm_1,key=lambda x:abs(x))[6:])
freq_cm_1_sorted


# 将剩下的3N-6个模式（对应分子内振动）按照频率数值从小到大排序，这样虚频对应的负值频率就会出现在输出的最前头，便于考察。

# In[115]:


freq_cm_1_sorted.sort()
freq_cm_1_sorted


# 大多数情况下，舍去绝对值最小的六个频率即可。
# Gaussian输出的分子频率为[-1968.6121,-480.7998,183.0363,346.1818,422.7588,467.5898,607.3214,669.8050,768.4299,918.8506,974.9894,1061.4280,1268.9356,1482.8253,1651.0721]。**直接舍去前6个振动频率后计算的分子谐振频率与 Gaussian 给出的结果有少许的不同。**
# 
# 简正坐标在这里我们暂时不进行更多说明；在叙述去除平动、转动的频率后，我们再讨论简正坐标的导出。  '''added'''

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
center_coord


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


rot_tmp = np.zeros((natm, 3, 3))
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
rot_eig


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
rot_coord.shape


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


proj_scr = np.zeros((natm, 3, 6))
proj_scr[:, (0, 1, 2), (0, 1, 2)] = 1
proj_scr[:, :, 3] = (rot_coord[:, 1, :, 2] - rot_coord[:, 2, :, 1])
proj_scr[:, :, 4] = (rot_coord[:, 2, :, 0] - rot_coord[:, 0, :, 2])
proj_scr[:, :, 5] = (rot_coord[:, 0, :, 1] - rot_coord[:, 1, :, 0])
proj_scr *= np.sqrt(mol_weight)[:, None, None]
proj_scr.shape = (-1, 6)
proj_scr /= np.linalg.norm(proj_scr, axis=0)
proj_scr


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


proj_inv = np.zeros((natm * 3, natm * 3))
proj_inv[:, :6] = proj_scr
cur = 6
for i in range(0, natm * 3):
    vec_i = np.einsum("Ai, i -> A", proj_inv[:, :cur], proj_inv[i, :cur])
    vec_i[i] -= 1
    if np.linalg.norm(vec_i) > 1e-8:
        proj_inv[:, cur] = vec_i / np.linalg.norm(vec_i)
        cur += 1
    if cur >= natm * 3:
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
freq_cm_1


# 下面的计算结果显示与Gaussian输出的分子频率之间的误差非常小。

# In[130]:


gaussian_freq = np.array([-1968.6121,-480.7998,183.0363,346.1818,422.7588,467.5898,607.3214,669.8050,768.4299,918.8506,974.9894,1061.4280,1268.9356,1482.8253,1651.0721])
print("Gaussian frequency:",gaussian_freq)
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


q_unnormed = np.einsum("AtQ, A -> AtQ", (proj_inv @ q).reshape(natm, 3, (proj_inv @ q).shape[-1]), 1 / np.sqrt(mol_weight))
q_unnormed = q_unnormed.reshape(-1, q_unnormed.shape[-1])
print(q_unnormed)

# 而将每一个简正坐标的振动强度归一化的矩阵称为 `q_normed` $q_{A_t q}$；它是我们的目标的简正坐标。

# In[132]:


q_normed = q_unnormed / np.linalg.norm(q_unnormed, axis=0)


# 我们可以以下述代码核对前三个简正坐标。这些坐标应当与 Gaussian 所输出的坐标几乎相同，或刚好相差正负号。

# In[133]:


q_normed.reshape(natm, 3, 3 * natm - 6)[:, :, :3].transpose((2, 0, 1))

print(q_normed)
# ## 修订记录

# - 2021-06-21：重新写了 Schmidt 正交化代码。我不太理解当时的代码到底为什么是对的 (>.<)

# In[ ]:




