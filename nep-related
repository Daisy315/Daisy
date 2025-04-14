import os
from pylab import *
from ase.io import read
from sklearn.decomposition import PCA
from calorine.nep import get_descriptors, \
                         get_potential_forces_and_virials, \
                         get_dipole, \
                         get_latent_space


### 设置matplotlib属性 ###
aw = 2
fs = 18
lw = 2
font = {'size': fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

# set_fig_properties(): 设置图像属性 
def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)

### 创建结果目录 ###
os.makedirs('./results-fit-current', exist_ok=True)

### 读取训练和测试数据 ### 
def extract(xyz_path):
    structures = read(xyz_path, ":")
    # structures = [structure for structure in structures 
                  # if "Fe" not in structure.get_chemical_symbols() 
                  # and "Co" not in structure.get_chemical_symbols()]
    return structures

train = extract("train.xyz")
test = extract("added_sphere.xyz")

print(len(train))
print(len(test))

### 统计原子数量 ###
# 遍历结构,记录C、H原子索引及总量
carbon_idxs_train = []
atom_offset = 0
for atoms in train:
    c_idxs_train = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s=='C']
    carbon_idxs_train.extend([idx + atom_offset for idx in c_idxs_train])  
    atom_offset += len(atoms)

hydrogen_idxs_train = []
atom_offset = 0
for atoms in train:
    h_idxs_train = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s=='H']
    hydrogen_idxs_train.extend([idx + atom_offset for idx in h_idxs_train])  
    atom_offset += len(atoms)

carbon_idxs_test = []
atom_offset = 0
for atoms in test:
    c_idxs_test = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s=='C']
    carbon_idxs_test.extend([idx + atom_offset for idx in c_idxs_test])  
    atom_offset += len(atoms)

hydrogen_idxs_test = []
atom_offset = 0
for atoms in test:
    h_idxs_test = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s=='H']
    hydrogen_idxs_test.extend([idx + atom_offset for idx in h_idxs_test])  
    atom_offset += len(atoms)

###分别记录训练集和测试集中的原子总数、碳原子数和氢原子数###
str_list = np.array([sum([len(structure) for structure in train]), 
                    sum([len(structure) for structure in test])])
str_sum = np.cumsum(str_list)
str_sum = str_sum.astype(int)

c_nums = np.array([sum([atoms.get_chemical_symbols().count('C') 
                        for atoms in train]),
                   sum([atoms.get_chemical_symbols().count('C')   
                        for atoms in test])])
c_sums = np.cumsum(c_nums)
c_sums = c_sums.astype(int)  

h_nums = np.array([sum([atoms.get_chemical_symbols().count('H') 
                        for atoms in train]),
                   sum([atoms.get_chemical_symbols().count('H')   
                        for atoms in test])])
h_sums = np.cumsum(h_nums)
h_sums = h_sums.astype(int) 
with open("./results-fit-current/str_sum.txt", 'a') as f:
    np.savetxt(f, str_sum, fmt='%d') 
    np.savetxt(f, c_sums, fmt='%d')
    np.savetxt(f, h_sums, fmt='%d')

###投影###
descriptors_train = []
descriptors_test = []
latent = []
i = 0
for subset in [train]:
    for structure in subset:
        descriptors_train.append(get_descriptors(structure, model_filename='nep.txt'))
        latent.append(get_latent_space(structure, model_filename='nep.txt'))  #list,每个元素有atom_num行，每行有n_neu列，但是同类元素的latent_value基本类似
        if i % 100 == 0:
            print(f"Now we have addressed {i} structures.\n")
        i += 1

for subset in [test]:
    for structure in subset:
        descriptors_test.append(get_descriptors(structure, model_filename='nep.txt'))

# Concatenate all descriptors
all_descriptors_train = np.concatenate(descriptors_train, axis=0)
all_descriptors_test = np.concatenate(descriptors_test, axis=0)
all_latent = np.concatenate(latent, axis=0)    #转化为数组,shape = (total_atom_num * n_neu)
print(f'Shape of descriptors_train: {all_descriptors_train.shape}')
print(f'Shape of descriptors_test: {all_descriptors_test.shape}')
print(f'Total number of atoms in training dataset: {all_descriptors_train.shape[0]}')
print(f'Total number of atoms in test dataset: {all_descriptors_test.shape[0]}')
print(f'Number of descriptor components:  {all_descriptors_test.shape[1]}')
print(f'Number of latent components:      {all_latent.shape[1]}')
pca = PCA(n_components=2)

pc = pca.fit(all_descriptors_test)  #descriptor_space
pc_train = pca.transform(all_descriptors_train)
pc_test = pca.transform(all_descriptors_test)
np.savetxt('./results-fit-current/pc_train.txt', pc_train)
np.savetxt("./results-fit-current/pc_test.txt", pc_test)

carbon_pcs_train = pc_train[carbon_idxs_train]  #descriptor_space for carbon
np.savetxt("./results-fit-current/pc_train_c.txt", carbon_pcs_train)
carbon_pcs_test = pc_test[carbon_idxs_test]  #latent_space for carbon
np.savetxt("./results-fit-current/pc_test_c.txt", carbon_pcs_test)

hydrogen_pcs_train = pc_train[hydrogen_idxs_train]  #descriptor_space for hydrogen
np.savetxt("./results-fit-current/pc_train_h.txt", hydrogen_pcs_train)
hydrogen_pcs_test = pc_test[hydrogen_idxs_test]    #latent_space for hydrogen
np.savetxt("./results-fit-current/pc_test_h.txt", hydrogen_pcs_test)

# descriptor_mean = all_descriptors.mean(axis=0)
# descriptor_std = all_descriptors.std(axis=0)
# normalized_descriptors = [(d - descriptor_mean) / descriptor_std for d in descriptors]
# # Perform PCA
# pca = PCA(n_components=2)
# pc = pca.fit_transform(np.concatenate(normalized_descriptors, axis=0))
# np.savetxt("pc_normalized.txt", pc)


###画图###
figure(figsize=(10, 10))
set_fig_properties([gca()])
scatter(pc_train[:, 0], pc_train[:, 1], label='Training dataset', alpha = 0.5, zorder = 1)
scatter(pc_test[:, 0], pc_test[:, 1], label='Testing dataset', alpha = 0.5, zorder = 2)
legend(frameon = False)
xlabel('PC1')
ylabel('PC2')
savefig('./results-fit-current/Descriptor_space.png', dpi=300, bbox_inches='tight')

figure(figsize=(10, 10))
set_fig_properties([gca()])
scatter(carbon_pcs_train[:, 0], carbon_pcs_train[:, 1], label='Training dataset', alpha = 0.5, zorder = 1)
scatter(carbon_pcs_test[:, 0], carbon_pcs_test[:, 1], label='Testing dataset', alpha = 0.5, zorder = 2)
legend(frameon = False)
xlabel('PC1')
ylabel('PC2')
savefig('./results-fit-current/Descriptor_space_c.png', dpi=300, bbox_inches='tight')

figure(figsize=(10, 10))
set_fig_properties([gca()])
scatter(hydrogen_pcs_train[:, 0], hydrogen_pcs_train[:, 1], label='Training dataset', alpha = 0.5, zorder = 1)
scatter(hydrogen_pcs_test[:, 0], hydrogen_pcs_test[:, 1], label='Testing dataset', alpha = 0.5, zorder = 2)
legend(frameon = False)
xlabel('PC1')
ylabel('PC2')
savefig('./results-fit-current/Descriptor_space_h.png', dpi=300, bbox_inches='tight')
