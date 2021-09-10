# Daisy
Daisy's daily dairy
Some scripts in the process of learning molecular simulation and machine learning potential.

MTP_retrain_cfg.py 根据MTP的calc_grade功能将γ值超过范围的结构重新投入计算；
bond_length_MDAnalysis 利用MDAnalysis分析cp2k的轨迹文件并计算平均键长；
cfg2poscar 将cfg转化为poscar；
data2poscar.py 将lammps的data文件转化为vasp的Poscar文件；
get_K.py 根据poscar生成Kpoints（计算K点密度）；
json2poscar.py json转poscar；
snapshot.py 从AIMD的XDATCAR每隔帧数取结构组成Poscars；
vasp_E.pbs 能量收敛测试；
vasp_K.pbs K点测试；
vasp_disturb_on_structure.py（下载后添加拓展名.py） 对完美晶胞微扰产生新的结构；
