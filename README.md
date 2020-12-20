# Daisy
Daisy's daily dairy
Some scripts in the process of learning molecular simulation and machine learning potential.

MTP_retrain_cfg.py 根据MTP的calc_grade功能将γ值超过范围的结构重新投入计算；
bond_length_MDAnalysis 利用MDAnalysis分析cp2k的轨迹文件并计算平均键长；
cfg2poscar 将cfg转化为poscar；
get_K.py 根据poscar生成Kpoints；
json2poscar.py json转poscar；
snapshot.py 从AIMD的XDATCAR每隔帧数取结构组成Poscars；
vasp_E.pbs 能量收敛测试；
vasp_K.pbs K点测试（K点密度不变）；
vasp_disturb_on_structure.py 对完美晶胞微扰产生新的结构；
