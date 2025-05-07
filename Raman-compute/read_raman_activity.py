import os
import math

def get_epsilon_from_OUTCAR(outcar_fh):
    """从OUTCAR文件中读取介电常数张量"""
    epsilon = []
    outcar_fh.seek(0)
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in line:
            outcar_fh.readline()  # 跳过标题行
            epsilon.append(list(map(float, outcar_fh.readline().split())))
            epsilon.append(list(map(float, outcar_fh.readline().split())))
            epsilon.append(list(map(float, outcar_fh.readline().split())))
            return epsilon
    raise ValueError("未找到介电常数张量")

def calculate_activity(base_path):
    # 常量定义
    step_size = 0.01
    vol = 3.9794829184134670 * 4.5055757901396056 * 12.7651373897190048
    coeff_map = {'1+': 0.5, '1-': -0.5}
    pi = math.pi

    # 读取norms.dat
    with open(os.path.join(base_path, 'norms.dat'), 'r') as f:
        norms = [float(line.strip()) for line in f]

    activities = []
    
    for dir_num in range(1, 88):
        dir_name = f"{dir_num:04d}"
        norm = norms[dir_num-1]
        total_ra = [[0.0]*3 for _ in range(3)]

        for polarity in ['1+', '1-']:
            # 构建文件路径
            outcar_path = os.path.join(base_path, dir_name, polarity, 'OUTCAR')
            
            # 读取介电常数
            with open(outcar_path, 'r') as f:
                epsilon = get_epsilon_from_OUTCAR(f)
            
            # 计算系数
            coeff = coeff_map[polarity]
            scaling = coeff / step_size * norm * vol / (4 * pi)
            
            # 计算当前ra并累加
            for i in range(3):
                for j in range(3):
                    total_ra[i][j] += epsilon[i][j] * scaling

        # 计算光学参数
        alpha = (total_ra[0][0] + total_ra[1][1] + total_ra[2][2]) / 3.0
        beta2 = ((total_ra[0][0] - total_ra[1][1])**2 + 
                (total_ra[0][0] - total_ra[2][2])**2 + 
                (total_ra[1][1] - total_ra[2][2])**2 + 
                6.0 * (total_ra[0][1]**2 + total_ra[0][2]**2 + total_ra[1][2]**2)) / 2.0
        activity = 45.0 * alpha**2 + 7.0 * beta2
        
        activities.append(activity)

    return activities

if __name__ == "__main__":
    base_directory = "./"  # 替换为实际路径
    results = calculate_activity(base_directory)
    
    # 输出结果
    for idx, activity in enumerate(results, 1):
        print(f"目录 {idx:04d}: 光学活性 = {activity:.6f}")
