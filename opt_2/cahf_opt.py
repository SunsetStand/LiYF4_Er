import numpy as np
from pyscf import gto, scf, qmmm
import cahf
import rdiis

def read_xyz(filename):
    """读取常规 xyz 文件，跳过前两行"""
    symbols = []
    coords = []
    with open(filename, 'r') as f:
        lines = f.readlines()[2:]
        for line in lines:
            if not line.strip(): continue
            parts = line.split()
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)

def read_dat(filename):
    """读取单列的 .dat 电荷文件"""
    return np.loadtxt(filename)

def run_cahf():
    # --- 1. 读取分离的坐标和电荷数据 ---
    print("正在加载团簇与点电荷文件...")
    
    # QM 层和 AIMP 层 (只有坐标)
    qm_syms, qm_coords = read_xyz('cluster.xyz')
    aimp_syms, aimp_coords = read_xyz('aimp.xyz')
    
    # PC 层 (坐标 + 电荷)
    pc_syms, pc_coords = read_xyz('rawChgs.xyz')
    pc_charges = read_dat('rawCharges.dat')
    
    # Fitted 层 (坐标 + 电荷)
    fit_syms, fit_coords = read_xyz('surfChgs.xyz')
    fit_charges = read_dat('surfaceCharges.dat')

    # 合并所有背景点电荷 (MM区)
    all_bg_coords = np.vstack([pc_coords, fit_coords])
    all_bg_charges = np.concatenate([pc_charges, fit_charges])

    # --- 2. 构建 Mole 对象 ---
    mol = gto.Mole()
    
    # --- 2. 构建 Mole 对象 ---
    mol = gto.Mole()
    
    atom_list = []
    for s, c in zip(qm_syms, qm_coords):
        atom_list.append([s, tuple(c)])
    for s, c in zip(aimp_syms, aimp_coords):
        # 核心修改：将后缀 _aimp 改为数字 99
        atom_list.append([f"{s}99", tuple(c)])
    
    mol.atom = atom_list

    # --- 3. 基组与 ECP 设定 ---
    # 定义一个“零权重”的基组：[角动量, [指数, 系数]]
    # 这里使用一个非常大的指数 (1e6) 和 0.0 的系数，使其在物理上不产生任何轨道贡献
    dummy_basis = [[0, [1e8, 1.0]]]

    mol.basis = {
        'Er': 'stuttgart_rsc',  # Er 使用 Stuttgart RSC ECP 基组
        'F': '6-31g*',
        'Li99': dummy_basis, # 使用 99 后缀匹配 AIMP 原子
        'Y99': dummy_basis,
        'F99': dummy_basis
    }
    
    mol.ecp = {
        'Er': 'stuttgart_rsc',  # Er 使用 Stuttgart RSC ECP
        'Li99': 'bfd',        # Li 和 F 属于极轻元素，使用 bfd 或 crenbl
        'F99': 'bfd',
        'Y99': 'lanl2dz'    # Y 属于重元素，使用 lanl2dz 或 Stuttgart RSC ECP
    }

    mol.charge = -5
    mol.spin = 3
    mol.verbose = 4
    mol.build()

    # 这一步极其关键，否则 AIMP 层的假电子会毁掉整个计算
    mol.nelectron = 117  # 145 (ErF8) - 28 (stuttgart_rsc) = 117

    # --- 4. 初始化 CAHF ---
    print(f"\n构建 CAHF: Er3+ 具有 7 个 4f 轨道，包含 11 个电子。")
    mf_pure = cahf.CAHF(mol, ncas=7, nelecas=11, spin=3)
    
    # 使用 qmmm 模块注入背景电荷，完美兼容 X2C
    print("注入长程背景点电荷场...")
    mf = qmmm.add_mm_charges(mf_pure, all_bg_coords, all_bg_charges)

    # --- 5. 绑定 RDIIS 与收敛参数 ---
    print("挂载 RDIIS 加速收敛，锁定 Er 的 4f 轨道...")
    # 搜索并锁定 Er 的 4f 轨道，给予 0.2 的 power 权重
    f_idx = mol.search_ao_label(['Er.*4f'])
    mf.diis = rdiis.RDIIS(rdiis_prop='dS', imp_idx=f_idx, power=0.5)

    # 检查点文件，方便后续读取轨道
    mf.chkfile = 'Er_CAHF_RDIIS.chk'

    # 对于包含重元素和大基组的体系，通常需要更好的初猜或收敛控制
    mf.init_guess = 'chk'  
    # RDIIS 的专属参数设定
    mf.diis_space = 15
    # mf.diis_start_cycle = 5        
    mf.max_cycle = 50     # 给足迭代空间
    mf.level_shift = 1.0  # 初始较大，帮助稳定收敛，后续可逐步减小
    mf.convergence_tol = 1e-8  # 先放宽标准，收敛后再以chk为初猜，加入rdiis加速收敛     
    

    # --- 5. 执行计算 ---
    print("\n--- 开始 CAHF 计算 ---")
    mf.kernel()

    if mf.converged:
        print(f"\n计算收敛！总能量: {mf.e_tot:.8f} Hartree")
    else:
        print("\n计算未收敛，可能需要调整 DIIS 或 Level Shift 策略。")

    return mf

if __name__ == "__main__":
    mf_results = run_cahf()