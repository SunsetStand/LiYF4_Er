import numpy as np
import matplotlib.pyplot as plt
from pyscf import lib
from pyscf.lib import chkfile

def plot_4f_splitting(chk_path, metal_label='Er'):
    # 1. 从 chk 文件加载数据
    if not lib.os.path.exists(chk_path):
        print(f"Error: 找不到文件 {chk_path}")
        return

    data = chkfile.load(chk_path, 'scf')
    mo_energy = data['mo_energy']
    mo_occ = data['mo_occ']
    
    # 单位转换常量
    Hartree2cm = 219474.63
    
    # 2. 识别 4f 轨道
    # 在 ROHF/CAHF 中，单占轨道（singly occupied）通常对应 4f
    # Er3+ 有 11 个 4f 电子，在 ROHF 下通常表现为 3 个未成对电子 (spin=3)
    # 这里我们提取所有单占轨道和能量最接近的占据/虚拟轨道
    f_indices = np.where((mo_occ > 0) & (mo_occ < 2))[0]
    
    # 如果是 CAHF 锁定后的轨道，我们需要寻找能量在特定范围内的 7 个 4f 轨道
    # 这里的逻辑是寻找 HOMO 附近的轨道能级
    all_f_energies = mo_energy[f_indices]
    
    if len(all_f_energies) == 0:
        print("未检测到单占轨道，请检查 chk 文件是否为 ROHF/CAHF 结果。")
        return

    # 以最低的 4f 轨道为零点
    relative_energies_cm = (all_f_energies - np.min(all_f_energies)) * Hartree2cm
    relative_energies_cm = np.sort(relative_energies_cm)

    # 3. 绘图部分
    plt.figure(figsize=(5, 8))
    x_position = np.ones_like(relative_energies_cm)
    
    # 绘制能级横线
    for i, e in enumerate(relative_energies_cm):
        plt.hlines(e, 0.7, 1.3, colors='blue', linewidth=2)
        plt.text(1.35, e, f"{e:.1f} $cm^{{-1}}$", verticalalignment='center')

    plt.xlim(0.5, 2.0)
    plt.xticks([])
    plt.ylabel("Relative Energy ($cm^{-1}$)", fontsize=12)
    plt.title(f"Stark Splitting of {metal_label} 4f Orbitals\n(from QMMM-CAHF)", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('4f_splitting.png', dpi=300)
    plt.show()

    print("--- 4f 能级数据 (cm^-1) ---")
    for i, e in enumerate(relative_energies_cm):
        print(f"Level {i+1}: {e:10.2f} cm^-1")

# 使用方法
plot_4f_splitting('/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/Er_CAHF_RDIIS.chk', metal_label='Er3+')