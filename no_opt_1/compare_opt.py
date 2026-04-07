import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
import os

def get_distances_to_center(atoms, center_index):
    """计算超胞中所有原子到中心原子的距离（考虑周期性边界条件）"""
    # mic=True 开启最小镜像约定
    distances = atoms.get_distances(center_index, range(len(atoms)), mic=True)
    return distances

def analyze_structure_change(prim_file, opt_file):
    # 1. 构建初始参考超胞 (2x2x2)
    print(f"正在构建参考超胞自: {prim_file}")
    prim = read(prim_file)
    super_init = prim * (2, 2, 2)
    
    # 按照你之前的排序逻辑：Er, Li, Y, F
    # 找到中心 Y 并替换为 Er (这里假设逻辑与之前一致)
    cell_vectors = super_init.get_cell()
    geo_center = np.sum(cell_vectors, axis=0) / 2.0
    y_indices = [i for i, a in enumerate(super_init) if a.symbol == 'Y']
    target_idx = y_indices[np.argmin([np.linalg.norm(super_init[i].position - geo_center) for i in y_indices])]
    super_init[target_idx].symbol = 'Er'
    
    # 统一元素排序以便一一对应 (Er, Li, Y, F)
    order = ['Er', 'Li', 'Y', 'F']
    new_indices = []
    for s in order:
        new_indices.extend([i for i, a in enumerate(super_init) if a.symbol == s])
    super_init = super_init[new_indices]
    
    # 2. 读取优化后的超胞
    print(f"正在读取优化后超胞: {opt_file}")
    super_opt = read(opt_file)
    
    # 检查原子数是否一致
    if len(super_init) != len(super_opt):
        print("错误：两个文件的原子总数不匹配！")
        return

    # 3. 计算距离 (中心 Er 的索引在排序后通常是 0)
    center_idx = 0 
    d_init = get_distances_to_center(super_init, center_idx)
    d_opt = get_distances_to_center(super_opt, center_idx)
    
    # 计算百分比变化
    # 排除 Er 自身 (距离为 0 的点)
    mask = d_init > 0.1
    d_init_filtered = d_init[mask]
    d_opt_filtered = d_opt[mask]
    symbols_filtered = np.array(super_init.get_chemical_symbols())[mask]
    
    diff_percent = (d_opt_filtered - d_init_filtered) / d_init_filtered * 100
    
    # 4. 按距离排序
    sort_idx = np.argsort(d_init_filtered)
    d_init_sorted = d_init_filtered[sort_idx]
    diff_sorted = diff_percent[sort_idx]
    symbols_sorted = symbols_filtered[sort_idx]

    # 5. 绘图
    plt.figure(figsize=(12, 6))
    
    # 定义颜色映射
    colors = {'F': '#30D5C8', 'Li': '#FF7F50', 'Y': '#6495ED'}
    bar_colors = [colors.get(s, 'gray') for s in symbols_sorted]
    
    bars = plt.bar(range(len(diff_sorted)), diff_sorted, color=bar_colors, alpha=0.8)
    
    # 装饰图表
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.xlabel(f"Atoms (Sorted by distance from Er, total {len(diff_sorted)} atoms)", fontsize=12)
    plt.ylabel("Distance Change Percentage (%)", fontsize=12)
    plt.title(f"Lattice Relaxation Profile: {os.path.basename(opt_file)} vs Ideal", fontsize=14)
    
    # 添加图例
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=c, lw=4, label=s) for s, c in colors.items()]
    plt.legend(handles=legend_elements, title="Elements")

    # 如果原子多，不显示 X 轴每一个刻度，改为显示距离区间
    tick_pos = np.linspace(0, len(d_init_sorted)-1, 10, dtype=int)
    plt.xticks(tick_pos, [f"{d_init_sorted[i]:.2f}Å" for i in tick_pos], rotation=45)
    
    plt.tight_layout()
    plt.savefig("relaxation_analysis.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    # 使用你上传的原始文件和优化后的文件
    analyze_structure_change('LiYF4 (1).poscar', 'POSCAR')