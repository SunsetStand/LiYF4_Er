import numpy as np
from pyscf import gto, scf, mcscf, lib
from pyscf.lib import chkfile
from src.AIMP3_DMET_SCEI import AIMPEnvLoader, AIMP_ROHF # 确保路径正确
# 如果你的环境里有师兄用的 sacasscf_mixer，也可以直接调用

def run_casci_from_chk(chk_path, input_yaml='input.yaml'):
    # 1. 环境与分子重建 (必须与原计算保持一致)
    # 建议直接读取你原有的配置文件来重建 CLUS_MOL，这里为示例逻辑
    # mol = ... (重建你的 Er 配合物分子对象)
    
    # 2. 加载 MF 对象并读入 chk 数据
    # 假设你使用的是 AIMP_ROHF
    # mf = AIMP_ROHF(mol, AIMP_LOADER) 
    # mf.chkfile = chk_path
    
    if not lib.os.path.exists(chk_path):
        print(f"找不到文件: {chk_path}")
        return

    print(f"正在从 {chk_path} 读取轨道数据...")
    scf_data = chkfile.load(chk_path, 'scf')
    mo_coeff = scf_data['mo_coeff']
    mo_occ = scf_data['mo_occ']

    # 3. 定义活性空间 (4f 轨道)
    ncas = 7      # 7个 4f 轨道
    nelecas = 11  # Er3+ 的 11个 4f 电子
    
    # 4. 轨道排序/选择 (关键步骤)
    # 在 CASCI 之前，必须确保活性空间轨道在 mo_coeff 的中心。
    # 建议使用 AVAS 或根据轨道占据数/能量手动筛选
    # 这里演示手动切片逻辑（假设 4f 轨道在 HOMO 附近）
    ncore = np.where(mo_occ == 2)[0].shape[0] 
    # 实际上对于 Er3+，建议使用师兄代码里的 myavas 自动寻找 4f
    
    # 5. 执行 CASCI
    # 我们计算前 2 个态以获得基态和第一激发态
    n_states = 2
    mycas = mcscf.CASCI(mf, ncas, nelecas)
    mycas.fcisolver.nroots = n_states
    
    print("开始 CASCI 计算...")
    e_tot = mycas.kernel(mo_coeff)[0]
    
    # 6. 能量差分析
    if n_states > 1:
        e_ground = mycas.e_states[0]
        e_excited = mycas.e_states[1]
        delta_e_hartree = e_excited - e_ground
        
        # 单位转换
        ha2ev = 27.2114
        ha2cm = 219474.63
        
        print("\n" + "="*30)
        print(f"基态能量: {e_ground:.8f} Hartree")
        print(f"第一激发态能量: {e_excited:.8f} Hartree")
        print(f"能量差 (ΔE):")
        print(f"  {delta_e_hartree:.8f} Hartree")
        print(f"  {delta_e_hartree * ha2ev:.6f} eV")
        print(f"  {delta_e_hartree * ha2cm:.2f} cm^-1")
        print("="*30)

# 使用示例
# run_casci_from_chk('Er_CAHF_RDIIS.chk')