import sys, os
import yaml
import numpy as np
from pyscf import gto, lib, mcscf
from pyscf.lib import chkfile
from embed_sim import rdiis, myavas
import basis_set_exchange as bse

# 导入 AIMP 相关的自定义类
from src.AIMP3_DMET_SCEI import AIMPEnvLoader,AIMP_CAHF
from src.pckit2 import OrganicPCLoader, PointChargeParams

# --- 1. 环境初始化与参数加载 ---
Ha2cm = 219474.63
workdir = "./"
inputdir = workdir + "input.yaml"

with open(inputdir, 'r') as f:
    inputdict = yaml.safe_load(f)

aimpdict = inputdict["aimp"]
clusterdict = inputdict["cluster"]
title = 'Er' # 修改为 Er

# 加载 AIMP 环境
aimpdir = workdir + aimpdict.get("dir", "aimp.xyz")
aimpdict['dir'] = aimpdir
aimpdict['workdir'] = workdir
AIMP_LOADER = AIMPEnvLoader(aimpdict)
orthoreg = aimpdict.get("orthoreg", 0)

# 加载点电荷环境 (PCPARAM) - 保持原有逻辑
try:
    pcdict = inputdict["pointcharge"]
    rawxyzdir = workdir + pcdict.get("rawxyzdir", "rawChgs.xyz")
    if not pcdict.get('organic', False):
        rawchgdir = workdir + pcdict.get("rawchgdir", "rawCharges.dat")
        pcparam_raw = PointChargeParams(rawxyzdir, rawchgdir)
    else:
        pc_loader = OrganicPCLoader({**aimpdict, "dir": rawxyzdir})
        pcparam_raw = pc_loader.make_param()
    
    surfxyzdir = workdir + pcdict.get("surfxyzdir", "surfChgs.xyz")
    surfchgdir = workdir + pcdict.get("surfchgdir", "surfaceCharges.dat")
    PCPARAM = pcparam_raw + PointChargeParams(surfxyzdir, surfchgdir)
except:
    PCPARAM = None

# --- 2. 定义计算核心函数 ---
def run_aimp_scf(CLUS_MOL, scfdict):
    # 使用 CAHF 处理 Er3+ (S=3/2, spin=3)
    MF = AIMP_CAHF(CLUS_MOL, AIMP_LOADER)
    
    if PCPARAM is not None: 
        MF.addPCParam2(PCPARAM)
    MF.set_orthoreg_param(orthoreg)
    
    # 关键修改：设置针对 Er 的收敛参数
    MF.chkfile = workdir + inputdict.get('chkfile', 'Er_CAHF.chk')
    MF.max_cycle = 500          # 给予足够的迭代空间
    MF.level_shift = 1.0        # 核心：使用 1.0 Hartree 的位移稳定 4f 轨道
    MF.conv_tol = 1e-7
    
    # 设置 RDIIS 锁定 Er 的 4f 轨道
    f_idx = CLUS_MOL.search_ao_label(['Er.*4f'])
    MF.diis = rdiis.RDIIS(MF, rdiis_prop='dS', imp_idx=f_idx, power=0.5)
    
    return MF

# --- 3. 执行 CAHF 计算 ---
clusterdir = workdir + clusterdict.get("dir", "cluster.xyz")
spin = clusterdict.get("spin", 3) # Er3+ 基态通常 spin=3

# 载入基组
for basis_key in clusterdict["basis"]:
    if isinstance(clusterdict["basis"][basis_key], str) and ('parse' in clusterdict["basis"][basis_key]):
        clusterdict["basis"][basis_key] = eval(clusterdict["basis"][basis_key])

CLUS_MOL = gto.M(atom=clusterdir, basis=clusterdict["basis"], 
                 charge=clusterdict["charge"], spin=spin, verbose=4)

MF = run_aimp_scf(CLUS_MOL, clusterdict["scf"])

if os.path.exists(MF.chkfile):
    print(f"\n--- 从检查点文件 {MF.chkfile} 读取轨道 ---")
    MF.init_guess = 'chk'
    # 仅加载不执行，直接进入 CAS 阶段
else:
    print("\n--- 开始 CAHF 计算 (Level Shift = 1.0) ---")
    MF.init_guess = 'atom'

MF.kernel()

# --- 4. 执行 CASCI 计算 ---
print("\n--- 开始 CASCI 计算分析 ---")

# 使用 AVAS 自动提取 Er 的 7 个 4f 轨道
# Er3+: 11个 4f 电子 -> nelec = 11, ncas = 7
ncas_set = 7
nelec_set = 11

ncas, nelec, mo = myavas.avas(MF, ['Er 4f'], 
                              minao=CLUS_MOL._basis['Er'], 
                              threshold=0.5, 
                              openshell_option=2)

print(f"AVAS 提取结果: ncas={ncas}, nelec={nelec}")

# 定义 CASCI (计算基态和第一激发态)
n_states = 2
mycas = mcscf.CASCI(MF, ncas_set, nelec_set)
mycas.fcisolver.nroots = n_states

print(f"正在求解前 {n_states} 个态...")
mycas.kernel(mo)

# --- 5. 输出结果与能级分析 ---
e_states = mycas.e_states
delta_e_cm = (e_states[1] - e_states[0]) * Ha2cm

print("\n" + "="*40)
print(f"{title} 基态能量: {e_states[0]:.8f} Hartree")
print(f"{title} 第一激发态能量: {e_states[1]:.8f} Hartree")
print(f"基态-激发态能级差: {delta_e_cm:.2f} cm^-1")
print("="*40)

# 保存能级数据
np.savetxt(f'{title}_casci_levels.txt', (e_states - np.min(e_states)) * Ha2cm, 
           fmt='%.6f', header='Energy levels in cm^-1')
print(f"能级数据已保存至 {title}_casci_levels.txt")