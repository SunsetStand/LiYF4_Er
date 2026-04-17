import sys, os
from embed_sim import rdiis, myavas, siso

from src.AIMP3_DMET_SCEI import AIMPEnvLoader, AIMP_RHF, AIMP_RKS, AIMP_ROHF, AIMP_ROKS, AIMP_CAHF
from src.pckit2 import OrganicPCLoader, PointChargeParams
from src.EnvGenerator import XYZParser
from pyscf import gto, lib, mcscf
from pyscf.lib import chkfile
import yaml, getopt
import numpy as np
import basis_set_exchange as bse







Ha2eV = 27.211324570273

workdir = "./"
inputdir = workdir + "input.yaml"
with open(inputdir, 'r') as f:
    inputdict = yaml.safe_load(f)
aimpdict = inputdict["aimp"]
clusterdict = inputdict["cluster"]

# Load AIMP environment
try: aimpdir = workdir + aimpdict["dir"]
except KeyError: aimpdir = workdir + "aimp.xyz"
aimpdict['dir'] = aimpdir
aimpdict['workdir'] = workdir
AIMP_LOADER = AIMPEnvLoader(aimpdict)
try: orthoreg = aimpdict["orthoreg"]
except: orthoreg = 0

# Load Point Charge Environment
try: 
    pcdict = inputdict["pointcharge"]
    try: is_organic = pcdict['organic']
    except KeyError: is_organic = False

    try: rawxyzdir = workdir + pcdict["rawxyzdir"]
    except KeyError: rawxyzdir = workdir + "rawChgs.xyz"

    if not is_organic:
        try: rawchgdir = workdir + pcdict["rawchgdir"]
        except KeyError: rawchgdir = workdir + "rawCharges.dat"
        pcparam_raw = PointChargeParams(rawxyzdir, rawchgdir)
    else:
        pcdict = aimpdict
        pcdict["dir"] = rawxyzdir
        pcdict['workdir'] = workdir
        pc_loader = OrganicPCLoader(pcdict)
        pcparam_raw = pc_loader.make_param()

    try: surfxyzdir = workdir + pcdict["surfxyzdir"]
    except KeyError: surfxyzdir = workdir + "surfChgs.xyz"
    try: surfchgdir = workdir + pcdict["surfchgdir"]
    except KeyError: surfchgdir = workdir + "surfaceCharges.dat"
    pcparam_surf = PointChargeParams(surfxyzdir, surfchgdir)

    PCPARAM = pcparam_raw + pcparam_surf

except: PCPARAM = None


# Load cluster molecule
try: ecp = clusterdict["ecp"]
except KeyError: ecp = {}
try: spin = clusterdict["spin"]
except KeyError: spin = 0

def aimp_calc(CLUS_MOL, scfdict, interpolation=False):

    scftype = scfdict["calc"].upper()

    if scftype in ['HF', 'RHF', 'TDHF', 'TDAHF']:
        MF = AIMP_RHF(CLUS_MOL, AIMP_LOADER)
    elif scftype in ['ROHF']:
        MF = AIMP_ROHF(CLUS_MOL, AIMP_LOADER).density_fit()
    elif scftype in ['DFT', 'KS', 'RKS', 'TDDFT', 'TDKS']:
        try: xc = clusterdict['scf']['xc']
        except KeyError: xc = 'b3lyp'
        MF = AIMP_RKS(CLUS_MOL, AIMP_LOADER, xc=xc)
    elif scftype in ['ROKS']:
        try: xc = clusterdict['scf']['xc']
        except KeyError: xc = "b3lyp"
        MF = AIMP_ROKS(CLUS_MOL, AIMP_LOADER, xc=xc)
    elif scftype in ['CAHF']:
        try: ortho = clusterdict['scf']['orthoreg']
        except KeyError: ortho = 0
        MF = AIMP_CAHF(CLUS_MOL, AIMP_LOADER)
        MF.set_orthoreg_param(ortho)
    else:
        raise NotImplementedError("Other tags have not been implemented yet!")

    if PCPARAM is not None: MF.addPCParam2(PCPARAM)
    MF.set_orthoreg_param(orthoreg)

    # checkpoint file for TDA method
    try: MF.chkfile = workdir + inputdict['chkfile']
    except KeyError: pass
    
    return MF

# AIMP Energy calculation
if inputdict["type"].upper() in ["GEO_OPT", "GEOM_OPT", "RELAX", "GEOMOPT", "ENERGY"]:
    try: clusterdir = workdir + clusterdict["dir"]
    except KeyError: clusterdir = workdir + "cluster.xyz"

    # added for custom basis sets in the NWChem format
    for basis in clusterdict["basis"]:
        if 'parse' in clusterdict["basis"][basis] or 'load' in clusterdict["basis"][basis]:
            clusterdict["basis"][basis] = eval(clusterdict["basis"][basis])
    
    CLUS_MOL = gto.M(atom=clusterdir, basis=clusterdict["basis"], charge=clusterdict["charge"], spin=spin, verbose=4)
    scfdict = clusterdict["scf"]
    MF = aimp_calc(CLUS_MOL, scfdict)
    # MF.diis = rdiis.RDIIS(rdiis_prop='dS', imp_idx=CLUS_MOL.search_ao_label(['Ce.*']),power=0.2)
    MF.max_cycle = 3000
    MF.conv_tol = 1e-07
    MF.level_shift = 1.0
    if os.path.exists(MF.chkfile):
        print("Load from chk file.")
        MF.init_guess = 'chk'
        scfdat = chkfile.load(MF.chkfile,'scf')
        MF.e_tot = scfdat['e_tot']
        MF.mo_coeff = scfdat['mo_coeff']
        MF.mo_occ = scfdat['mo_occ']
        MF.mo_energy = scfdat['mo_energy']
        #MF.kernel()
    else:
        MF.init_guess = 'atom'
        print()
        print("Ready for MF.kernel().")
        MF.kernel()

print()

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
n_states = 5
mycas = mcscf.CASCI(MF, ncas_set, nelec_set)
mycas.fcisolver.spin = spin
mycas.fcisolver.nroots = n_states

print(f"正在求解前 {n_states} 个态...")
mycas.kernel(mo)

title = "LiYF4:Er3+"

mysiso = siso.SISO(title, mycas, amfi=True, verbose=6)
mysiso.kernel()

mysiso.analyze()

# --- 5. 输出结果与能级分析 ---
if hasattr(mysiso, 'e_states') and mysiso.e_states is not None:
    e_states = np.array(mysiso.e_states)
else:
    # PySCF 执行多态计算后，e_tot 本身就是一个数组
    e_states = np.atleast_1d(mysiso.e_tot)

print(f"成功获取到 {len(e_states)} 个态的能量")

Ha2cm = 219474.63
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


# # All-electron CASSCF+NEVPT2
# from embed_sim import myavas, sacasscf_mixer, siso

# title = 'Ce'
# ncas, nelec, mo = myavas.avas(MF, ['Ce 4f','Ce 5d'], minao=CLUS_MOL._basis['Ce'], threshold=0.5, openshell_option=2)
# print("This is ncas", ncas)
# print("This is nelec", nelec)
# ncas = 12
# nelec = 1


# mycas = sacasscf_mixer.sacasscf_mixer(MF, ncas, nelec, statelis=[0, 12, 0])
# mycas.kernel(mo)
# Ha2cm = 219474.63
# np.savetxt(title+'_cas_NO_SOC.txt',(mycas.fcisolver.e_states-np.min(mycas.fcisolver.e_states))*Ha2cm,fmt='%.6f')


# #NVEPT2

# ecorr = sacasscf_mixer.sacasscf_nevpt2(mycas, method='SC')
# mycas.fcisolver.e_states = mycas.fcisolver.e_states + ecorr
# np.savetxt(title+'_nevpt2.txt',ecorr)

# Ha2cm = 219474.63
# np.savetxt(title+'_opt.txt',(mycas.fcisolver.e_states-np.min(mycas.fcisolver.e_states))*Ha2cm,fmt='%.6f')


# #mysiso = siso.SISO(title, mycas, amfi=True, verbose=6).density_fit()
# #mysiso.kernel()