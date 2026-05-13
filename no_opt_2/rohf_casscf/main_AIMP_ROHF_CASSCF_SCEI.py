import sys, os
from embed_sim import myavas, sacasscf_mixer, siso

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

    # checkpoint file
    try: MF.chkfile = workdir + inputdict['chkfile']
    except KeyError: pass
    
    return MF

# ===================================================================
# AIMP Energy calculation — ROHF initial guess, then CASSCF
# ===================================================================
if inputdict["type"].upper() in ["GEO_OPT", "GEOM_OPT", "RELAX", "GEOMOPT", "ENERGY"]:
    try: clusterdir = workdir + clusterdict["dir"]
    except KeyError: clusterdir = workdir + "cluster.xyz"

    # Parse custom basis sets in NWChem format
    for basis in clusterdict["basis"]:
        if 'parse' in clusterdict["basis"][basis] or 'load' in clusterdict["basis"][basis]:
            clusterdict["basis"][basis] = eval(clusterdict["basis"][basis])
    
    CLUS_MOL = gto.M(atom=clusterdir, basis=clusterdict["basis"], 
                     charge=clusterdict["charge"], spin=spin, verbose=4)
    scfdict = clusterdict["scf"]
    MF = aimp_calc(CLUS_MOL, scfdict)
    
    # ---- ROHF SCF: simple setup, no RDIIS ----
    MF.max_cycle = 3000
    MF.conv_tol = 1e-07
    
    # 使用 atom 初猜（纯 ROHF，不依赖 CAHF chk）
    if os.path.exists(MF.chkfile):
        print("Load from chk file.")
        MF.init_guess = 'chk'
        scfdat = chkfile.load(MF.chkfile, 'scf')
        MF.e_tot = scfdat['e_tot']
        MF.mo_coeff = scfdat['mo_coeff']
        MF.mo_occ = scfdat['mo_occ']
        MF.mo_energy = scfdat['mo_energy']
    else:
        MF.init_guess = 'atom'
        print()
        print("Ready for ROHF MF.kernel().")
        MF.kernel()

print()

# ===================================================================
# All-electron CASSCF + NEVPT2 + SISO
# ===================================================================
from embed_sim import myavas, sacasscf_mixer, siso

title = 'LiYF4:Er3+'
ncas, nelec, mo = myavas.avas(MF, ['Er 4f'], 
                              minao=CLUS_MOL._basis['Er'], 
                              threshold=0.5, 
                              openshell_option=2)
print(f"AVAS 提取结果: ncas={ncas}, nelec={nelec}")

ncas_set = 7
nelec_set = 11

# CASSCF: 35 roots for all S=3/2 states in the 4f11 CAS(11,7) manifold
statelis = [0, 0, 0, 35]
mycas = sacasscf_mixer.sacasscf_mixer(MF, ncas_set, nelec_set, statelis=statelis)
mycas.kernel(mo)
Ha2cm = 219474.63
np.savetxt(title+'_cas_NO_SOC.txt',
           (mycas.fcisolver.e_states - np.min(mycas.fcisolver.e_states)) * Ha2cm,
           fmt='%.6f')
print(f"CASSCF 完成，能级已保存至 {title}_cas_NO_SOC.txt")

# NEVPT2
ecorr = sacasscf_mixer.sacasscf_nevpt2(mycas, method='SC')
mycas.fcisolver.e_states = mycas.fcisolver.e_states + ecorr
np.savetxt(title+'_nevpt2.txt', ecorr)
Ha2cm = 219474.63
np.savetxt(title+'_opt.txt',
           (mycas.fcisolver.e_states - np.min(mycas.fcisolver.e_states)) * Ha2cm,
           fmt='%.6f')
print(f"NEVPT2 完成，结果已保存")

# SISO (spin-orbit coupling)
mysiso = siso.SISO(title, mycas, amfi=True, verbose=6)
mysiso.kernel()
mysiso.analyze()
print(f"SISO 完成")
