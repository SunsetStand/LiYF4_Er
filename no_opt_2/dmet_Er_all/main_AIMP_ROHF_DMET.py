import sys, os
from embed_sim import myavas, sacasscf_mixer, siso
from embed_sim.ssdmet import SSDMET

from src.AIMP3_DMET_SCEI import AIMPEnvLoader, AIMP_RHF, AIMP_RKS, AIMP_ROHF, AIMP_ROKS, AIMP_CAHF
from src.pckit2 import OrganicPCLoader, PointChargeParams
from pyscf import gto, lib
from pyscf.lib import chkfile
import yaml
import numpy as np
import basis_set_exchange as bse

Ha2eV = 27.211324570273
Ha2cm = 219474.63

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
    try: MF.chkfile = workdir + inputdict['chkfile']
    except KeyError: pass
    return MF

# ===================================================================
# ROHF — load from chk file
# ===================================================================
if inputdict["type"].upper() in ["GEO_OPT", "GEOM_OPT", "RELAX", "GEOMOPT", "ENERGY"]:
    try: clusterdir = workdir + clusterdict["dir"]
    except KeyError: clusterdir = workdir + "cluster.xyz"
    for basis in clusterdict["basis"]:
        if 'parse' in clusterdict["basis"][basis] or 'load' in clusterdict["basis"][basis]:
            clusterdict["basis"][basis] = eval(clusterdict["basis"][basis])
    CLUS_MOL = gto.M(atom=clusterdir, basis=clusterdict["basis"],
                     charge=clusterdict["charge"], spin=spin, verbose=4)
    scfdict = clusterdict["scf"]
    MF = aimp_calc(CLUS_MOL, scfdict)

    # Load ROHF from chk
    if not os.path.exists(MF.chkfile):
        raise FileNotFoundError(f"ROHF chk file not found: {MF.chkfile}")
    print(f"Loading ROHF from chk: {MF.chkfile}")
    scfdat = chkfile.load(MF.chkfile, 'scf')
    MF.e_tot = scfdat['e_tot']
    MF.mo_coeff = scfdat['mo_coeff']
    MF.mo_occ = scfdat['mo_occ']
    MF.mo_energy = scfdat['mo_energy']
    print(f"Loaded: e_tot={MF.e_tot:.8f}")

print()

# ===================================================================
# DMET — impurity: Er (all orbitals)
# ===================================================================
IMP_NAME = 'Er_all'
IMP_LABELS = ['Er']
IMP_DESC = 'All Er orbitals'

title_base = 'LiYF4:Er3+'
ncas_set = 7
nelec_set = 11

title = f'{title_base}_DMET_{IMP_NAME}'
print("=" * 70)
print(f"  DMET impurity: {IMP_DESC}")
print(f"  Title: {title}")
print("=" * 70)

# Step 1: Build DMET embedding
print(f"\n[Step 1] Building DMET embedding")
dmet = SSDMET(MF, title=title, imp_idx=IMP_LABELS, threshold=1e-13)
es_mf = dmet.build()
print(f"  Embedded space: {dmet.nes} orbitals")
print(f"  impurity: {len(dmet.imp_idx)} orbitals")
print(f"  frozen occupied: {dmet.nfo},  frozen virtual: {dmet.nfv}")
print(f"  Frozen orbital energy: {dmet.fo_ene():.8f} Hartree")

# Step 2: AVAS
print(f"\n[Step 2] AVAS: Er 4f from embedded space")
ncas, nelec, es_mo = dmet.avas(['Er 4f'],
                                minao=CLUS_MOL._basis['Er'],
                                threshold=0.5,
                                openshell_option=2)
print(f"  AVAS: ncas={ncas}, nelec={nelec}")

# Step 3: SACASSCF
print(f"\n[Step 3] SACASSCF({ncas_set},{nelec_set})")
statelis = [0, 0, 0, 35]
es_mycas = sacasscf_mixer.sacasscf_mixer(es_mf, ncas_set, nelec_set, statelis=statelis)
es_mycas.kernel(es_mo)
e_cas = es_mycas.fcisolver.e_states
e_min = np.min(e_cas)
print(f"  CASSCF E_min = {e_min:.8f} H")
np.savetxt(f'{title}_cas_NO_SOC.txt', (e_cas - e_min) * Ha2cm, fmt='%.6f')

# Step 4: NEVPT2
print(f"\n[Step 4] NEVPT2")
ecorr = sacasscf_mixer.sacasscf_nevpt2(es_mycas, method='SC')
es_mycas.fcisolver.e_states = e_cas + ecorr
np.savetxt(f'{title}_nevpt2.txt', ecorr)
np.savetxt(f'{title}_opt.txt',
           (es_mycas.fcisolver.e_states - np.min(es_mycas.fcisolver.e_states)) * Ha2cm,
           fmt='%.6f')
print(f"  NEVPT2 corrections: {np.min(ecorr)*Ha2cm:.2f} to {np.max(ecorr)*Ha2cm:.2f} cm^-1")

# Step 5: SISO
print(f"\n[Step 5] SISO")
total_cas = dmet.total_cas(es_mycas)
mysiso = siso.SISO(title, total_cas, amfi=True, verbose=6)
mysiso.kernel()
mysiso.analyze()
print(f"  SISO done → {title}_mag.txt")

print(f"\n  [{IMP_NAME}] finished!\n")
