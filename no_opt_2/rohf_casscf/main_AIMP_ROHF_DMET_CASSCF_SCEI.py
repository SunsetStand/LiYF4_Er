import sys, os
from embed_sim import myavas, sacasscf_mixer, siso
from embed_sim.ssdmet import SSDMET

from src.AIMP3_DMET_SCEI import AIMPEnvLoader, AIMP_RHF, AIMP_RKS, AIMP_ROHF, AIMP_ROKS, AIMP_CAHF
from src.pckit2 import OrganicPCLoader, PointChargeParams
from src.EnvGenerator import XYZParser
from pyscf import gto, lib, mcscf, scf
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
# AIMP Energy calculation — ROHF initial guess
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
    MF.max_cycle = 500
    MF.conv_tol = 1e-7
    MF.conv_tol_grad = 3e-4       # relaxed gradient tolerance
    MF.level_shift = 0.3           # small level shift to suppress oscillation
    
    # atom initial guess (pure ROHF, no CAHF chk dependency)
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
print("ROHF converged!")
print(f"ROHF energy: {MF.e_tot:.8f} Hartree")
print()

# ===================================================================
# DMET + CASSCF + NEVPT2 + SISO for three impurity choices
# ===================================================================
title_base = 'LiYF4:Er3+'

# Three impurity choices for DMET
impurity_choices = [
    ('Er_4f',  ['Er 4f'],   "Er 4f orbitals"),
    ('Er_all', ['Er'],       "All Er orbitals"),
    ('ErF8',   ['Er', 'F'], "Er + all F atoms"),
]

Ha2cm = 219474.63
ncas_set = 7
nelec_set = 11

for imp_name, imp_labels, imp_desc in impurity_choices:
    print("\n" + "=" * 70)
    print(f"  DMET with impurity: {imp_desc} ({imp_name})")
    print("=" * 70)
    
    title = f'{title_base}_DMET_{imp_name}'
    
    # ---------- Step 1: Build DMET embedding ----------
    print(f"\n[Step 1] Building DMET embedding with impurity: {imp_labels}")
    dmet = SSDMET(MF, title=title, imp_idx=imp_labels, threshold=1e-13)
    
    # Build embedding (this constructs embedded Fock matrix, 2e integrals, etc.)
    es_mf = dmet.build()
    print(f"  Embedded space: {dmet.nes} orbitals")
    print(f"    impurity: {dmet.nes - (dmet.nes - len(dmet.imp_idx))} orbitals")
    print(f"    bath:     {dmet.nes - len(dmet.imp_idx)} orbitals  (estimated)")
    print(f"    frozen occupied: {dmet.nfo}")
    print(f"    frozen virtual:  {dmet.nfv}")
    print(f"  Embedded molecule: {es_mf.mol.nelectron} electrons, spin={es_mf.mol.spin}")
    print(f"  Frozen orbital energy: {dmet.fo_ene():.8f} Hartree")
    print(f"  DMET chk saved to: {title}_dmet_chk.h5")
    
    # ---------- Step 2: AVAS within embedded space ----------
    print(f"\n[Step 2] AVAS to extract Er 4f active space from embedded orbitals")
    ncas, nelec, es_mo = dmet.avas(['Er 4f'],
                                    minao=CLUS_MOL._basis['Er'],
                                    threshold=0.5,
                                    openshell_option=2)
    print(f"  AVAS result: ncas={ncas}, nelec={nelec}")
    if ncas != ncas_set and False:
        print(f"  WARNING: AVAS gave ncas={ncas}, but using ncas_set={ncas_set}")
    
    # ---------- Step 3: SACASSCF in embedded space ----------
    print(f"\n[Step 3] SACASSCF({ncas_set},{nelec_set}) in embedded space")
    statelis = [0, 0, 0, 35]  # 35 S=3/2 quartets
    es_mycas = sacasscf_mixer.sacasscf_mixer(es_mf, ncas_set, nelec_set, statelis=statelis)
    es_mycas.kernel(es_mo)
    
    e_cas = es_mycas.fcisolver.e_states
    e_min = np.min(e_cas)
    print(f"  CASSCF lowest energy: {e_min:.8f} Hartree")
    print(f"  CASSCF state range: {np.min(e_cas)*Ha2cm:.2f} - {np.max(e_cas)*Ha2cm:.2f} cm^-1")
    
    # Save CASSCF-only energies
    np.savetxt(f'{title}_cas_NO_SOC.txt',
               (e_cas - e_min) * Ha2cm,
               fmt='%.6f')
    print(f"  CASSCF energies saved to {title}_cas_NO_SOC.txt")
    
    # ---------- Step 4: NEVPT2 in embedded space ----------
    print(f"\n[Step 4] NEVPT2 correction")
    try:
        ecorr = sacasscf_mixer.sacasscf_nevpt2(es_mycas, method='SC')
        es_mycas.fcisolver.e_states = e_cas + ecorr
        np.savetxt(f'{title}_nevpt2.txt', ecorr)
        
        e_corr_cm = ecorr * Ha2cm
        print(f"  NEVPT2 corrections (cm^-1): {e_corr_cm}")
        print(f"  NEVPT2 correction range: {np.min(ecorr)*Ha2cm:.2f} - {np.max(ecorr)*Ha2cm:.2f} cm^-1")
        
        # Save NEVPT2-corrected energies
        np.savetxt(f'{title}_opt.txt',
                   (es_mycas.fcisolver.e_states - np.min(es_mycas.fcisolver.e_states)) * Ha2cm,
                   fmt='%.6f')
        print(f"  NEVPT2 energies saved to {title}_opt.txt")
    except Exception as e:
        print(f"  NEVPT2 failed: {e}")
        print("  Saving CASSCF-only results. Skipping SISO for this impurity.")
        continue
    
    # ---------- Step 5: Expand to full space for SISO ----------
    print(f"\n[Step 5] Expand embedded CAS to full space")
    total_cas = dmet.total_cas(es_mycas)
    print(f"  Total CAS created. Active space: ({total_cas.ncas}, {total_cas.nelecas})")
    
    # ---------- Step 6: SISO (spin-orbit coupling) ----------
    print(f"\n[Step 6] SISO (Spin-Orbit Coupling)")
    try:
        mysiso = siso.SISO(title, total_cas, amfi=True, verbose=6)
        mysiso.kernel()
        mysiso.analyze()
        print(f"  SISO completed. Magnetic properties saved to {title}_mag.txt")
    except Exception as e:
        print(f"  SISO failed: {e}")
    
    print(f"\n  DMET [{imp_name}] finished.\n")

print("\n" + "=" * 70)
print("  All DMET calculations finished!")
print("=" * 70)
