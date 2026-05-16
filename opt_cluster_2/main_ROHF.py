import sys, os
from src.AIMP3_DMET_SCEI import AIMPEnvLoader, AIMP_RHF, AIMP_RKS, AIMP_ROHF, AIMP_ROKS, AIMP_CAHF
from src.pckit2 import OrganicPCLoader, PointChargeParams
from pyscf import gto, lib
from pyscf.lib import chkfile
import yaml
import numpy as np
import basis_set_exchange as bse

workdir = "./"
inputdir = workdir + "input.yaml"
with open(inputdir, "r") as f:
    inputdict = yaml.safe_load(f)
aimpdict = inputdict["aimp"]
clusterdict = inputdict["cluster"]

try: aimpdir = workdir + aimpdict["dir"]
except KeyError: aimpdir = workdir + "aimp.xyz"
aimpdict["dir"] = aimpdir
aimpdict["workdir"] = workdir
AIMP_LOADER = AIMPEnvLoader(aimpdict)
try: orthoreg = aimpdict["orthoreg"]
except: orthoreg = 0

try:
    pcdict = inputdict["pointcharge"]
    try: is_organic = pcdict["organic"]
    except KeyError: is_organic = False
    try: rawxyzdir = workdir + pcdict["rawxyzdir"]
    except KeyError: rawxyzdir = workdir + "rawChgs.xyz"
    if not is_organic:
        try: rawchgdir = workdir + pcdict["rawchgdir"]
        except KeyError: rawchgdir = workdir + "rawCharges.dat"
        pcparam_raw = PointChargeParams(rawxyzdir, rawchgdir)
    else:
        pcdict2 = dict(aimpdict)
        pcdict2["dir"] = rawxyzdir
        pcdict2["workdir"] = workdir
        pc_loader = OrganicPCLoader(pcdict2)
        pcparam_raw = pc_loader.make_param()
    try: surfxyzdir = workdir + pcdict["surfxyzdir"]
    except KeyError: surfxyzdir = workdir + "surfChgs.xyz"
    try: surfchgdir = workdir + pcdict["surfchgdir"]
    except KeyError: surfchgdir = workdir + "surfaceCharges.dat"
    pcparam_surf = PointChargeParams(surfxyzdir, surfchgdir)
    PCPARAM = pcparam_raw + pcparam_surf
except: PCPARAM = None

try: ecp = clusterdict["ecp"]
except KeyError: ecp = {}
try: spin = clusterdict["spin"]
except KeyError: spin = 0

def aimp_calc(CLUS_MOL, scfdict):
    scftype = scfdict["calc"].upper()
    if scftype in ["ROHF"]:
        MF = AIMP_ROHF(CLUS_MOL, AIMP_LOADER).density_fit()
    else:
        raise NotImplementedError
    if PCPARAM is not None: MF.addPCParam2(PCPARAM)
    MF.set_orthoreg_param(orthoreg)
    try: MF.chkfile = workdir + inputdict["chkfile"]
    except KeyError: pass
    return MF

if inputdict["type"].upper() in ["GEO_OPT", "GEOM_OPT", "RELAX", "GEOMOPT", "ENERGY"]:
    try: clusterdir = workdir + clusterdict["dir"]
    except KeyError: clusterdir = workdir + "cluster.xyz"
    for basis in clusterdict["basis"]:
        val = clusterdict["basis"][basis]
        if isinstance(val, str) and ("parse" in val or "load" in val):
            clusterdict["basis"][basis] = eval(val)
    CLUS_MOL = gto.M(atom=clusterdir, basis=clusterdict["basis"],
                     charge=clusterdict["charge"], spin=spin, verbose=4)
    scfdict = clusterdict["scf"]
    MF = aimp_calc(CLUS_MOL, scfdict)
    MF.max_cycle = 500
    MF.conv_tol = 1e-7
    MF.conv_tol_grad = 3e-4
    MF.level_shift = 0.3
    MF.init_guess = "atom"
    print("Starting ROHF...")
    MF.kernel()
    print(f"ROHF converged: e_tot = {MF.e_tot:.8f} Hartree")
