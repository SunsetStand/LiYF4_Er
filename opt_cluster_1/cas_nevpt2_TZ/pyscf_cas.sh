#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_cas.out
#SBATCH -e pyscf_cas.err
#SBATCH -J ErF8_TZ_cas
#SBATCH -n 16
#SBATCH --time=48:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "===================================================================="
echo "  opt_cluster_1 — 9-atom CASSCF+NEVPT2+SISO, F=DEF2-TZVPP"
echo "  Host: $(hostname)  Date: $(date)"
echo "===================================================================="

$PYEXEC -c "import pyscf; from embed_sim import myavas, sacasscf_mixer, siso; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "CASSCF TZ Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_CASSCF_SCEI.py

echo "CASSCF TZ Finished at $(date)" >> JobProcessing.state
