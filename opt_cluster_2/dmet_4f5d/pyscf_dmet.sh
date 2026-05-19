#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_dmet.out
#SBATCH -e pyscf_dmet.err
#SBATCH -J optbig_dmet_4f5d
#SBATCH -n 16
#SBATCH --time=72:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "===================================================================="
echo "  opt_cluster_2 — DMET Er_all, 4f+5d active space"
echo "  Host: $(hostname)  Date: $(date)"
echo "===================================================================="

$PYEXEC -c "import pyscf; from embed_sim import myavas, sacasscf_mixer, siso; from embed_sim.ssdmet import SSDMET; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "DMET 4f+5d Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_DMET.py

echo "DMET 4f+5d Finished at $(date)" >> JobProcessing.state
