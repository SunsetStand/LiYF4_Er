#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_cas.out
#SBATCH -e pyscf_cas.err
#SBATCH -J optbig_casci
#SBATCH -n 32
#SBATCH --time=72:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

# Use 32 OpenMP threads for PySCF parallelism
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export OPENBLAS_NUM_THREADS=32

echo "===================================================================="
echo "  opt_cluster_2 — 25-atom CASSCF+NEVPT2+SISO (NO DMET control)"
echo "  Host: $(hostname)"
echo "  Date: $(date)"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  Convergence tuning: max_cycle_macro=150, level_shift=0.3, ah_level_shift=1e-3"
echo "===================================================================="

$PYEXEC -c "import pyscf; from embed_sim import myavas, sacasscf_mixer, siso; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "CASSCF+NEVPT2+SISO Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_CASSCF_SCEI.py

echo "CASSCF+NEVPT2+SISO Finished at $(date)" >> JobProcessing.state
