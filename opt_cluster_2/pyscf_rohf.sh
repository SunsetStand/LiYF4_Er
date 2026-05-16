#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_rohf.out
#SBATCH -e pyscf_rohf.err
#SBATCH -J optbig_rohf
#SBATCH -n 16
#SBATCH --time=24:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "opt_cluster_2 ROHF"

$PYEXEC -c 'import pyscf; from src.AIMP3_DMET_SCEI import AIMP_ROHF; print("Modules OK")' || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "ROHF Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_ROHF.py -i ./

echo "ROHF Finished at $(date)" >> JobProcessing.state
