#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_dmet.out
#SBATCH -e pyscf_dmet.err
#SBATCH -J opt_dmet_ErF8
#SBATCH -n 16
#SBATCH --time=48:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "Running on host: $(hostname)"
echo "opt_cluster_1 DMET ErF8"

$PYEXEC -c "import pyscf; from embed_sim.ssdmet import SSDMET; from src.AIMP3_DMET_SCEI import AIMP_ROHF; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "opt DMET ErF8 Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_DMET.py -i ./

echo "opt DMET ErF8 Finished at $(date)" >> $HOME/finish
echo "opt DMET ErF8 Finished at $(date)" >> JobProcessing.state
