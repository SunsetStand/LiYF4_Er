#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_dmet.out
#SBATCH -e pyscf_dmet.err
#SBATCH -J opt2_dmet
#SBATCH -n 16
#SBATCH --time=48:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "Running on host: $(hostname)"
echo "opt_cluster_2 DMET Er_all (25-atom big cluster)"

$PYEXEC -c "import pyscf; from embed_sim.ssdmet import SSDMET; from src.AIMP3_DMET_SCEI import AIMP_ROHF; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "opt2 DMET Er_all Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_DMET.py -i ./

echo "opt2 DMET Er_all Finished at $(date)" >> $HOME/finish
echo "opt2 DMET Er_all Finished at $(date)" >> JobProcessing.state
