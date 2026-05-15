#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_dmet.out
#SBATCH -e pyscf_dmet.err
#SBATCH -J dmet_Er4f
#SBATCH -n 16
#SBATCH --time=24:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "Running on host: $(hostname)"
echo "DMET impurity: Er 4f"

$PYEXEC -c "import pyscf; from embed_sim.ssdmet import SSDMET; from src.AIMP3_DMET_SCEI import AIMP_ROHF; print('Modules OK')" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "DMET Er_4f Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_DMET.py -i ./

echo "DMET Er_4f Finished at $(date)" >> $HOME/finish
echo "DMET Er_4f Finished at $(date)" >> JobProcessing.state
