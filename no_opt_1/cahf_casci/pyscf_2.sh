#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_3.out
#SBATCH -e pyscf_3.err
#SBATCH -J pyscf_3
#SBATCH -n 16

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "Running on host: $(hostname)"
echo "Python Executable: $PYEXEC"
echo "PYTHONPATH: $PYTHONPATH"

$PYEXEC -c "import sys; print(\"Python OK\")" || { echo "FATAL"; exit 1; }
$PYEXEC -c "import pyscf; from embed_sim import cahf; from src.AIMP3_DMET_SCEI import AIMP_CAHF; print(\"Modules OK\")" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "Job Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_CAHF_CASCI_SCEI.py -i ./

echo "Job $SLURM_JOB_ID $SLURM_JOB_NAME Finished !" >> $HOME/finish
echo "Job Finished at $(date)" >> JobProcessing.state
