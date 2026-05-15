#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_dmet.out
#SBATCH -e pyscf_dmet.err
#SBATCH -J pyscf_dmet_rohf
#SBATCH -n 16
#SBATCH --time=48:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

echo "Running on host: $(hostname)"
echo "Python Executable: $PYEXEC"
echo "PYTHONPATH: $PYTHONPATH"

$PYEXEC -c "import sys; print(\"Python OK\")" || { echo "FATAL"; exit 1; }
$PYEXEC -c "import pyscf; from embed_sim import myavas, sacasscf_mixer, siso; from embed_sim.ssdmet import SSDMET; from src.AIMP3_DMET_SCEI import AIMP_ROHF; print(\"Modules OK\")" || { echo "FATAL"; exit 1; }

touch JobProcessing.state
echo "DMET Job Running at $(date)" >> JobProcessing.state

$PYEXEC -u main_AIMP_ROHF_DMET_CASSCF_SCEI.py -i ./

echo "DMET Job $SLURM_JOB_ID $SLURM_JOB_NAME Finished !" >> $HOME/finish
echo "DMET Job Finished at $(date)" >> JobProcessing.state
