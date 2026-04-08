#!/bin/bash
#SBATCH -p amd
#SBATCH -o scei.out
#SBATCH -J scei-F
#SBATCH -e scei.err
#SBATCH -n 16

#Doing job
touch JobProcessing.state
echo "Job Running at time" >> JobProcessing.state
echo `date` >> JobProcessing.state

#Task
source /data/home/wangcx/LiYF4_Er3+/miniforge3/bin/activate
conda activate /data/home/wangcx/LiYF4_Er3+/env

python /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/SCEI_r_w.py -i ./

#when job finished
echo "Job $SLURM_JOB_ID $SLURM_JOB_NAME Finished !" >> $HOME/finish
echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo "Job Finished at time" >> JobProcessing.state
echo `date` >> JobProcessing.state  