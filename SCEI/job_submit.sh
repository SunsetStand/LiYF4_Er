#!/bin/bash
#SBATCH -p amd
#SBATCH -o model.out
#SBATCH -J build_model
#SBATCH -e model.err
#SBATCH -n 8

#Doing job
touch JobProcessing.state
echo "Job Running at time" >> JobProcessing.state
echo `date` >> JobProcessing.state

#Task
source /data/home/wangcx/LiYF4_Er3+/miniforge3/bin/activate
conda activate /data/home/wangcx/LiYF4_Er3+/env

cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/SCEI/01F
python /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/main.py -i input_file.yaml

cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/SCEI/02Li
python /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/main.py -i input_file.yaml

cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/SCEI/03Y
python /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/main.py -i input_file.yaml

#when job finished
echo "Job $SLURM_JOB_ID $SLURM_JOB_NAME Finished !" >> $HOME/finish
echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo "Job Finished at time" >> JobProcessing.state
echo `date` >> JobProcessing.state  