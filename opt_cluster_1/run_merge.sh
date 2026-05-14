#!/bin/bash
#SBATCH -p amd
#SBATCH -o merge.out
#SBATCH -J merge_cluster
#SBATCH -e merge.err

source /data/home/wangcx/LiYF4_Er3+/miniforge3/bin/activate
conda activate /data/home/wangcx/LiYF4_Er3+/env
cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_2
python build_merged_cluster.py
