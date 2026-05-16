#!/bin/bash
#SBATCH -p amd
#SBATCH -o merge.log
#SBATCH -e merge.log
#SBATCH -J merge_cluster
#SBATCH -n 1
#SBATCH --time=01:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_1
mkdir -p ../opt_cluster_2
$PYEXEC build_merged_cluster.py

echo "=== OUTPUT ==="
ls -la ../opt_cluster_2/
