#!/bin/bash
#SBATCH -p amd
#SBATCH -o build.log
#SBATCH -e build.log
#SBATCH -J build_opt_cluster
#SBATCH -n 1
#SBATCH --time=02:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh

# Step 1: Extract big cluster from CONTCAR
$PYEXEC /tmp/build_big_cluster2.py

# Step 2: Run the surface charge fitting (from build_merged_cluster.py)
cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_1
$PYEXEC build_merged_cluster.py

echo "=== DONE ==="
ls -la /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_2/
