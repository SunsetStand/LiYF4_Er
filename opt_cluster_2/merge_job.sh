#!/bin/bash
#SBATCH -p amd
#SBATCH -o merge.log
#SBATCH -e merge.log
#SBATCH -J merge_big2
#SBATCH -n 1
#SBATCH --time=02:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh
$PYEXEC /tmp/build_merged_big.py
