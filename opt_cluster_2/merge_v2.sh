#!/bin/bash
#SBATCH -p amd
#SBATCH -o merge_v2.log
#SBATCH -e merge_v2.log
#SBATCH -J merge_big3
#SBATCH -n 4
#SBATCH --time=04:00:00

source /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/setup_env.sh
$PYEXEC /tmp/build_merged_big_v2.py
