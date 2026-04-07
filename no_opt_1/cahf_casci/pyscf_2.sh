#!/bin/bash
#SBATCH -p amd
#SBATCH -o pyscf_2.out
#SBATCH -e pyscf_2.err
#SBATCH -J pyscf_2
#SBATCH -n 16

PYEXEC="/data/home/wangcx/LiYF4_Er3+/env/bin/python"
# 3. 【核心修复】加入额外的拓展包路径
# 假设你的拓展包在 /data/home/wangcx/custom_packages 目录下
export LAMP_ROOT=/data/home/wangcx/LAMP_emb
export PYTHONPATH=$PYTHONPATH:$LAMP_ROOT
export PYTHONPATH=$PYTHONPATH:/data/home/wangcx/LAMP_emb/embed_sim
export PYTHONPATH=$PYTHONPATH:/data/home/wangcx/.local/lib/python3.10/site-packages
# 如果还有其他特定的源码目录，继续累加：
# export PYTHONPATH=$PYTHONPATH:/你的/拓展包/绝对路径

# 4. 线程设置
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16

# 5. 进入目录
cd /data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/no_opt || exit 1

echo "Running on host: $(hostname)"
echo "Python Executable: $PYEXEC"
echo "PYTHONPATH: $PYTHONPATH"

# 测试 Python 是否能启动且不报错
$PYEXEC -c "import sys; print('Python starts successfully!')"
if [ $? -ne 0 ]; then
    echo "Fatal Error: Python failed to start."
    exit 1
fi

# 测试 PySCF 和 LAMP_emb 是否能导入
$PYEXEC -c "import pyscf; from embed_sim import cahf; print('Modules loaded successfully!')"
if [ $? -ne 0 ]; then
    echo "Fatal Error: Modules not found. Check PYTHONPATH."
    exit 1
fi

# 6. 任务开始
touch JobProcessing.state
echo "Job Running at $(date)" >> JobProcessing.state

/data/home/wangcx/LiYF4_Er3+/env/bin/python -u cahf_no_opt.py

# 7. 任务结束
echo "Job $SLURM_JOB_ID $SLURM_JOB_NAME Finished !" >> $HOME/finish
echo "Job Finished at $(date)" >> JobProcessing.state