#!/bin/bash

mkdir cahf_casci

cp ../no_opt_1/cahf_casci/input.yaml cahf_casci/
cp ../no_opt_1/cahf_casci/pyscf_2.sh cahf_casci/
cp ../no_opt_1/cahf_casci/main_AIMP_CAHF_CASCI_SCEI.py cahf_casci/
cp ../no_opt_1/cahf_casci/MO_* cahf_casci/

mv cluster.xyz aimp.xyz rawChgs.xyz surfChgs.xyz rawCharges.dat surfaceCharges.dat cahf_casci/

total_charge=$(python calculate_the_total_charge.py)

#将cluster.xyz中的第一个Y原子替换为Er原子
sed -i '0,/Y/{s/Y/Er/}' cahf_casci/cluster.xyz

sed -i "s/charge: -5/charge: $total_charge/g" cahf_casci/input.yaml

cd cahf_casci
sbatch pyscf_2.sh