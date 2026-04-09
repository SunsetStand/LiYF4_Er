#!/bin/bash

mkdir cahf_casci

cp ../no_opt_1/cahf_casci/input.yaml cahf_casci/
cp ../no_opt_1/cahf_casci/pyscf_2.sh cahf_casci/
cp ../no_opt_1/cahf_casci/main_AIMP_CAHF_CASCI_SCEI.py cahf_casci/
cp ../no_opt_1/cahf_casci/MO_* cahf_casci/

mv cluster.xyz aimp.xyz rawChgs.xyz surfChgs.xyz rawCharges.dat surfaceCharges.dat cahf_casci/

cd cahf_casci
sbatch pyscf_2.sh