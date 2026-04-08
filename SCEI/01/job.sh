#!/bin/bash


#F
cd 01F/
sbatch job_submit.sh
cd ../

#Li
cd 02Li/
sbatch job_submit.sh
cd ../

#Y
cd 03Y/
sbatch job_submit.sh
cd ../