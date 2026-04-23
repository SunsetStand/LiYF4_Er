#!/bin/bash

# get the current directory "no_opt_x", x is the current round number
current_dir=$(basename "$PWD")
# extract the number x from "no_opt_x"
current_round=$(echo "$current_dir" | grep -oP '(?<=no_opt_)\d+')
# calculate the next round number
next_round=$((current_round + 1))
# create the next round directory "no_opt_next_round"
next_dir="no_opt_$next_round"
cd ..
mkdir "$next_dir"

# copy the files needed for the next round to the new directory
cp "$current_dir/next_round.sh" "$next_dir/"
cp "$current_dir/mv.sh" "$next_dir/"
cp "$current_dir/calculate_the_total_charge.py" "$next_dir/"
cp "$current_dir/input_file.yaml" "$next_dir/"
cp "$current_dir/job_submit.sh" "$next_dir/"
cp "$current_dir/LiYF4 (1).poscar" "$next_dir/"

# update the radii of cluster in input_file.yaml for the next round
# eg: "rCluster: 4.0" -> "rCluster: 4.5"
# extract the current rCluster value from input_file.yaml
current_rCluster=$(grep -oP '(?<=rCluster: )\d+(\.\d+)?' "$current_dir/input_file.yaml")
# calculate the next rCluster value by adding 0.5 to the current value
next_rCluster=$(echo "$current_rCluster + 0.5" | bc)
# update the rCluster value in input_file.yaml for the next round
sed -i "s/rCluster: $current_rCluster/rCluster: $next_rCluster/" "$next_dir/input_file.yaml"

# submit the job
cd "$next_dir"
sbatch job_submit.sh
