#!/bin/bash

system=$1
pull_id=$2

conf_dir="conf_files"
base_dir="${system}/system"
pull_dir="${system}/pulling/pull_${pull_id}/pull_structs"

for i in {1..26}; do

    mkdir -p us_dir/window_${i}
    
    gmx grompp -f ${conf_dir}/us_${system}.mdp -c ${pull_dir}/structure_${i}.gro -n ${base_dir}/index.ndx -p ${base_dir}/prepared_complex.top -o us_dir/window_${i}/umbrella_sampling.tpr
    gmx mdrun -deffnm us_dir/window_${i}/umbrella_sampling

done
