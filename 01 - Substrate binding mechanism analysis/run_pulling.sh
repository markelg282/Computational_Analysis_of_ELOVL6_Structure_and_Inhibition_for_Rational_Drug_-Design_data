#!/bin/bash

system=$1

base_dir="${system}/system"
conf_dir="conf_files"
equil_dir="${system}/equilibration"
pull_dir="${system}/pulling"

# mkdir -p $equil_dir $pull_dir

# gmx grompp -f ${conf_dir}/minimization.mdp -c ${base_dir}/prepared_complex.gro -r ${base_dir}/prepared_complex.gro -n ${base_dir}/index.ndx -p ${base_dir}/prepared_complex.top -o ${equil_dir}/min.tpr
# gmx mdrun -deffnm ${equil_dir}/min
# gmx grompp -f ${conf_dir}/nvt_equilibration.mdp -c ${equil_dir}/min.gro -r ${equil_dir}/min.gro -n ${base_dir}/index.ndx -p ${base_dir}/prepared_complex.top -o ${equil_dir}/nvt_equil.tpr
# gmx mdrun -deffnm ${equil_dir}/nvt_equil
# gmx grompp -f ${conf_dir}/npt_equilibration.mdp -c ${equil_dir}/nvt_equil.gro -r ${equil_dir}/nvt_equil.gro -n ${base_dir}/index.ndx -p ${base_dir}/prepared_complex.top -o ${equil_dir}/npt_equil.tpr
# gmx mdrun -deffnm ${equil_dir}/npt_equil

for i in {1..10}; do
    mkdir -p ${pull_dir}/pull_${i}
    gmx grompp -f ${conf_dir}/pulling_${system}.mdp -c ${equil_dir}/npt_equil.gro -n ${base_dir}/index.ndx -p ${base_dir}/prepared_complex.top -o ${pull_dir}/pull_${i}/pulling.tpr
    gmx mdrun -deffnm ${pull_dir}/pull_${i}/pulling
    python utils/extract_snapshots.py -s ${equil_dir}/npt_equil.gro -t ${pull_dir}/pull_${i}/pulling.xtc -n 26 -o ${pull_dir}/pull_${i}/pull_structs
done
