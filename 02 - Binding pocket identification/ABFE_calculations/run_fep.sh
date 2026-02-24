#!/bin/bash

family=$1
model=$2

complex_sys=Equilibration/Family_${family}/${model}/systems/complex
solvent_sys=Equilibration/Family_${family}/${model}/systems/solvent

#############################################
##             COMPLEX FEP                 ##
#############################################
source /home/markelgi/MD_Programs/gromacs-2024.4/build/scripts/GMXRC
num_complex_lambdas=$(more conf/fep_complex.mdp | grep vdw-lambdas | cut -d "=" -f 2 | wc -w)
num_solvent_lambdas=$(more conf/fep_solvent.mdp | grep vdw-lambdas | cut -d "=" -f 2 | wc -w)

complex_prod=Production/Family_${family}/${model}/production/complex
min_dir=${complex_prod}/minimization
mkdir -p ${min_dir}

# Minimization
gmx grompp -f conf/minimization_complex.mdp -c ${complex_sys}/restraint_system.gro -p ${complex_sys}/system.top -r ${complex_sys}/restraint_system.gro -n ${complex_sys}/index.ndx -o ${min_dir}/min.tpr
gmx mdrun -v -deffnm ${min_dir}/min -ntomp 20 -ntmpi 1

#for i in $(seq 0 $(( ${num_complex_lambdas} - 1 ))); do
for i in 0; do
    equil_dir=${complex_prod}/equilibration/lambda_${i}
    prod_dir=${complex_prod}/production/lambda_${i}
    mkdir -p ${equil_dir}
    mkdir -p ${prod_dir}

    # NPT Equilibration
    cp conf/npt_equil_complex.mdp ${equil_dir}/npt.mdp
    sed -i "s/^init_lambda_state .*/init_lambda_state        = ${i}/" "${equil_dir}/npt.mdp"

    gmx grompp -f ${equil_dir}/npt.mdp -c ${min_dir}/min.gro -p ${complex_sys}/system.top -r ${min_dir}/min.gro -n ${complex_sys}/index.ndx -o ${equil_dir}/npt.tpr
    gmx mdrun -v -deffnm ${equil_dir}/npt -ntomp 20 -ntmpi 1 -pin on

    # Production
    cp conf/fep_complex.mdp ${prod_dir}/fep.mdp
    sed -i "s/^init_lambda_state .*/init_lambda_state        = ${i}/" "${prod_dir}/fep.mdp"

    gmx grompp -f ${prod_dir}/fep.mdp -c ${equil_dir}/npt.gro -p ${complex_sys}/system.top -r ${equil_dir}/npt.gro -n ${complex_sys}/index.ndx -o ${prod_dir}/fep.tpr
    gmx mdrun -v -deffnm ${prod_dir}/fep -ntomp 20 -ntmpi 1 -pin on

done

#############################################
##             SOLVENT FEP                 ##
#############################################
solvent_prod=Production/Family_${family}/${model}/production/solution
min_dir=${solvent_prod}/minimization
mkdir -p ${min_dir}

# Minimization
gmx grompp -f conf/minimization_solvent.mdp -c ${solvent_sys}/system.gro -p ${solvent_sys}/system.top -r ${solvent_sys}/system.gro -n ${solvent_sys}/index.ndx -o ${min_dir}/min.tpr
gmx mdrun -v -deffnm ${min_dir}/min -ntomp 20 -ntmpi 1

#for i in $(seq 0 $(( ${num_solvent_lambdas} - 1 ))); do
for i in 0 20; do
    equil_dir=${solvent_prod}/equilibration/lambda_${i}
    prod_dir=${solvent_prod}/production/lambda_${i}
    mkdir -p ${equil_dir}
    mkdir -p ${prod_dir}

    # NVT Equilibration
    cp conf/nvt_equil_solvent.mdp ${equil_dir}/nvt.mdp
    sed -i "s/^init_lambda_state .*/init_lambda_state        = ${i}/" "${equil_dir}/nvt.mdp"

    gmx grompp -f ${equil_dir}/nvt.mdp -c ${min_dir}/min.gro -p ${solvent_sys}/system.top -r ${min_dir}/min.gro -n ${solvent_sys}/index.ndx -o ${equil_dir}/nvt.tpr
    gmx mdrun -v -deffnm ${equil_dir}/nvt -ntomp 20 -ntmpi 1 -pin on

    # NPT Equilibration
    cp conf/npt_equil_solvent.mdp ${equil_dir}/npt.mdp
    sed -i "s/^init_lambda_state .*/init_lambda_state        = ${i}/" "${equil_dir}/npt.mdp"

    gmx grompp -f ${equil_dir}/npt.mdp -c ${equil_dir}/nvt.gro -p ${solvent_sys}/system.top -r ${equil_dir}/nvt.gro -n ${solvent_sys}/index.ndx -o ${equil_dir}/npt.tpr
    gmx mdrun -v -deffnm ${equil_dir}/npt -ntomp 20 -ntmpi 1 -pin on

    # Production
    cp conf/fep_solvent.mdp ${prod_dir}/fep.mdp
    sed -i "s/^init_lambda_state .*/init_lambda_state        = ${i}/" "${prod_dir}/fep.mdp"

    gmx grompp -f ${prod_dir}/fep.mdp -c ${equil_dir}/npt.gro -p ${solvent_sys}/system.top -n ${solvent_sys}/index.ndx -o ${prod_dir}/fep.tpr
    gmx mdrun -v -deffnm ${prod_dir}/fep -ntomp 20 -ntmpi 1 -pin on

done

