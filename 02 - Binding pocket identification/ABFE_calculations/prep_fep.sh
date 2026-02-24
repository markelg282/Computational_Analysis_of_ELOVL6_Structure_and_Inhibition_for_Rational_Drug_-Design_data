#!/bin/bash

family=$1
system=$2
sim_dir=Equilibration/Family_${family}/${system}
python utils/build_ABFE_systems.py -d Structures/Family_${family}/Systems/${system} -l Structures/Family_${family}/Ligands/${system}_Ligand.sdf -o ${sim_dir}/systems

complex_sys=${sim_dir}/systems/complex
solvent_sys=${sim_dir}/systems/solvent

###########################################
##         MEMBRANE EQUILIBRATION        ##
###########################################
equil_complex=${sim_dir}/membrane_equilibration
mkdir ${equil_complex}

# Complex Minimization
gmx grompp -f conf_files/minimization_complex.mdp -c ${complex_sys}/system.gro -r ${complex_sys}/system.gro -p ${complex_sys}/system.top -n ${complex_sys}/index.ndx -o ${equil_complex}/min.tpr
gmx mdrun -v -deffnm ${equil_complex}/min

# Complex Equibration
for i in $(seq 1 1 6); do

    j=$(( $i - 1 ))
    
    if [ $i == 1 ]; then
        struct=min
    else
        struct=equil_${j}
    fi

    gmx grompp -f conf_files/step_${i}_equil.mdp -c ${equil_complex}/${struct}.gro -r ${equil_complex}/${struct}.gro -p ${complex_sys}/system.top -n ${complex_sys}/index.ndx -o ${equil_complex}/equil_${i}.tpr
    gmx mdrun -v -deffnm ${equil_complex}/equil_${i}

done

# Production run for restraint assignment
gmx grompp -f conf_files/find_restraints.mdp -c ${equil_complex}/equil_6.gro -r ${equil_complex}/equil_6.gro -p ${complex_sys}/system.top -n ${complex_sys}/index.ndx -o ${equil_complex}/find_restr.tpr
gmx mdrun -v -deffnm ${equil_complex}/find_restr

restraints_dir=${sim_dir}/restraints
mkdir ${restraints_dir}
python utils/find_restraints.py --top ${equil_complex}/find_restr.tpr --traj ${equil_complex}/find_restr.xtc --ligand_selection "resname UNK" --force_constant 20 --outpath ${restraints_dir}

cp ${restraints_dir}/BoreschRestraint.top ${complex_sys}/toppar/BoreschRestraint.itp
cp ${restraints_dir}/ClosestRestraintFrame.gro ${complex_sys}/restraint_system.gro
echo "#include "toppar/BoreschRestraint.itp"" >> ${complex_sys}/system.top
