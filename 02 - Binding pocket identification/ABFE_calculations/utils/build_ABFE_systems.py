from pathlib import Path
from shutil import copy
import mdtraj as md
from openff.toolkit import Molecule
from openmm import app, unit
from openmmforcefields.generators import GAFFTemplateGenerator
import parmed as pm

def ndx_add_group(index_file, group_name, group_indices):
    
    with open(index_file, 'a') as ndx_file:
        
        ndx_file.write(f'\n[ {group_name} ]\n')
        
        n_lines = len(group_indices) // 15
        n_lines += 0 if len(group_indices) % 15 == 0 else 1
        
        for i in range(n_lines):
            ndx_file.write(' '.join([str(i).rjust(5) for i in group_indices[i*15:(i+1)*15]]) + '\n')

def add_ligand_and_binding_pocket(system, ligand_name, index_file):
    
    traj = md.load(system)
    protein = traj.top.select(f'protein and not resname {ligand_name}')
    ligand = traj.top.select(f'resname {ligand_name}')
    
    ndx_add_group(index_file, ligand_name, ligand + 1)
    
    binding_pocket = md.compute_neighbors(traj, cutoff=0.5, query_indices=ligand, haystack_indices=protein)[0] + 1
    ndx_add_group(index_file, 'Binding_Pocket', binding_pocket)

def prepare_complex_directory(gmx_dir, new_directory):
        
    gro_file = Path(gmx_dir, 'step5_input.gro')
    top_file = Path(gmx_dir, 'topol.top')
    index_file = Path(gmx_dir, 'index.ndx')
    topologies = Path(gmx_dir, 'toppar')
    
    copy(gro_file, Path(new_directory, 'system.gro'))
    copy(top_file, Path(new_directory, 'system.top'))
    copy(index_file, Path(new_directory, 'index.ndx'))
    
    toppor_dir = Path(new_directory, 'toppar')
    toppor_dir.mkdir(exist_ok=True)
    
    for top in topologies.glob('*'):
        copy(top, Path(toppor_dir, top.name))

    add_ligand_and_binding_pocket(gro_file.as_posix(), 'UNK', Path(new_directory, 'index.ndx').as_posix())

def _add_ligand_restraints(system, topology, restraint_file_dir):
    
    traj = md.load(system)
    ligand = traj.top.select('resname UNK') + 1
    
    ligand_restraint = Path(restraint_file_dir, 'ligres.itp')
    with open(ligand_restraint, 'w') as lig_rest:
        
        lig_rest.write('[ position_restraints ]\n')
        
        for idx in ligand:
            
            lig_rest.write(str(idx).rjust(6) + '     1  1000  1000  1000\n')
    
    with open(topology, 'r') as old_top:
        lines = old_top.readlines()
        
    with open(topology, 'w') as new_top:
        
        write_ligand_rest = False
        
        for n, line in enumerate(lines):
            
            if line.startswith('[ moleculetype ]'):
                
                if write_ligand_rest:
                    new_top.write('#ifdef LIGRES\n')
                    new_top.write(f'#include "{ligand_restraint.absolute().as_posix()}"\n')
                    new_top.write('#endif\n\n')
                    write_ligand_rest = False
                
                if lines[n+1].startswith('UNK') or lines[n+2].startswith('UNK'):
                    write_ligand_rest = True
            
            new_top.write(line)

def prepare_ligand_solution(ligand_file, output_directory):
    
    ligand = Molecule.from_file(ligand_file.as_posix())
    
    forcefield = app.ForceField('amber/tip3p_standard.xml')
    gaff = GAFFTemplateGenerator(ligand, forcefield='gaff-2.11')
    forcefield.registerTemplateGenerator(gaff.generator)
    model = app.Modeller(ligand.to_topology().to_openmm(), ligand.to_topology().get_positions().to_openmm())
    model.addSolvent(forcefield, padding=2 * unit.nanometer, neutralize=True, ionicStrength=0.15 * unit.molar, positiveIon='K+')
    
    system = forcefield.createSystem(model.topology, constraints=None, rigidWater=False)
    
    solvated_structure = pm.openmm.load_topology(model.topology, system, model.positions)
    solvated_structure.save(f'{output_directory}/system.top')
    solvated_structure.save(f'{output_directory}/system.gro')
    
    _add_ligand_restraints(f'{output_directory}/system.gro', f'{output_directory}/system.top', output_directory)

def main(gmx_dir, ligand_file, output_dir):
    
    complex_dir = Path(output_dir, 'complex')
    complex_dir.mkdir(exist_ok=True, parents=True)
    ligand_dir = Path(output_dir, 'solution')
    ligand_dir.mkdir(exist_ok=True)
    
    prepare_complex_directory(gmx_dir, complex_dir)
    prepare_ligand_solution(ligand_file, ligand_dir)

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-d', '--complex_dir', type=str, required=True)
    parser.add_argument('-l', '--ligand_file', type=str, required=True)
    parser.add_argument('-o', '--output_dir', type=str, required=True)
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    ligand_file = Path(args.ligand_file)
    complex_dir = Path(args.complex_dir)
    
    main(complex_dir, ligand_file, output_dir)
    