#!/usr/bin/env python3
"""
Run AutoDock Vina docking - Simple version using command line
"""
import subprocess
import os

BASE_DIR = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling"
DOCKING_DIR = f"{BASE_DIR}/docking"
RECEPTOR_PDB = f"{BASE_DIR}/output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"
VENV_PYTHON = f"{BASE_DIR}/.venv/bin/python"

# Binding site center
CENTER = [32.28, 11.17, 22.14]
BOX_SIZE = [25, 25, 25]

def create_proper_pdbqt(input_pdb, output_pdbqt):
    """Create proper PDBQT format for receptor"""
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    with open(output_pdbqt, 'w') as f:
        atom_num = 0
        for line in lines:
            if line.startswith('ATOM'):
                atom_num += 1
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                # Determine atom type for AutoDock
                element = atom_name[0]
                if atom_name.startswith(('1H', '2H', '3H', 'H')):
                    ad_type = 'HD'  # Hydrogen donor
                elif element == 'C':
                    ad_type = 'C'
                elif element == 'N':
                    # Check if it's a hydrogen bond donor (has H attached)
                    ad_type = 'NA'  # Nitrogen acceptor
                elif element == 'O':
                    ad_type = 'OA'  # Oxygen acceptor  
                elif element == 'S':
                    ad_type = 'SA'  # Sulfur acceptor
                else:
                    ad_type = element
                
                # Proper PDBQT format with partial charge
                # Columns: 1-6, 7-11, 13-16, 17, 18-20, 22, 23-26, 27, 31-38, 39-46, 47-54, 55-60, 61-66, 67-76, 77-78
                f.write(f"ATOM  {atom_num:5d} {atom_name:4s} {res_name:3s} {chain}{res_num}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +0.000 {ad_type:2s}\n")
            elif line.startswith('TER'):
                f.write(f"TER   {atom_num+1:5d}      {line[17:26]}\n")
            elif line.startswith('END'):
                break
        f.write("END\n")
    print(f"Created receptor PDBQT: {output_pdbqt}")

def run_vina_docking():
    """Run docking using Vina Python API with proper format"""
    
    print("="*60)
    print(" AUTODOCK VINA DOCKING")
    print(" M. abscessus DNA Gyrase + Fluoroquinolones")  
    print("="*60)
    print(f"\nBinding site: center={CENTER}, size={BOX_SIZE}")
    
    # Create receptor PDBQT
    receptor_pdbqt = f"{DOCKING_DIR}/receptor_fixed.pdbqt"
    create_proper_pdbqt(RECEPTOR_PDB, receptor_pdbqt)
    
    # Check if format is correct
    with open(receptor_pdbqt) as f:
        first_line = f.readline()
    print(f"First receptor line: {first_line.strip()}")
    
    # Now try docking with the properly formatted files
    from vina import Vina
    from openbabel import openbabel as ob
    
    ligands = ['moxifloxacin', 'ciprofloxacin', 'levofloxacin']
    results = {}
    
    for lig_name in ligands:
        sdf_file = f"{DOCKING_DIR}/ligands/{lig_name}.sdf"
        if not os.path.exists(sdf_file):
            print(f"Skipping {lig_name} - SDF not found")
            continue
            
        print(f"\n{'='*60}")
        print(f" Processing: {lig_name}")
        print('='*60)
        
        # Convert ligand with OpenBabel
        lig_pdbqt = f"{DOCKING_DIR}/{lig_name}_lig.pdbqt"
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "pdbqt")
        mol = ob.OBMol()
        obConversion.ReadFile(mol, sdf_file)
        mol.AddHydrogens()
        
        # Optimize
        ff = ob.OBForceField.FindForceField("mmff94")
        if ff:
            ff.Setup(mol)
            ff.ConjugateGradients(500)
            ff.GetCoordinates(mol)
        
        obConversion.WriteFile(mol, lig_pdbqt)
        print(f"Ligand prepared: {lig_pdbqt}")
        
        try:
            v = Vina(sf_name='vina')
            v.set_receptor(receptor_pdbqt)
            v.set_ligand_from_file(lig_pdbqt)
            v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
            
            print("Running docking...")
            v.dock(exhaustiveness=16, n_poses=5)
            
            energies = v.energies()
            
            print(f"\nDocking Results for {lig_name}:")
            print("-"*50)
            for i, e in enumerate(energies):
                print(f"Pose {i+1}: {e[0]:.2f} kcal/mol (RMSD: {e[1]:.1f}/{e[2]:.1f})")
            
            # Save poses
            out_pdbqt = f"{DOCKING_DIR}/docked_{lig_name}.pdbqt"
            v.write_poses(out_pdbqt, n_poses=5, overwrite=True)
            
            results[lig_name] = energies[0][0]
            
        except Exception as e:
            print(f"Error docking {lig_name}: {e}")
            results[lig_name] = None
    
    # Print summary
    print("\n" + "="*60)
    print(" DOCKING SCORE SUMMARY")
    print("="*60)
    print(f"{'Ligand':<20} {'Score (kcal/mol)':<20}")
    print("-"*40)
    for name, score in results.items():
        if score is not None:
            print(f"{name:<20} {score:<20.2f}")
        else:
            print(f"{name:<20} {'FAILED'}")
            
    return results

if __name__ == "__main__":
    run_vina_docking()
