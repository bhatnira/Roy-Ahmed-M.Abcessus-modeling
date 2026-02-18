#!/usr/bin/env python3
"""Dock fluoroquinolones to the template (M. tuberculosis) gyrase for comparison"""

import os
from vina import Vina

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_PDB = os.path.join(BASE_DIR, "../input/templates/5bs8_protein_only.pdb")
LIGAND_DIR = os.path.join(BASE_DIR, "ligands")

# Binding site center (same region as M. abscessus - coordinates from template)
CENTER = [32.28, 11.17, 22.14]
BOX_SIZE = [25, 25, 25]

def pdb_to_pdbqt(pdb_file, pdbqt_file):
    """Convert PDB to PDBQT format for receptor"""
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    with open(pdbqt_file, 'w') as f:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Only keep protein atoms
                resname = line[17:20].strip()
                if resname in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                              'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                              'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    atom_name = line[12:16].strip()
                    # Determine atom type
                    if atom_name.startswith('C'):
                        atom_type = 'C'
                    elif atom_name.startswith('N'):
                        atom_type = 'N'
                    elif atom_name.startswith('O'):
                        atom_type = 'O'
                    elif atom_name.startswith('S'):
                        atom_type = 'S'
                    elif atom_name.startswith('H'):
                        atom_type = 'H'
                    else:
                        atom_type = atom_name[0]
                    
                    # Format PDBQT line
                    pdbqt_line = line[:54].ljust(54) + "  0.00  0.00"
                    pdbqt_line = pdbqt_line.ljust(77) + f" {atom_type}\n"
                    f.write(pdbqt_line)
            elif line.startswith('END'):
                break

def sdf_to_pdbqt(sdf_file, pdbqt_file):
    """Convert SDF to PDBQT using OpenBabel"""
    import subprocess
    obabel_path = os.path.join(BASE_DIR, "../.venv/bin/obabel")
    cmd = [obabel_path, sdf_file, '-O', pdbqt_file, '--gen3d']
    subprocess.run(cmd, capture_output=True)

def main():
    print("="*50)
    print("DOCKING TO TEMPLATE (M. tuberculosis - 5BS8)")
    print("="*50)
    
    # Convert template to PDBQT
    receptor_pdbqt = os.path.join(BASE_DIR, "template_receptor.pdbqt")
    print(f"\nConverting template to PDBQT...")
    pdb_to_pdbqt(TEMPLATE_PDB, receptor_pdbqt)
    print(f"  Saved: {receptor_pdbqt}")
    
    # Initialize Vina
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
    
    # Ligands to dock
    ligands = ['moxifloxacin', 'ciprofloxacin', 'levofloxacin']
    
    results = {}
    
    for ligand_name in ligands:
        print(f"\n--- Docking {ligand_name} to template ---")
        
        sdf_file = os.path.join(LIGAND_DIR, f"{ligand_name}.sdf")
        pdbqt_in = os.path.join(BASE_DIR, f"{ligand_name}_temp.pdbqt")
        pdbqt_out = os.path.join(BASE_DIR, f"template_docked_{ligand_name}.pdbqt")
        
        # Convert ligand
        sdf_to_pdbqt(sdf_file, pdbqt_in)
        
        # Dock
        v.set_ligand_from_file(pdbqt_in)
        v.dock(exhaustiveness=16, n_poses=5)
        
        # Get score
        energies = v.energies()
        best_score = energies[0][0]
        results[ligand_name] = best_score
        
        print(f"  Best score: {best_score:.2f} kcal/mol")
        
        # Save
        v.write_poses(pdbqt_out, n_poses=1, overwrite=True)
        print(f"  Saved: {pdbqt_out}")
        
        # Cleanup
        if os.path.exists(pdbqt_in):
            os.remove(pdbqt_in)
    
    # Summary
    print("\n" + "="*50)
    print("TEMPLATE DOCKING SCORE SUMMARY")
    print("="*50)
    print(f"{'Ligand':<20} {'Score (kcal/mol)':<20}")
    print("-"*40)
    for ligand, score in results.items():
        print(f"{ligand:<20} {score:<20.2f}")
    
    return results

if __name__ == "__main__":
    main()
