#!/usr/bin/env python3
"""
Run AutoDock Vina docking of fluoroquinolones to M. abscessus gyrase
"""

from vina import Vina
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Paths
BASE_DIR = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling"
DOCKING_DIR = f"{BASE_DIR}/docking"
RECEPTOR_PDB = f"{BASE_DIR}/output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"

# Binding site center (from setup_docking.py)
CENTER = [32.28, 11.17, 22.14]
BOX_SIZE = [25, 25, 25]

def prepare_receptor_simple(pdb_file, pdbqt_file):
    """Simple receptor preparation - just add charges"""
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    with open(pdbqt_file, 'w') as f:
        for line in lines:
            if line.startswith('ATOM'):
                # Add dummy charge and atom type
                atom_name = line[12:16].strip()
                element = atom_name[0]
                # Write PDBQT format with AutoDock atom type
                new_line = line[:54] + "  0.00" + line[60:76] + f"  {element}\n" if len(line) < 78 else line[:54] + "  0.00" + line[60:76] + f"  {element}\n"
                f.write(line[:66] + f"  0.00     {element:>2}\n")
            elif line.startswith(('TER', 'END')):
                f.write(line)
    print(f"Receptor saved to: {pdbqt_file}")

def prepare_ligand(sdf_file, name):
    """Prepare ligand using Meeko"""
    # Read molecule
    mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
    if mol is None:
        print(f"Error reading {sdf_file}")
        return None
    
    # Add hydrogens and generate 3D coordinates if needed
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Prepare for AutoDock using Meeko
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    
    # Write PDBQT
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
    
    pdbqt_file = f"{DOCKING_DIR}/{name}.pdbqt"
    with open(pdbqt_file, 'w') as f:
        f.write(pdbqt_string)
    
    print(f"Ligand saved to: {pdbqt_file}")
    return pdbqt_file

def run_docking(receptor_pdbqt, ligand_pdbqt, ligand_name):
    """Run Vina docking"""
    print(f"\n{'='*60}")
    print(f" DOCKING: {ligand_name}")
    print(f"{'='*60}")
    
    v = Vina(sf_name='vina')
    
    # Set receptor
    v.set_receptor(receptor_pdbqt)
    
    # Set ligand
    v.set_ligand_from_file(ligand_pdbqt)
    
    # Compute docking box
    v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
    
    # Run docking
    print("Running docking (exhaustiveness=32)...")
    v.dock(exhaustiveness=32, n_poses=10)
    
    # Get results
    energies = v.energies()
    
    print(f"\nResults for {ligand_name}:")
    print("-" * 50)
    print(f"{'Pose':<6} {'Affinity (kcal/mol)':<20} {'RMSD l.b.':<12} {'RMSD u.b.':<12}")
    print("-" * 50)
    
    for i, e in enumerate(energies):
        print(f"{i+1:<6} {e[0]:<20.2f} {e[1]:<12.2f} {e[2]:<12.2f}")
    
    # Save poses
    output_file = f"{DOCKING_DIR}/docked_{ligand_name}.pdbqt"
    v.write_poses(output_file, n_poses=10, overwrite=True)
    print(f"\nDocked poses saved to: {output_file}")
    
    return energies[0][0]  # Return best score

def main():
    os.makedirs(DOCKING_DIR, exist_ok=True)
    
    print("="*60)
    print(" AUTODOCK VINA DOCKING")
    print(" M. abscessus DNA Gyrase + Fluoroquinolones")
    print("="*60)
    print(f"\nBinding site center: {CENTER}")
    print(f"Box size: {BOX_SIZE} Angstroms")
    
    # Prepare receptor
    print("\n1. Preparing receptor...")
    receptor_pdbqt = f"{DOCKING_DIR}/receptor.pdbqt"
    prepare_receptor_simple(RECEPTOR_PDB, receptor_pdbqt)
    
    # Ligands to dock
    ligands = {
        'moxifloxacin': f"{DOCKING_DIR}/ligands/moxifloxacin.sdf",
        'ciprofloxacin': f"{DOCKING_DIR}/ligands/ciprofloxacin.sdf",
        'levofloxacin': f"{DOCKING_DIR}/ligands/levofloxacin.sdf"
    }
    
    # Prepare and dock each ligand
    results = {}
    print("\n2. Preparing and docking ligands...")
    
    for name, sdf_file in ligands.items():
        if os.path.exists(sdf_file):
            print(f"\nProcessing {name}...")
            ligand_pdbqt = prepare_ligand(sdf_file, name)
            if ligand_pdbqt:
                try:
                    score = run_docking(receptor_pdbqt, ligand_pdbqt, name)
                    results[name] = score
                except Exception as e:
                    print(f"Docking failed for {name}: {e}")
                    results[name] = None
    
    # Summary
    print("\n" + "="*60)
    print(" DOCKING RESULTS SUMMARY")
    print("="*60)
    print(f"\n{'Ligand':<20} {'Best Score (kcal/mol)':<25}")
    print("-"*45)
    for name, score in results.items():
        if score is not None:
            print(f"{name:<20} {score:<25.2f}")
        else:
            print(f"{name:<20} {'FAILED':<25}")
    
    print("\n" + "="*60)
    print(" INTERPRETATION")
    print("="*60)
    print("""
  More negative = stronger binding
  
  Typical interpretation:
    < -10 kcal/mol : Excellent binding
    -8 to -10      : Very good binding
    -6 to -8       : Good binding
    -4 to -6       : Moderate binding
    > -4           : Weak binding
  
  Note: These are predicted values and should be
  validated experimentally (e.g., MIC assays)
""")
    print("="*60)

if __name__ == "__main__":
    main()
