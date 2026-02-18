#!/usr/bin/env python3
"""
Run AutoDock Vina docking of fluoroquinolones to M. abscessus gyrase
Using OpenBabel with proper rigid receptor settings
"""

from vina import Vina
from openbabel import openbabel as ob
import os
import subprocess

# Paths
BASE_DIR = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling"
DOCKING_DIR = f"{BASE_DIR}/docking"
RECEPTOR_PDB = f"{BASE_DIR}/output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"

# Binding site center (from setup_docking.py)
CENTER = [32.28, 11.17, 22.14]
BOX_SIZE = [25, 25, 25]

def pdb_to_rigid_pdbqt(input_pdb, output_pdbqt):
    """Convert PDB to rigid receptor PDBQT"""
    # Use OpenBabel with rigid receptor option
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")
    # Set option for rigid receptor (no torsions)
    obConversion.AddOption("r", ob.OBConversion.OUTOPTIONS)  
    obConversion.AddOption("x", ob.OBConversion.OUTOPTIONS)  # No hydrogens
    
    mol = ob.OBMol()
    obConversion.ReadFile(mol, input_pdb)
    
    # Don't add hydrogens for receptor to keep it simple
    obConversion.WriteFile(mol, output_pdbqt)
    print(f"Receptor converted: {output_pdbqt}")

def create_simple_receptor_pdbqt(input_pdb, output_pdbqt):
    """Create a simple rigid receptor PDBQT by direct conversion"""
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    with open(output_pdbqt, 'w') as f:
        for line in lines:
            if line.startswith('ATOM'):
                # Parse atom info
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                
                # Determine element/atom type
                element = atom_name[0] if atom_name else 'C'
                if atom_name.startswith(('1H', '2H', '3H')):
                    element = 'H'
                elif len(atom_name) > 1 and atom_name[1].isalpha() and atom_name[0] in 'CNOSH':
                    if atom_name[:2] in ['CL', 'BR', 'CA', 'MG', 'FE', 'ZN']:
                        element = atom_name[:2]
                    else:
                        element = atom_name[0]
                
                # Map to AutoDock atom types
                ad_type = element
                if element == 'C':
                    ad_type = 'C'
                elif element == 'N':
                    ad_type = 'N' if 'H' not in [l.strip()[12:16].strip() for l in lines if l.startswith('ATOM') and l[22:26] == line[22:26]] else 'NA'
                elif element == 'O':
                    ad_type = 'OA'
                elif element == 'S':
                    ad_type = 'SA'
                elif element == 'H':
                    ad_type = 'HD'
                
                # Simple partial charge (0.0 for now)
                charge = 0.000
                
                # Format PDBQT line
                new_line = f"{line[:54]}{charge:>8.3f} {ad_type:>2}\n"
                f.write(new_line)
            elif line.startswith(('TER', 'END')):
                f.write(line)
    
    print(f"Receptor saved to: {output_pdbqt}")

def sdf_to_pdbqt(input_file, output_file):
    """Convert SDF to PDBQT using OpenBabel"""
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "pdbqt")
    
    mol = ob.OBMol()
    obConversion.ReadFile(mol, input_file)
    
    # Add hydrogens
    mol.AddHydrogens()
    
    # Optimize geometry
    ff = ob.OBForceField.FindForceField("mmff94")
    if ff:
        ff.Setup(mol)
        ff.ConjugateGradients(500)
        ff.GetCoordinates(mol)
    
    obConversion.WriteFile(mol, output_file)
    print(f"Ligand saved to: {output_file}")

def run_docking(receptor_pdbqt, ligand_pdbqt, ligand_name):
    """Run Vina docking"""
    print(f"\n{'='*60}")
    print(f" DOCKING: {ligand_name}")
    print(f"{'='*60}")
    
    v = Vina(sf_name='vina')
    
    # Set receptor (positional argument)
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
    print("-" * 60)
    print(f"{'Pose':<6} {'Affinity (kcal/mol)':<20} {'RMSD l.b.':<12} {'RMSD u.b.':<12}")
    print("-" * 60)
    
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
    print("\n1. Preparing rigid receptor...")
    receptor_pdbqt = f"{DOCKING_DIR}/receptor.pdbqt"
    create_simple_receptor_pdbqt(RECEPTOR_PDB, receptor_pdbqt)
    
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
            ligand_pdbqt = f"{DOCKING_DIR}/{name}.pdbqt"
            try:
                sdf_to_pdbqt(sdf_file, ligand_pdbqt)
                score = run_docking(receptor_pdbqt, ligand_pdbqt, name)
                results[name] = score
            except Exception as e:
                print(f"Error with {name}: {e}")
                import traceback
                traceback.print_exc()
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
  
  Typical values for fluoroquinolone-gyrase binding:
    -10 to -12 kcal/mol : Expected range for good binders
  
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
