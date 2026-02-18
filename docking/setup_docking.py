#!/usr/bin/env python3
"""
Fluoroquinolone Docking Setup for M. abscessus Gyrase Model
Docks moxifloxacin/ciprofloxacin to the DNA cleavage site

This script:
1. Identifies the binding site from a reference structure (5CDQ with moxifloxacin)
2. Transfers binding site to M. abscessus model
3. Sets up docking with AutoDock Vina (if available) or provides coordinates for other tools
"""

import os
from pathlib import Path

# Base directories
BASE_DIR = Path("/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling")
DOCKING_DIR = BASE_DIR / "docking"
LIGAND_DIR = DOCKING_DIR / "ligands"

def mean_coords(coords):
    """Calculate mean of coordinates without numpy"""
    if not coords:
        return None
    n = len(coords)
    x = sum(c[0] for c in coords) / n
    y = sum(c[1] for c in coords) / n
    z = sum(c[2] for c in coords) / n
    return [x, y, z]

def read_pdb_coordinates(pdb_file, selection=None):
    """Read coordinates from PDB file"""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if selection is None or selection in line:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    return coords if coords else None

def get_binding_site_center(pdb_file, ligand_name="MFX"):
    """Get center of ligand binding site from reference structure"""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM') and ligand_name in line:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if coords:
        return mean_coords(coords)
    return None

def find_equivalent_binding_site(model_pdb, ref_pdb, ref_ligand_center):
    """
    Find equivalent binding site in model by structural alignment
    The binding site is at the GyrA-GyrB interface near the DNA cleavage site
    """
    # Key residues in fluoroquinolone binding (conserved across species):
    # GyrA: Ser83 (QRDR), Asp87 (QRDR), Tyr129 (catalytic)
    # GyrB: Interactions with DNA
    
    # For M. abscessus, we'll estimate based on the structural alignment
    # The binding site is typically at the DNA gate between GyrA subunits
    
    # Read model coordinates
    model_ca = []
    with open(model_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                chain = line[21]
                resi = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                model_ca.append({'chain': chain, 'resi': resi, 'coord': [x, y, z]})
    
    # Get the center of GyrA chains (A and C) around the QRDR region
    # QRDR typically spans residues 67-106 in E. coli numbering
    qrdr_coords = []
    for atom in model_ca:
        if atom['chain'] in ['A', 'C'] and 70 <= atom['resi'] <= 110:
            qrdr_coords.append(atom['coord'])
    
    if qrdr_coords:
        return mean_coords(qrdr_coords)
    
    # Fallback: center of the tetramer
    all_coords = [atom['coord'] for atom in model_ca]
    return mean_coords(all_coords)

def create_vina_config(center, size=(25, 25, 25)):
    """Create AutoDock Vina configuration file"""
    config = f"""# AutoDock Vina Configuration for M. abscessus Gyrase
# Fluoroquinolone binding site (DNA cleavage region)

receptor = receptor.pdbqt
ligand = ligand.pdbqt

center_x = {center[0]:.3f}
center_y = {center[1]:.3f}
center_z = {center[2]:.3f}

size_x = {size[0]}
size_y = {size[1]}
size_z = {size[2]}

exhaustiveness = 32
num_modes = 10
energy_range = 5

out = docked_poses.pdbqt
log = docking_log.txt
"""
    return config

def create_pymol_visualization():
    """Create PyMOL script to visualize binding site and docking"""
    
    pml = '''# PyMOL script to visualize fluoroquinolone binding site
# Run: pymol visualize_binding_site.pml

# Load structures
load ../output/tetramer/mabs_gyrase_tetramer_protein_only.pdb, protein
load reference_gyrase_fq.pdb, reference

# Load ligands
load ligands/moxifloxacin.sdf, moxifloxacin
load ligands/ciprofloxacin.sdf, ciprofloxacin

# Hide everything initially
hide everything

# Show protein as cartoon
show cartoon, protein
color cyan, protein and chain A+C
color lightblue, protein and chain B+D

# Superimpose reference to our model (GyrA chains)
super reference and chain A, protein and chain A

# Show reference ligand (MFX from 5CDQ)
show sticks, reference and resn MFX
show spheres, reference and resn MFX and name MG
color yellow, reference and resn MFX
set sphere_scale, 0.3

# QRDR region (Quinolone Resistance Determining Region)
# Residues ~83-87 in E. coli numbering, adjust for M. abscessus
select qrdr, protein and chain A and resi 80-95
color red, qrdr
show sticks, qrdr

# Catalytic tyrosine (Y129 equivalent)
select cat_tyr, protein and chain A and resi 120-135 and resn TYR
color magenta, cat_tyr
show sticks, cat_tyr

# Show DNA if present
show cartoon, reference and chain E+F
color orange, reference and chain E+F

# Binding site sphere
# pseudoatom binding_center, pos=[BINDING_CENTER_X, BINDING_CENTER_Y, BINDING_CENTER_Z]
# show spheres, binding_center
# set sphere_scale, 1.0, binding_center
# color green, binding_center

# Zoom to binding site
zoom reference and resn MFX, 15

# Set nice rendering
bg_color white
set ray_shadows, 0
set cartoon_fancy_helices, 1

# Labels
set label_size, 14
label qrdr and name CA, "QRDR"

print("Binding site visualization loaded")
print("Yellow sticks: Moxifloxacin from reference (PDB 5CDQ)")
print("Red: QRDR region (fluoroquinolone resistance mutations)")
print("Magenta: Catalytic tyrosine")
'''
    return pml

def main():
    print("="*70)
    print(" FLUOROQUINOLONE DOCKING SETUP")
    print(" M. abscessus DNA Gyrase Model")
    print("="*70)
    print()
    
    # Check files exist
    ref_pdb = DOCKING_DIR / "reference_gyrase_fq.pdb"
    model_pdb = BASE_DIR / "output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"
    
    if not ref_pdb.exists():
        print(f"ERROR: Reference structure not found: {ref_pdb}")
        return
    
    if not model_pdb.exists():
        print(f"ERROR: Model not found: {model_pdb}")
        return
    
    # Get binding site from reference
    print("1. Extracting binding site from reference (PDB 5CDQ with moxifloxacin)...")
    ref_center = get_binding_site_center(ref_pdb, "MFX")
    if ref_center is not None:
        print(f"   Reference binding site center: {ref_center}")
    
    # Find equivalent site in model
    print("\n2. Finding equivalent binding site in M. abscessus model...")
    model_center = find_equivalent_binding_site(model_pdb, ref_pdb, ref_center)
    print(f"   Model binding site center: {model_center}")
    
    # Create output files
    print("\n3. Creating docking configuration files...")
    
    # Vina config
    vina_config = create_vina_config(model_center)
    config_file = DOCKING_DIR / "vina_config.txt"
    with open(config_file, 'w') as f:
        f.write(vina_config)
    print(f"   Created: {config_file}")
    
    # PyMOL visualization
    pml_script = create_pymol_visualization()
    pml_file = DOCKING_DIR / "visualize_binding_site.pml"
    with open(pml_file, 'w') as f:
        f.write(pml_script)
    print(f"   Created: {pml_file}")
    
    # Print summary
    print()
    print("="*70)
    print(" BINDING SITE INFORMATION")
    print("="*70)
    print(f"""
  Center coordinates (for docking box):
    X: {model_center[0]:.2f}
    Y: {model_center[1]:.2f}
    Z: {model_center[2]:.2f}
  
  Box size: 25 x 25 x 25 Angstroms
  
  Key binding site residues (approximate):
    - QRDR: GyrA residues 80-95 (Ser83, Asp87 equivalents)
    - Catalytic Tyr: GyrA ~residue 129
    - Mg2+ coordination site
    - DNA intercalation site
""")
    
    print("="*70)
    print(" LIGANDS AVAILABLE")
    print("="*70)
    for lig in LIGAND_DIR.glob("*.sdf"):
        print(f"   - {lig.name}")
    
    print()
    print("="*70)
    print(" NEXT STEPS")
    print("="*70)
    print("""
  Option A: Visualize binding site in PyMOL:
    cd docking
    pymol visualize_binding_site.pml
    
  Option B: Run docking with AutoDock Vina (if installed):
    # 1. Convert receptor to PDBQT:
    obabel -ipdb ../output/tetramer/mabs_gyrase_tetramer_protein_only.pdb -opdbqt -O receptor.pdbqt
    
    # 2. Convert ligand to PDBQT:
    obabel -isdf ligands/moxifloxacin.sdf -opdbqt -O ligand.pdbqt --gen3D
    
    # 3. Run Vina:
    vina --config vina_config.txt
    
  Option C: Use online docking server:
    - SwissDock: http://www.swissdock.ch/
    - DockThor: https://dockthor.lncc.br/
    Upload: receptor PDB + ligand SDF
    Use binding site center: ({model_center[0]:.1f}, {model_center[1]:.1f}, {model_center[2]:.1f})
""")
    print("="*70)

if __name__ == "__main__":
    main()
