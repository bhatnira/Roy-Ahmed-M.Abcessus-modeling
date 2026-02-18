#!/usr/bin/env python3
"""
Structural alignment of M. abscessus gyrase model with M. tuberculosis template (5BS8)
Aligns the model to the template and outputs a combined PDB for visualization
"""
import numpy as np
import os

def read_pdb_atoms(filename):
    """Read ATOM records from PDB file"""
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append(line)
    return atoms

def get_ca_coords_and_resids(pdb_lines, chain_id=None):
    """Extract CA coordinates and residue IDs"""
    coords = []
    resids = []
    for line in pdb_lines:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            if chain_id is None or line[21] == chain_id:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resid = int(line[22:26].strip())
                coords.append([x, y, z])
                resids.append(resid)
    return np.array(coords), resids

def kabsch_align(P, Q):
    """
    Kabsch algorithm to find optimal rotation matrix
    P: mobile points (model)
    Q: reference points (template)
    Returns rotation matrix R and translation t such that P' = R @ P + t aligns to Q
    """
    # Center both point sets
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    
    # Compute covariance matrix
    H = P_centered.T @ Q_centered
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Ensure proper rotation (not reflection)
    d = np.linalg.det(Vt.T @ U.T)
    if d < 0:
        Vt[-1, :] *= -1
    
    # Rotation matrix
    R = Vt.T @ U.T
    
    # Translation
    t = centroid_Q - R @ centroid_P
    
    return R, t

def apply_transform(pdb_lines, R, t):
    """Apply rotation and translation to PDB coordinates"""
    transformed = []
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            coord = np.array([x, y, z])
            new_coord = R @ coord + t
            
            new_line = (line[:30] + 
                       f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}" + 
                       line[54:])
            transformed.append(new_line)
        else:
            transformed.append(line)
    return transformed

def calculate_rmsd(P, Q):
    """Calculate RMSD between two sets of points"""
    diff = P - Q
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def change_chain_id(pdb_lines, old_chain, new_chain):
    """Change chain ID in PDB lines"""
    modified = []
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM', 'TER')) and len(line) > 21:
            if line[21] == old_chain:
                line = line[:21] + new_chain + line[22:]
        modified.append(line)
    return modified

def main():
    base_dir = '/Users/nb/Desktop/rosetta-cm'
    
    # Input files
    model_file = os.path.join(base_dir, 'output/tetramer/mabs_gyrase_tetramer_protein_only.pdb')
    template_file = os.path.join(base_dir, 'input/templates/5bs8.pdb')
    
    # Output file
    output_file = os.path.join(base_dir, 'output/tetramer/structural_alignment_model_vs_template.pdb')
    
    print("=" * 60)
    print("Structural Alignment: M. abscessus Model vs M. tuberculosis Template")
    print("=" * 60)
    
    # Read PDB files
    print("\nReading files...")
    model_atoms = read_pdb_atoms(model_file)
    template_atoms = read_pdb_atoms(template_file)
    
    print(f"  Model atoms: {len(model_atoms)}")
    print(f"  Template atoms: {len(template_atoms)}")
    
    # Get CA coordinates for chain A (GyrA)
    model_ca_A, model_resids_A = get_ca_coords_and_resids(model_atoms, 'A')
    template_ca_A, template_resids_A = get_ca_coords_and_resids(template_atoms, 'A')
    
    print(f"\nChain A (GyrA):")
    print(f"  Model CA atoms: {len(model_ca_A)} (residues {min(model_resids_A)}-{max(model_resids_A)})")
    print(f"  Template CA atoms: {len(template_ca_A)} (residues {min(template_resids_A)}-{max(template_resids_A)})")
    
    # Find common residues for alignment
    common_resids = sorted(set(model_resids_A) & set(template_resids_A))
    print(f"  Common residues: {len(common_resids)}")
    
    # Extract matching CA coordinates
    model_common = []
    template_common = []
    for resid in common_resids:
        if resid in model_resids_A and resid in template_resids_A:
            model_idx = model_resids_A.index(resid)
            template_idx = template_resids_A.index(resid)
            model_common.append(model_ca_A[model_idx])
            template_common.append(template_ca_A[template_idx])
    
    model_common = np.array(model_common)
    template_common = np.array(template_common)
    
    print(f"  Aligned residues: {len(model_common)}")
    
    # Calculate initial RMSD
    initial_rmsd = calculate_rmsd(model_common, template_common)
    print(f"\n  Initial RMSD: {initial_rmsd:.2f} Å")
    
    # Perform Kabsch alignment
    print("\nPerforming structural alignment (Kabsch algorithm)...")
    R, t = kabsch_align(model_common, template_common)
    
    # Apply transformation to model
    aligned_model = apply_transform(model_atoms, R, t)
    
    # Recalculate RMSD after alignment
    aligned_model_ca, _ = get_ca_coords_and_resids(aligned_model, 'A')
    aligned_common = []
    for resid in common_resids:
        if resid in model_resids_A:
            idx = model_resids_A.index(resid)
            aligned_common.append(aligned_model_ca[idx])
    aligned_common = np.array(aligned_common)
    
    final_rmsd = calculate_rmsd(aligned_common, template_common)
    print(f"  Final RMSD (Chain A): {final_rmsd:.2f} Å")
    
    # Also check Chain B alignment
    model_ca_B, model_resids_B = get_ca_coords_and_resids(aligned_model, 'B')
    template_ca_B, template_resids_B = get_ca_coords_and_resids(template_atoms, 'B')
    
    common_resids_B = sorted(set(model_resids_B) & set(template_resids_B))
    if common_resids_B:
        model_common_B = []
        template_common_B = []
        for resid in common_resids_B:
            if resid in model_resids_B and resid in template_resids_B:
                model_idx = model_resids_B.index(resid)
                template_idx = template_resids_B.index(resid)
                model_common_B.append(model_ca_B[model_idx])
                template_common_B.append(template_ca_B[template_idx])
        
        if model_common_B:
            rmsd_B = calculate_rmsd(np.array(model_common_B), np.array(template_common_B))
            print(f"  Final RMSD (Chain B): {rmsd_B:.2f} Å")
    
    # Rename template chains for visualization
    # Template: A->E, B->F, C->G, D->H (to avoid chain ID conflicts)
    print("\nPreparing combined PDB for visualization...")
    
    # Filter template to protein chains only (A, B, C, D)
    template_protein = [line for line in template_atoms 
                       if line.startswith('ATOM') and line[21] in 'ABCD']
    
    # Change template chain IDs: A->W, B->X, C->Y, D->Z
    template_renamed = change_chain_id(template_protein, 'A', 'W')
    template_renamed = change_chain_id(template_renamed, 'B', 'X')
    template_renamed = change_chain_id(template_renamed, 'C', 'Y')
    template_renamed = change_chain_id(template_renamed, 'D', 'Z')
    
    # Write combined PDB
    print(f"\nWriting aligned structure to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("REMARK   Structural Alignment: M. abscessus Model vs M. tuberculosis Template\n")
        f.write("REMARK   \n")
        f.write("REMARK   MODEL (M. abscessus):\n")
        f.write("REMARK     Chain A: GyrA subunit 1 (486 residues)\n")
        f.write("REMARK     Chain B: GyrB subunit 1 (245 residues)\n")
        f.write("REMARK     Chain C: GyrA subunit 2\n")
        f.write("REMARK     Chain D: GyrB subunit 2\n")
        f.write("REMARK   \n")
        f.write("REMARK   TEMPLATE (M. tuberculosis, PDB 5BS8):\n")
        f.write("REMARK     Chain W: GyrA subunit 1 (original chain A)\n")
        f.write("REMARK     Chain X: GyrB subunit 1 (original chain B)\n")
        f.write("REMARK     Chain Y: GyrA subunit 2 (original chain C)\n")
        f.write("REMARK     Chain Z: GyrB subunit 2 (original chain D)\n")
        f.write("REMARK   \n")
        f.write(f"REMARK   Alignment RMSD (Chain A/W GyrA): {final_rmsd:.2f} Angstroms\n")
        f.write("REMARK   \n")
        f.write("REMARK   Visualization tip: Color by chain to distinguish model vs template\n")
        f.write("REMARK     - Chains A-D: M. abscessus homology model (cyan/blue)\n")
        f.write("REMARK     - Chains W-Z: M. tuberculosis template (orange/red)\n")
        f.write("REMARK   \n")
        
        # Write aligned model (chains A, B, C, D)
        atom_num = 1
        for line in aligned_model:
            if line.startswith('ATOM'):
                f.write(f"ATOM  {atom_num:5d}" + line[11:])
                atom_num += 1
        f.write("TER\n")
        
        # Write template (chains W, X, Y, Z)
        for line in template_renamed:
            if line.startswith('ATOM'):
                f.write(f"ATOM  {atom_num:5d}" + line[11:])
                atom_num += 1
        f.write("TER\n")
        f.write("END\n")
    
    # Count final atoms
    with open(output_file, 'r') as f:
        total_atoms = sum(1 for line in f if line.startswith('ATOM'))
    
    print(f"\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Output file: {output_file}")
    print(f"Total atoms: {total_atoms}")
    print(f"Alignment RMSD: {final_rmsd:.2f} Å")
    print("\nChain mapping for visualization:")
    print("  MODEL (M. abscessus):   Chains A, B, C, D")
    print("  TEMPLATE (5BS8):        Chains W, X, Y, Z")
    print("\nOpen in PyMOL, ChimeraX, or VMD to visualize the alignment.")
    print("=" * 60)

if __name__ == '__main__':
    main()
