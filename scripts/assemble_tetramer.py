#!/usr/bin/env python3
"""
Assemble M. abscessus DNA gyrase heterotetramer (GyrA₂GyrB₂) with DNA and Mg²⁺
Using relaxed homology models and DNA/ions from 5BS8 template
"""
import os
import numpy as np

def read_pdb(filename):
    """Read PDB file and return list of lines"""
    with open(filename, 'r') as f:
        return [line for line in f if line.startswith(('ATOM', 'HETATM', 'TER'))]

def get_atoms_by_chain(pdb_lines, chain_id):
    """Extract atoms for a specific chain"""
    return [line for line in pdb_lines if len(line) > 21 and line[21] == chain_id]

def change_chain_id(pdb_lines, new_chain_id):
    """Change chain ID in PDB lines"""
    modified = []
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM', 'TER')):
            modified.append(line[:21] + new_chain_id + line[22:])
        else:
            modified.append(line)
    return modified

def get_ca_coords(pdb_lines):
    """Extract CA atom coordinates"""
    coords = []
    for line in pdb_lines:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
    return np.array(coords)

def extract_transformation(template_chain_a, template_chain_c):
    """
    Calculate transformation matrix from chain A to chain C
    This gives us the C2 symmetry operation
    """
    ca_A = get_ca_coords(template_chain_a)
    ca_C = get_ca_coords(template_chain_c)
    
    # Use overlapping residues (they may have different lengths due to gaps)
    min_len = min(len(ca_A), len(ca_C))
    ca_A = ca_A[:min_len]
    ca_C = ca_C[:min_len]
    
    # Center both sets
    center_A = np.mean(ca_A, axis=0)
    center_C = np.mean(ca_C, axis=0)
    
    ca_A_centered = ca_A - center_A
    ca_C_centered = ca_C - center_C
    
    # Use Kabsch algorithm to find rotation
    H = ca_A_centered.T @ ca_C_centered
    U, S, Vt = np.linalg.svd(H)
    
    # Ensure proper rotation (not reflection)
    d = np.linalg.det(Vt.T @ U.T)
    if d < 0:
        Vt[-1, :] *= -1
    
    R = Vt.T @ U.T
    t = center_C - R @ center_A
    
    return R, t

def apply_transformation(pdb_lines, R, t):
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

def renumber_atoms(pdb_lines, start_num=1):
    """Renumber atoms sequentially"""
    renumbered = []
    atom_num = start_num
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            renumbered.append(line[:6] + f"{atom_num:5d}" + line[11:])
            atom_num += 1
        else:
            renumbered.append(line)
    return renumbered, atom_num

def main():
    # Paths
    base_dir = '/Users/nb/Desktop/rosetta-cm'
    template_file = os.path.join(base_dir, 'input/templates/5bs8.pdb')
    gyrA_file = os.path.join(base_dir, 'output/mabs_gyrA_new_relaxed_mabs_gyrA_threaded_new_0001.pdb')
    gyrB_file = os.path.join(base_dir, 'output/mabs_gyrB_new_relaxed_mabs_gyrB_threaded_new_0001.pdb')
    output_file = os.path.join(base_dir, 'output/mabs_gyrase_tetramer_dna_mg.pdb')
    
    print("Reading input files...")
    template = read_pdb(template_file)
    gyrA = read_pdb(gyrA_file)
    gyrB = read_pdb(gyrB_file)
    
    # Get template chains for transformation calculation
    template_chain_A = get_atoms_by_chain(template, 'A')
    template_chain_B = get_atoms_by_chain(template, 'B')
    template_chain_C = get_atoms_by_chain(template, 'C')
    template_chain_D = get_atoms_by_chain(template, 'D')
    
    print(f"Template chain A: {len(template_chain_A)} lines")
    print(f"Template chain C: {len(template_chain_C)} lines")
    
    # Calculate transformation from chain A to chain C
    print("Calculating symmetry transformation...")
    R, t = extract_transformation(template_chain_A, template_chain_C)
    print(f"Rotation matrix:\n{R}")
    print(f"Translation vector: {t}")
    
    # Create GyrA chain A (original)
    gyrA_chainA = change_chain_id(gyrA, 'A')
    
    # Create GyrA chain C (symmetry copy)
    gyrA_chainC_raw = apply_transformation(gyrA, R, t)
    gyrA_chainC = change_chain_id(gyrA_chainC_raw, 'C')
    
    # Create GyrB chain B (original)
    gyrB_chainB = change_chain_id(gyrB, 'B')
    
    # Create GyrB chain D (symmetry copy)
    gyrB_chainD_raw = apply_transformation(gyrB, R, t)
    gyrB_chainD = change_chain_id(gyrB_chainD_raw, 'D')
    
    # Extract DNA chains (E, F, G, H) from template
    print("Extracting DNA and Mg²⁺ from template...")
    dna_chains = []
    for chain_id in ['E', 'F', 'G', 'H']:
        dna_chain = get_atoms_by_chain(template, chain_id)
        if dna_chain:
            dna_chains.extend(dna_chain)
            print(f"  DNA chain {chain_id}: {len(dna_chain)} atoms")
    
    # Extract Mg²⁺ ions
    mg_ions = [line for line in template if 'MG    MG' in line or ' MG ' in line[17:20]]
    print(f"  Mg²⁺ ions: {len(mg_ions)}")
    
    # Assemble the full complex
    print("Assembling tetramer complex...")
    complex_lines = []
    
    # Header
    complex_lines.append("REMARK   M. abscessus DNA gyrase heterotetramer with DNA and Mg2+\n")
    complex_lines.append("REMARK   GyrA: UniProt B1ME58, GyrB: UniProt B1ME45\n")
    complex_lines.append("REMARK   Template: PDB 5BS8 (M. tuberculosis DNA gyrase)\n")
    complex_lines.append("REMARK   Chain A: GyrA subunit 1\n")
    complex_lines.append("REMARK   Chain B: GyrB subunit 1\n")
    complex_lines.append("REMARK   Chain C: GyrA subunit 2 (symmetry copy)\n")
    complex_lines.append("REMARK   Chain D: GyrB subunit 2 (symmetry copy)\n")
    complex_lines.append("REMARK   Chains E-H: DNA from template\n")
    complex_lines.append("REMARK   Mg2+ ions from template\n")
    
    # Add protein chains
    complex_lines.extend(gyrA_chainA)
    complex_lines.append("TER\n")
    complex_lines.extend(gyrB_chainB)
    complex_lines.append("TER\n")
    complex_lines.extend(gyrA_chainC)
    complex_lines.append("TER\n")
    complex_lines.extend(gyrB_chainD)
    complex_lines.append("TER\n")
    
    # Add DNA chains
    complex_lines.extend(dna_chains)
    complex_lines.append("TER\n")
    
    # Add Mg²⁺ ions
    complex_lines.extend(mg_ions)
    
    complex_lines.append("END\n")
    
    # Renumber atoms
    print("Renumbering atoms...")
    final_lines = []
    atom_num = 1
    for line in complex_lines:
        if line.startswith(('ATOM', 'HETATM')):
            final_lines.append(line[:6] + f"{atom_num:5d}" + line[11:])
            atom_num += 1
        else:
            final_lines.append(line)
    
    # Write output
    print(f"Writing output to {output_file}...")
    with open(output_file, 'w') as f:
        f.writelines(final_lines)
    
    # Count atoms per chain in output
    print("\nFinal complex composition:")
    chain_counts = {}
    for line in final_lines:
        if line.startswith(('ATOM', 'HETATM')) and len(line) > 21:
            chain = line[21]
            chain_counts[chain] = chain_counts.get(chain, 0) + 1
    
    for chain, count in sorted(chain_counts.items()):
        print(f"  Chain {chain}: {count} atoms")
    
    print(f"\nTotal atoms: {sum(chain_counts.values())}")
    print(f"Output file: {output_file}")
    print("\nTetramer assembly complete!")

if __name__ == '__main__':
    main()
