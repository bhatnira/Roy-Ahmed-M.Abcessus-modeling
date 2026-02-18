#!/usr/bin/env python3
"""
Extract protein-only chains from the tetramer PDB
"""
import os

base_dir = '/Users/nb/Desktop/rosetta-cm'
input_file = os.path.join(base_dir, 'output/mabs_gyrase_tetramer_dna_mg.pdb')
output_file = os.path.join(base_dir, 'output/mabs_gyrase_tetramer_protein_only.pdb')

protein_chains = ['A', 'B', 'C', 'D']

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    atom_num = 1
    for line in f_in:
        if line.startswith('REMARK'):
            f_out.write(line)
        elif line.startswith('ATOM'):
            chain_id = line[21]
            if chain_id in protein_chains:
                # Renumber atoms
                new_line = line[:6] + f"{atom_num:5d}" + line[11:]
                f_out.write(new_line)
                atom_num += 1
        elif line.startswith('TER'):
            if atom_num > 1:  # Only write TER after atoms
                f_out.write(f"TER   {atom_num:5d}\n")
                atom_num += 1
        elif line.startswith('END'):
            f_out.write(line)

print(f"Protein-only PDB written to: {output_file}")

# Count atoms
with open(output_file, 'r') as f:
    atom_count = sum(1 for line in f if line.startswith('ATOM'))
print(f"Total protein atoms: {atom_count}")
