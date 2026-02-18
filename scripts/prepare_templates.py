#!/usr/bin/env python3
"""
Prepare template structures for RosettaCM dimeric homology modeling.

This script:
1. Cleans PDB files (removes waters, ligands, alternate conformations)
2. Renumbers residues
3. Extracts sequences for alignment
4. Prepares templates for threading

Usage:
    python prepare_templates.py --templates input/templates/ --output input/prepared_templates/
"""

import os
import argparse
from pathlib import Path

def clean_pdb(pdb_file, output_file, keep_chains=None):
    """
    Clean a PDB file for use as a template.
    
    Args:
        pdb_file: Input PDB file path
        output_file: Output cleaned PDB file path
        keep_chains: List of chain IDs to keep (None = keep all)
    """
    clean_lines = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            # Keep only ATOM records
            if line.startswith('ATOM'):
                # Skip alternate conformations (keep only 'A' or ' ')
                alt_loc = line[16]
                if alt_loc not in [' ', 'A']:
                    continue
                
                # Filter by chain if specified
                chain_id = line[21]
                if keep_chains and chain_id not in keep_chains:
                    continue
                
                # Clean the alternate location indicator
                clean_line = line[:16] + ' ' + line[17:]
                clean_lines.append(clean_line)
            
            elif line.startswith('TER'):
                clean_lines.append(line)
    
    clean_lines.append('END\n')
    
    with open(output_file, 'w') as f:
        f.writelines(clean_lines)
    
    print(f"Cleaned: {pdb_file} -> {output_file}")

def extract_sequence_from_pdb(pdb_file):
    """
    Extract amino acid sequence from a PDB file.
    
    Returns:
        Dictionary mapping chain ID to sequence
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'MSE': 'M',  # Selenomethionine
    }
    
    sequences = {}
    last_residue = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain_id = line[21]
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                
                if chain_id not in sequences:
                    sequences[chain_id] = []
                    last_residue[chain_id] = None
                
                if res_num != last_residue[chain_id]:
                    if res_name in three_to_one:
                        sequences[chain_id].append(three_to_one[res_name])
                    last_residue[chain_id] = res_num
    
    return {chain: ''.join(seq) for chain, seq in sequences.items()}

def write_template_fasta(sequences, template_name, output_file):
    """Write template sequences to FASTA format."""
    with open(output_file, 'w') as f:
        for chain_id, sequence in sequences.items():
            f.write(f">{template_name}_Chain{chain_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')

def main():
    parser = argparse.ArgumentParser(description='Prepare templates for RosettaCM')
    parser.add_argument('--templates', default='input/templates/',
                        help='Directory containing template PDB files')
    parser.add_argument('--output', default='input/prepared_templates/',
                        help='Output directory for cleaned templates')
    parser.add_argument('--chains', nargs='+', default=None,
                        help='Chains to keep (e.g., A B)')
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    os.makedirs('input/alignments', exist_ok=True)
    
    # Process each template
    template_dir = Path(args.templates)
    
    if not template_dir.exists():
        print(f"Template directory not found: {template_dir}")
        print("Please create the directory and add template PDB files.")
        return
    
    pdb_files = list(template_dir.glob('*.pdb'))
    
    if not pdb_files:
        print(f"No PDB files found in {template_dir}")
        print("Please add template PDB files to the directory.")
        return
    
    all_sequences = {}
    
    for pdb_file in pdb_files:
        template_name = pdb_file.stem
        output_pdb = Path(args.output) / f"{template_name}_clean.pdb"
        
        # Clean the PDB
        clean_pdb(pdb_file, output_pdb, keep_chains=args.chains)
        
        # Extract sequences
        sequences = extract_sequence_from_pdb(output_pdb)
        all_sequences[template_name] = sequences
        
        # Write template FASTA
        fasta_output = Path('input/alignments') / f"{template_name}.fasta"
        write_template_fasta(sequences, template_name, fasta_output)
        print(f"Wrote sequences to: {fasta_output}")
    
    print("\n" + "="*50)
    print("Template preparation complete!")
    print("="*50)
    print("\nNext steps:")
    print("1. Generate alignments: python scripts/generate_alignment.py")
    print("2. Review alignments in input/alignments/")

if __name__ == '__main__':
    main()
