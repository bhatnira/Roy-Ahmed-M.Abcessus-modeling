#!/usr/bin/env python3
"""
Generate sequence alignments for RosettaCM dimeric homology modeling.

This script generates GRISHIN format alignments required by Rosetta's
partial_thread application.

For dimeric modeling, we need to align both chains of the target
to both chains of each template.

Usage:
    python generate_alignment.py --target input/target.fasta \
                                 --templates input/prepared_templates/ \
                                 --output input/alignments/
"""

import os
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

def read_fasta(fasta_file: str) -> Dict[str, str]:
    """Read sequences from a FASTA file."""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences

def simple_align(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Simple global sequence alignment using dynamic programming.
    For production use, consider using BioPython or external tools.
    
    Returns aligned sequences with gaps represented as '-'.
    """
    # Scoring parameters
    match_score = 2
    mismatch_score = -1
    gap_penalty = -2
    
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    score = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column
    for i in range(m + 1):
        score[i][0] = i * gap_penalty
    for j in range(n + 1):
        score[0][j] = j * gap_penalty
    
    # Fill scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score[i-1][j] + gap_penalty
            insert = score[i][j-1] + gap_penalty
            score[i][j] = max(match, delete, insert)
    
    # Traceback
    aligned1, aligned2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            current = score[i][j]
            diagonal = score[i-1][j-1]
            
            if seq1[i-1] == seq2[j-1]:
                expected_diagonal = diagonal + match_score
            else:
                expected_diagonal = diagonal + mismatch_score
            
            if current == expected_diagonal:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        
        if i > 0 and score[i][j] == score[i-1][j] + gap_penalty:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))

def write_grishin_alignment(target_name: str, target_seq: str,
                           template_name: str, template_seq: str,
                           aligned_target: str, aligned_template: str,
                           output_file: str):
    """
    Write alignment in GRISHIN format for Rosetta.
    
    Format:
    ## target template
    # 
    scores_from_program: 0
    0 aligned_target_sequence
    0 aligned_template_sequence
    --
    """
    with open(output_file, 'w') as f:
        f.write(f"## {target_name} {template_name}\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target}\n")
        f.write(f"0 {aligned_template}\n")
        f.write("--\n")

def create_dimer_alignment(target_seqs: Dict[str, str],
                          template_seqs: Dict[str, str],
                          template_name: str,
                          output_dir: str) -> str:
    """
    Create alignment file for dimeric protein.
    
    For a homodimer, we align:
    - Target Chain A to Template Chain A
    - Target Chain B to Template Chain B
    
    Returns the path to the combined alignment file.
    """
    # Get target sequences (assuming first two are Chain A and B)
    target_names = list(target_seqs.keys())
    target_seq_a = target_seqs[target_names[0]]
    target_seq_b = target_seqs[target_names[1]] if len(target_names) > 1 else target_seqs[target_names[0]]
    
    # Get template sequences
    template_names = list(template_seqs.keys())
    template_seq_a = template_seqs[template_names[0]]
    template_seq_b = template_seqs[template_names[1]] if len(template_names) > 1 else template_seqs[template_names[0]]
    
    # Align chains
    aligned_target_a, aligned_template_a = simple_align(target_seq_a, template_seq_a)
    aligned_target_b, aligned_template_b = simple_align(target_seq_b, template_seq_b)
    
    # Write combined alignment file
    output_file = os.path.join(output_dir, f"alignment_{template_name}.grishin")
    
    with open(output_file, 'w') as f:
        # Chain A alignment
        f.write(f"## Target_ChainA {template_name}_ChainA\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_a}\n")
        f.write(f"0 {aligned_template_a}\n")
        f.write("--\n")
        
        # Chain B alignment
        f.write(f"## Target_ChainB {template_name}_ChainB\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_b}\n")
        f.write(f"0 {aligned_template_b}\n")
        f.write("--\n")
    
    print(f"Created alignment: {output_file}")
    return output_file

def calculate_identity(aligned1: str, aligned2: str) -> float:
    """Calculate sequence identity from aligned sequences."""
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')
    length = max(len(aligned1.replace('-', '')), len(aligned2.replace('-', '')))
    return (matches / length * 100) if length > 0 else 0

def main():
    parser = argparse.ArgumentParser(description='Generate alignments for RosettaCM')
    parser.add_argument('--target', default='input/target.fasta',
                        help='Target sequence FASTA file')
    parser.add_argument('--templates', default='input/prepared_templates/',
                        help='Directory with prepared template PDB files')
    parser.add_argument('--output', default='input/alignments/',
                        help='Output directory for alignments')
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Read target sequence
    if not os.path.exists(args.target):
        print(f"Target FASTA not found: {args.target}")
        return
    
    target_seqs = read_fasta(args.target)
    print(f"Loaded {len(target_seqs)} target sequence(s)")
    
    # Find template FASTA files
    template_dir = Path(args.templates)
    alignment_dir = Path('input/alignments')
    
    template_fastas = list(alignment_dir.glob('*.fasta'))
    
    if not template_fastas:
        print("No template FASTA files found.")
        print("Please run prepare_templates.py first.")
        return
    
    alignment_files = []
    
    for template_fasta in template_fastas:
        template_name = template_fasta.stem
        template_seqs = read_fasta(template_fasta)
        
        print(f"\nProcessing template: {template_name}")
        print(f"  Template chains: {list(template_seqs.keys())}")
        
        # Create alignment
        alignment_file = create_dimer_alignment(
            target_seqs, template_seqs, template_name, args.output
        )
        alignment_files.append(alignment_file)
    
    # Write alignment list file
    list_file = os.path.join(args.output, 'alignment_list.txt')
    with open(list_file, 'w') as f:
        for af in alignment_files:
            f.write(os.path.basename(af) + '\n')
    
    print("\n" + "="*50)
    print("Alignment generation complete!")
    print("="*50)
    print(f"\nAlignment files written to: {args.output}")
    print(f"Alignment list: {list_file}")
    print("\nNext steps:")
    print("1. Review alignments and manually refine if needed")
    print("2. Run threading: ./scripts/run_rosettacm.sh")

if __name__ == '__main__':
    main()
