#!/usr/bin/env python3
"""
Setup RosettaCM for M. abscessus DNA Gyrase modeling using 5BS8 as template.

5BS8 structure (M. tuberculosis):
- Chains A, C: GyrA Tower + C-gate domains (residues ~15-501 of full GyrA)
- Chains B, D: GyrB TOPRIM domain (residues ~425-675 of full GyrB)
- Chains E, F, G, H: DNA (excluded from protein modeling)

This script:
1. Extracts template sequences from 5BS8
2. Extracts corresponding regions from M. abscessus GyrA/GyrB
3. Creates proper GRISHIN alignments for RosettaCM
"""

import os
from pathlib import Path

# Three-letter to one-letter amino acid code
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'PTR': 'Y', 'SEP': 'S', 'TPO': 'T',  # Modified residues
}

def extract_sequence_from_pdb(pdb_file, chain_id):
    """Extract sequence from PDB ATOM records for a specific chain."""
    sequence = []
    last_resnum = None
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                if chain != chain_id:
                    continue
                
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                
                if resnum != last_resnum:
                    if resname in THREE_TO_ONE:
                        sequence.append((resnum, THREE_TO_ONE[resname]))
                    last_resnum = resnum
    
    return sequence

def get_seqres_sequence(pdb_file, chain_id):
    """Extract sequence from SEQRES records."""
    sequence = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('SEQRES') and line[11] == chain_id:
                residues = line[19:].split()
                for res in residues:
                    if res in THREE_TO_ONE:
                        sequence.append(THREE_TO_ONE[res])
    
    return ''.join(sequence)

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Simple Needleman-Wunsch global alignment."""
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    score = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        score[i][0] = i * gap
    for j in range(n + 1):
        score[0][j] = j * gap
    
    # Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(match_score, delete, insert)
    
    # Traceback
    aligned1, aligned2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            current = score[i][j]
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            if current == diag:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        
        if i > 0 and score[i][j] == score[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))

def calculate_identity(aligned1, aligned2):
    """Calculate sequence identity."""
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')
    length = len(aligned1.replace('-', ''))
    return (matches / length * 100) if length > 0 else 0

def main():
    # Paths
    pdb_file = 'input/templates/5bs8.pdb'
    target_fasta = 'input/target.fasta'
    output_dir = 'input/alignments'
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Full M. abscessus sequences
    gyra_full = (
        "MTDTTLPPGGDDAVDRIEPVDIQQEMQRSYIDYAMSVIVGRALPEVRDGLKPVHRRVLYAMYDSGFRPDR"
        "SHAKSARSVAETMGNYHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGSPGNDPAAAMRYTEARLTPLA"
        "MEMLREIDEETVDFIPNYDGRVMEPTVLPSRFPNLLANGSGGIAVGMATNMPPHNLRELAEAVYWALDNH"
        "EADEETTLKAVCEKITGPDFPTSGLIVGTQGIHDAYTTGRGSIRMRGVAEIEEDSKGRTSLVITELPYQV"
        "NHDNFITSIAEQVRDGKIAGISNIEDQSSDRVGLRIVVVLKRDAVAKVVLNNLYKHTQLQTSFGANMLSI"
        "VDGVPRTLRLDQLIRLYVNHQLDVIIRRTRYRLRKANERAHILRGLVKALDALDEVIALIRASQTVDIAR"
        "TGLIELLDVDEIQAQAILDMQLRRLAALERQKIIDDLAKIEAEIADLEDILAKPERQRAIVKDELAEITE"
        "KYGDDRRTRIISADGDVADEDLIAREDVVVTITETGYAKRTKTDLYRSQKRGGKGVQGAGLKQDDIVKHF"
        "FVCSTHDWILFFTTKGRVYRAKAYDLPEAARTARGQHVANLLAFQPEERIAQVIQIKSYEDAPYLVLATK"
        "NGLVKKSKLTEFDSNRSGGLVAVNLRDGDELVGAVLCSAEDDLLLVSAHGQSIRFSATDEALRPMGRATS"
        "GVQGMRFNGEDDLLSLNVVREGTYLLVATSGGYSKRTAIEEYPVQGRGGKGVLTVQYDPRRGSLVGALVV"
        "DEESELYAITSGGGVIRTIAKQVRKAGRQTKGVRLMNLGEGDTLLAIAHNADEGDADPDEDAAGTTAGE"
    )
    
    gyrb_full = (
        "MAAQKKSAKSEYSADSITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNSVDEAMAGYATTVEVTML"
        "ADGGIQVKDDGRGIPVAMHASGIPTVDVVMTQLHAGGKFDSDSYAVSGGLHGVGISVVNALSTKLEVEIL"
        "RDGFEWQQVYTRSEPGTLQKGAATKKTGTTVRFWADPEIFETTTYDFETVARRLQEQAFLNKGLTIKLTD"
        "ERVSDSEVTDEVVSDTAEAPKNAEEQAAESSAPHKVKNRVFHYPDGLVDFVKHINRTKSAIHTTIVDFSG"
        "KGEGHEVEIAMQWNAGYSESVHTFANTINTHEGGTHEEGFRAALTSVVNKYAKEKKLLKEKDSNLTGDDI"
        "REGLAAVISVKVGEPQFEGQTKTKLGNTEVKSFVQKVCNEQLQHWFDSNPADAKTVVNKAVSSAQARIAA"
        "RKARELVRRKSATDIGGLPGKLADCRSTDPSKSELYVVEGDSAGGSAKSGRDSMFQAILPLRGKIINVEK"
        "ARIDRVLKNTEVQAIITALGTGIHDEFDIAKLRYHKIVLMADADVDGQHISTLLLTLLFRFMRPLVEHGH"
        "IFLAQPPLYKLKWQRTQPEFAYSDRERDGLMEAGLKAGKKINKDDGIQRYKGLGEMDAKELWETTMDPSV"
        "RVLRQVTLDDAAAADELFSILMGEDVEARRSFITRNAKDVRFLDV"
    )
    
    print("=" * 70)
    print("RosettaCM Setup for M. abscessus DNA Gyrase")
    print("Template: 5BS8 (M. tuberculosis)")
    print("=" * 70)
    
    # Check if template exists
    if not os.path.exists(pdb_file):
        print(f"\nTemplate not found at {pdb_file}")
        print("Copying 5bs8.pdb to input/templates/...")
        os.makedirs('input/templates', exist_ok=True)
        os.system('cp 5bs8.pdb input/templates/')
        pdb_file = 'input/templates/5bs8.pdb'
    
    # Extract template sequences
    print("\n--- Template Sequences (5BS8) ---")
    
    template_gyra_a = get_seqres_sequence(pdb_file, 'A')
    template_gyrb_b = get_seqres_sequence(pdb_file, 'B')
    
    print(f"Chain A (GyrA): {len(template_gyra_a)} residues")
    print(f"Chain B (GyrB): {len(template_gyrb_b)} residues")
    
    # Align and find corresponding regions
    print("\n--- Sequence Alignments ---")
    
    # GyrA alignment
    aligned_target_a, aligned_template_a = needleman_wunsch(gyra_full, template_gyra_a)
    identity_a = calculate_identity(aligned_target_a, aligned_template_a)
    print(f"\nGyrA alignment:")
    print(f"  Target (M. abscessus): {len(gyra_full)} aa")
    print(f"  Template (M. tuberculosis): {len(template_gyra_a)} aa")
    print(f"  Sequence identity: {identity_a:.1f}%")
    
    # GyrB alignment  
    aligned_target_b, aligned_template_b = needleman_wunsch(gyrb_full, template_gyrb_b)
    identity_b = calculate_identity(aligned_target_b, aligned_template_b)
    print(f"\nGyrB alignment:")
    print(f"  Target (M. abscessus): {len(gyrb_full)} aa")
    print(f"  Template (M. tuberculosis): {len(template_gyrb_b)} aa")
    print(f"  Sequence identity: {identity_b:.1f}%")
    
    # Write GRISHIN alignment files
    print("\n--- Writing Alignment Files ---")
    
    # Combined alignment file for all chains
    alignment_file = os.path.join(output_dir, 'gyrase_5bs8.grishin')
    
    with open(alignment_file, 'w') as f:
        # Chain A (GyrA copy 1)
        f.write("## GyrA_chain1 5bs8A\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_a}\n")
        f.write(f"0 {aligned_template_a}\n")
        f.write("--\n")
        
        # Chain B (GyrB copy 1)
        f.write("## GyrB_chain1 5bs8B\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_b}\n")
        f.write(f"0 {aligned_template_b}\n")
        f.write("--\n")
        
        # Chain C (GyrA copy 2)
        f.write("## GyrA_chain2 5bs8C\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_a}\n")
        f.write(f"0 {aligned_template_a}\n")
        f.write("--\n")
        
        # Chain D (GyrB copy 2)
        f.write("## GyrB_chain2 5bs8D\n")
        f.write("#\n")
        f.write("scores_from_program: 0\n")
        f.write(f"0 {aligned_target_b}\n")
        f.write(f"0 {aligned_template_b}\n")
        f.write("--\n")
    
    print(f"  Written: {alignment_file}")
    
    # Write target FASTA matching template structure
    target_file = os.path.join('input', 'target_for_5bs8.fasta')
    with open(target_file, 'w') as f:
        f.write("# M. abscessus DNA Gyrase - matching 5BS8 template structure\n")
        f.write("# Template: 5BS8 (M. tuberculosis gyrase-DNA complex)\n")
        f.write("# Complex: GyrA2-GyrB2 heterotetramer\n")
        f.write("#\n")
        f.write(">GyrA_chain1\n")
        for i in range(0, len(gyra_full), 70):
            f.write(gyra_full[i:i+70] + "\n")
        f.write(">GyrB_chain1\n")
        for i in range(0, len(gyrb_full), 70):
            f.write(gyrb_full[i:i+70] + "\n")
        f.write(">GyrA_chain2\n")
        for i in range(0, len(gyra_full), 70):
            f.write(gyra_full[i:i+70] + "\n")
        f.write(">GyrB_chain2\n")
        for i in range(0, len(gyrb_full), 70):
            f.write(gyrb_full[i:i+70] + "\n")
    
    print(f"  Written: {target_file}")
    
    # Write template sequence FASTA
    template_fasta = os.path.join(output_dir, '5bs8_template.fasta')
    with open(template_fasta, 'w') as f:
        f.write(">5bs8A_GyrA\n")
        for i in range(0, len(template_gyra_a), 70):
            f.write(template_gyra_a[i:i+70] + "\n")
        f.write(">5bs8B_GyrB\n")
        for i in range(0, len(template_gyrb_b), 70):
            f.write(template_gyrb_b[i:i+70] + "\n")
    
    print(f"  Written: {template_fasta}")
    
    print("\n" + "=" * 70)
    print("Setup Complete!")
    print("=" * 70)
    print(f"""
Next steps for RosettaCM:

1. Prepare template (remove DNA chains if modeling protein only):
   grep -v "^ATOM.* [EFGH] " 5bs8.pdb > input/templates/5bs8_protein.pdb

2. Run partial_thread for each chain:
   partial_thread.default.linuxgccrelease \\
     -in:file:fasta input/target_for_5bs8.fasta \\
     -in:file:alignment input/alignments/gyrase_5bs8.grishin \\
     -in:file:template_pdb input/templates/5bs8_protein.pdb

3. Run RosettaCM hybridize:
   rosetta_scripts.default.linuxgccrelease \\
     -parser:protocol xml/rosettacm_gyrase.xml \\
     -in:file:fasta input/target_for_5bs8.fasta \\
     -nstruct 100

Template coverage:
  - GyrA: models residues 1-503 (of 839) - {len(template_gyra_a)} template residues
  - GyrB: models residues 423-675 (of 675) - {len(template_gyrb_b)} template residues
  
Note: 5BS8 contains domain fragments, not full-length proteins.
For full-length modeling, consider using additional templates.
""")

if __name__ == '__main__':
    main()
