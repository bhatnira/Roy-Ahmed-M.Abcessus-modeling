#!/usr/bin/env python3
"""
Check conservation of fluoroquinolone binding site residues
between M. tuberculosis and M. abscessus gyrase

Key binding site residues (E. coli numbering):
- GyrA Ser83 (QRDR) - hydrogen bonding with fluoroquinolone
- GyrA Asp87 (QRDR) - Mg2+ coordination
- GyrA Ala90 - hydrophobic contact
- GyrA Tyr129 - catalytic tyrosine (DNA cleavage)
"""

def parse_grishin_alignment(filepath):
    """Parse Grishin alignment file"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    sequences = {}
    for i, line in enumerate(lines):
        if line.startswith('0 '):
            seq = line[2:].strip()
            if i == 3:  # Target sequence (M. abscessus)
                sequences['mabs'] = seq
            elif i == 4:  # Template sequence (M. tuberculosis)
                sequences['mtb'] = seq
    
    return sequences

def align_residue_mapping(target_seq, template_seq):
    """
    Create mapping between target and template residue numbers
    Returns: dict mapping target position -> template position (1-indexed)
    """
    target_pos = 0  # 1-indexed position in target
    template_pos = 0  # 1-indexed position in template
    
    mapping = {}
    
    for i in range(len(target_seq)):
        target_char = target_seq[i] if i < len(target_seq) else '-'
        template_char = template_seq[i] if i < len(template_seq) else '-'
        
        if target_char != '-':
            target_pos += 1
        if template_char != '-':
            template_pos += 1
        
        if target_char != '-' and template_char != '-':
            mapping[target_pos] = {
                'target_res': target_char,
                'template_res': template_char,
                'template_pos': template_pos
            }
    
    return mapping

def find_conserved_binding_residues():
    """Analyze conservation of key fluoroquinolone binding residues"""
    
    print("=" * 70)
    print("FLUOROQUINOLONE BINDING SITE CONSERVATION ANALYSIS")
    print("M. tuberculosis vs M. abscessus Gyrase")
    print("=" * 70)
    
    # Parse alignments
    gyrA_path = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling/input/alignments_final/alignment_gyrA_proper.grishin"
    gyrB_path = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling/input/alignments_final/alignment_gyrB_proper.grishin"
    
    # Read raw alignments
    with open(gyrA_path) as f:
        gyrA_lines = f.readlines()
    with open(gyrB_path) as f:
        gyrB_lines = f.readlines()
    
    # Extract aligned sequences
    mabs_gyrA = gyrA_lines[3][2:].strip()  # M. abscessus
    mtb_gyrA = gyrA_lines[4][2:].strip()   # M. tuberculosis (5BS8)
    
    mabs_gyrB = gyrB_lines[3][2:].strip()
    mtb_gyrB = gyrB_lines[4][2:].strip()
    
    print("\n" + "-" * 70)
    print("GyrA Analysis")
    print("-" * 70)
    
    # Key binding residues in M. tuberculosis GyrA (5BS8 numbering)
    # From literature: fluoroquinolone resistance mutations
    mtb_key_residues = {
        'GyrA': {
            # QRDR region (Quinolone Resistance Determining Region)
            # M. tuberculosis numbering based on 5BS8 structure
            90: ('A', 'QRDR - hydrophobic contact'),  # Ala90
            94: ('D', 'QRDR - Mg2+ coordination (equivalent to E. coli D87)'),  # Asp94
            91: ('S', 'QRDR - H-bond with FQ (equivalent to E. coli S83)'),  # Ser91
            88: ('G', 'QRDR - Gly88'),
        }
    }
    
    # Create position mapping from alignment
    print("\nAligned sequences (first 150 chars):")
    print(f"M.abscessus: {mabs_gyrA[:150]}...")
    print(f"M.tuberc:    {mtb_gyrA[:150]}...")
    
    # Find positions of key residues in alignment
    print("\n" + "=" * 70)
    print("KEY BINDING SITE RESIDUES - GyrA QRDR Region")
    print("=" * 70)
    
    # Walk through alignment to find equivalent positions
    mabs_pos = 0
    mtb_pos = 0
    
    residue_comparison = []
    
    for i in range(min(len(mabs_gyrA), len(mtb_gyrA))):
        if mabs_gyrA[i] != '-':
            mabs_pos += 1
        if mtb_gyrA[i] != '-':
            mtb_pos += 1
        
        # Check if this is a key MTB position
        if mtb_pos in mtb_key_residues['GyrA'] and mtb_gyrA[i] != '-':
            expected_res, description = mtb_key_residues['GyrA'][mtb_pos]
            mabs_res = mabs_gyrA[i] if mabs_gyrA[i] != '-' else 'GAP'
            mtb_res = mtb_gyrA[i]
            
            conserved = "CONSERVED" if mabs_res == mtb_res else "DIFFERENT"
            
            residue_comparison.append({
                'mtb_pos': mtb_pos,
                'mabs_pos': mabs_pos if mabs_gyrA[i] != '-' else '-',
                'mtb_res': mtb_res,
                'mabs_res': mabs_res,
                'conserved': conserved,
                'description': description
            })
    
    print(f"\n{'MTB Pos':<10} {'MTB Res':<10} {'MAbs Pos':<10} {'MAbs Res':<10} {'Status':<12} Description")
    print("-" * 90)
    for r in residue_comparison:
        print(f"{r['mtb_pos']:<10} {r['mtb_res']:<10} {r['mabs_pos']:<10} {r['mabs_res']:<10} {r['conserved']:<12} {r['description']}")
    
    # Analyze full QRDR region (typically positions 70-110)
    print("\n" + "=" * 70)
    print("FULL QRDR REGION ALIGNMENT (positions ~70-110)")
    print("=" * 70)
    
    mabs_pos = 0
    mtb_pos = 0
    qrdr_start_mtb = 70
    qrdr_end_mtb = 110
    
    qrdr_mabs = []
    qrdr_mtb = []
    
    for i in range(min(len(mabs_gyrA), len(mtb_gyrA))):
        if mabs_gyrA[i] != '-':
            mabs_pos += 1
        if mtb_gyrA[i] != '-':
            mtb_pos += 1
        
        if qrdr_start_mtb <= mtb_pos <= qrdr_end_mtb:
            qrdr_mabs.append(mabs_gyrA[i])
            qrdr_mtb.append(mtb_gyrA[i])
    
    qrdr_mabs_str = ''.join(qrdr_mabs)
    qrdr_mtb_str = ''.join(qrdr_mtb)
    
    print(f"\nM. abscessus QRDR: {qrdr_mabs_str}")
    print(f"M. tuberculosis:   {qrdr_mtb_str}")
    
    # Calculate identity
    matches = sum(1 for a, b in zip(qrdr_mabs_str, qrdr_mtb_str) if a == b and a != '-')
    total = sum(1 for a, b in zip(qrdr_mabs_str, qrdr_mtb_str) if a != '-' and b != '-')
    identity = matches / total * 100 if total > 0 else 0
    
    print(f"\nQRDR Identity: {matches}/{total} = {identity:.1f}%")
    
    # Show differences
    print("\n" + "=" * 70)
    print("DIFFERENCES IN QRDR REGION")
    print("=" * 70)
    
    mabs_pos = 0
    mtb_pos = 0
    differences = []
    
    for i in range(min(len(mabs_gyrA), len(mtb_gyrA))):
        if mabs_gyrA[i] != '-':
            mabs_pos += 1
        if mtb_gyrA[i] != '-':
            mtb_pos += 1
        
        if qrdr_start_mtb <= mtb_pos <= qrdr_end_mtb:
            if mabs_gyrA[i] != mtb_gyrA[i] and mabs_gyrA[i] != '-' and mtb_gyrA[i] != '-':
                differences.append({
                    'mtb_pos': mtb_pos,
                    'mabs_pos': mabs_pos,
                    'mtb_res': mtb_gyrA[i],
                    'mabs_res': mabs_gyrA[i]
                })
    
    if differences:
        print(f"\n{'MTB Pos':<10} {'MTB Res':<10} {'MAbs Pos':<10} {'MAbs Res':<10}")
        print("-" * 50)
        for d in differences:
            print(f"{d['mtb_pos']:<10} {d['mtb_res']:<10} {d['mabs_pos']:<10} {d['mabs_res']:<10}")
    else:
        print("\nNo differences in QRDR region - binding site is FULLY CONSERVED!")
    
    # GyrB Analysis
    print("\n" + "=" * 70)
    print("GyrB Analysis - Additional binding site residues")
    print("=" * 70)
    
    print("\nAligned sequences (relevant region):")
    # Find where alignment starts (skip gaps at beginning)
    mtb_start = 0
    for i, c in enumerate(mtb_gyrB):
        if c != '-':
            mtb_start = i
            break
    
    print(f"M.abscessus: ...{mabs_gyrB[mtb_start:mtb_start+100]}...")
    print(f"M.tuberc:    ...{mtb_gyrB[mtb_start:mtb_start+100]}...")
    
    # Calculate overall identity for aligned region
    matches = sum(1 for a, b in zip(mabs_gyrB, mtb_gyrB) if a == b and a != '-')
    total = sum(1 for a, b in zip(mabs_gyrB, mtb_gyrB) if a != '-' and b != '-')
    identity = matches / total * 100 if total > 0 else 0
    
    print(f"\nGyrB aligned region identity: {matches}/{total} = {identity:.1f}%")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
Key findings for fluoroquinolone binding site conservation:

1. QRDR (Quinolone Resistance Determining Region) in GyrA:
   - This region contains the critical residues for FQ binding
   - High conservation suggests similar drug sensitivity
   
2. Important residues for FQ binding:
   - Ser91 (MTB) / equivalent in M.abs - H-bonding with FQ C3/C4 keto group
   - Asp94 (MTB) / equivalent in M.abs - Mg2+ coordination
   - Ala90 (MTB) / equivalent in M.abs - hydrophobic interactions
   
3. Clinical relevance:
   - Mutations in QRDR are associated with fluoroquinolone resistance
   - S91P (MTB) and D94G/N/Y (MTB) are common resistance mutations
   - Conservation suggests M. abscessus may have similar resistance mechanisms
""")
    
    return residue_comparison

if __name__ == "__main__":
    find_conserved_binding_residues()
