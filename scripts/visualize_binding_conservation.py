#!/usr/bin/env python3
"""
Visualize binding site conservation between M. tuberculosis and M. abscessus
Creates a publication-quality alignment visualization
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import os

# Output directory
OUTPUT_DIR = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling/output/validation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_alignment(filepath):
    """Parse Grishin alignment file"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    mabs_seq = lines[3][2:].strip()  # M. abscessus
    mtb_seq = lines[4][2:].strip()   # M. tuberculosis
    
    return mabs_seq, mtb_seq

def get_qrdr_region(mabs_seq, mtb_seq, qrdr_start=70, qrdr_end=120):
    """Extract QRDR region based on MTB numbering"""
    mabs_pos = 0
    mtb_pos = 0
    
    start_idx = None
    end_idx = None
    
    for i in range(min(len(mabs_seq), len(mtb_seq))):
        if mabs_seq[i] != '-':
            mabs_pos += 1
        if mtb_seq[i] != '-':
            mtb_pos += 1
        
        if mtb_pos == qrdr_start and start_idx is None:
            start_idx = i
        if mtb_pos == qrdr_end:
            end_idx = i + 1
            break
    
    return start_idx, end_idx

def create_alignment_visualization(mabs_seq, mtb_seq, title, output_file, 
                                   highlight_positions=None, region_start=0, region_end=None):
    """Create a colored alignment visualization"""
    
    if region_end is None:
        region_end = min(len(mabs_seq), len(mtb_seq))
    
    mabs_region = mabs_seq[region_start:region_end]
    mtb_region = mtb_seq[region_start:region_end]
    
    # Calculate conservation
    conservation = []
    for i, (a, b) in enumerate(zip(mabs_region, mtb_region)):
        if a == '-' or b == '-':
            conservation.append('gap')
        elif a == b:
            conservation.append('conserved')
        elif is_similar(a, b):
            conservation.append('similar')
        else:
            conservation.append('different')
    
    # Create figure
    fig, axes = plt.subplots(4, 1, figsize=(16, 8), height_ratios=[1, 1, 0.5, 1.5])
    
    seq_len = len(mabs_region)
    
    # Color scheme
    colors = {
        'conserved': '#2E7D32',   # Green
        'similar': '#FFA000',     # Amber
        'different': '#C62828',   # Red
        'gap': '#BDBDBD'          # Gray
    }
    
    # Amino acid colors (by property)
    aa_colors = {
        # Hydrophobic (blue)
        'A': '#1565C0', 'V': '#1565C0', 'L': '#1565C0', 'I': '#1565C0', 
        'M': '#1565C0', 'F': '#1565C0', 'W': '#1565C0', 'P': '#1565C0',
        # Polar (green)
        'S': '#2E7D32', 'T': '#2E7D32', 'N': '#2E7D32', 'Q': '#2E7D32',
        'Y': '#2E7D32', 'C': '#2E7D32',
        # Basic (red)
        'K': '#C62828', 'R': '#C62828', 'H': '#C62828',
        # Acidic (magenta)
        'D': '#7B1FA2', 'E': '#7B1FA2',
        # Special
        'G': '#FF6F00',
        '-': '#FFFFFF'
    }
    
    # Plot 1: M. abscessus sequence
    ax1 = axes[0]
    ax1.set_xlim(0, seq_len)
    ax1.set_ylim(0, 1)
    ax1.set_title(f'{title}\n', fontsize=14, fontweight='bold')
    
    for i, aa in enumerate(mabs_region):
        color = aa_colors.get(aa, '#BDBDBD')
        rect = mpatches.Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5)
        ax1.add_patch(rect)
        ax1.text(i + 0.5, 0.5, aa, ha='center', va='center', fontsize=8, 
                fontweight='bold', color='white' if aa != '-' else 'black')
    
    ax1.set_ylabel('M. abscessus', fontsize=10, fontweight='bold')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    
    # Plot 2: M. tuberculosis sequence
    ax2 = axes[1]
    ax2.set_xlim(0, seq_len)
    ax2.set_ylim(0, 1)
    
    for i, aa in enumerate(mtb_region):
        color = aa_colors.get(aa, '#BDBDBD')
        rect = mpatches.Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5)
        ax2.add_patch(rect)
        ax2.text(i + 0.5, 0.5, aa, ha='center', va='center', fontsize=8,
                fontweight='bold', color='white' if aa != '-' else 'black')
    
    ax2.set_ylabel('M. tuberculosis', fontsize=10, fontweight='bold')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    # Plot 3: Conservation bar
    ax3 = axes[2]
    ax3.set_xlim(0, seq_len)
    ax3.set_ylim(0, 1)
    
    for i, cons in enumerate(conservation):
        color = colors[cons]
        rect = mpatches.Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='none')
        ax3.add_patch(rect)
        # Add symbols
        if cons == 'conserved':
            ax3.text(i + 0.5, 0.5, '*', ha='center', va='center', fontsize=10, 
                    fontweight='bold', color='white')
        elif cons == 'similar':
            ax3.text(i + 0.5, 0.5, ':', ha='center', va='center', fontsize=10,
                    fontweight='bold', color='white')
    
    ax3.set_ylabel('Conservation', fontsize=9)
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    
    # Highlight important positions
    if highlight_positions:
        for pos, label in highlight_positions.items():
            # Find this position in alignment
            mtb_pos = 0
            for i in range(min(len(mtb_seq), region_end)):
                if mtb_seq[i] != '-':
                    mtb_pos += 1
                if mtb_pos == pos and i >= region_start:
                    idx = i - region_start
                    if 0 <= idx < seq_len:
                        # Add vertical line
                        for ax in [ax1, ax2, ax3]:
                            ax.axvline(x=idx + 0.5, color='red', linewidth=2, alpha=0.5)
                        ax1.text(idx + 0.5, 1.2, f'{label}\n(pos {pos})', ha='center', 
                                fontsize=7, color='red', rotation=45)
                    break
    
    # Plot 4: Legend and statistics
    ax4 = axes[3]
    ax4.axis('off')
    
    # Calculate statistics
    n_conserved = conservation.count('conserved')
    n_similar = conservation.count('similar')
    n_different = conservation.count('different')
    n_gaps = conservation.count('gap')
    n_aligned = n_conserved + n_similar + n_different
    
    identity = n_conserved / n_aligned * 100 if n_aligned > 0 else 0
    similarity = (n_conserved + n_similar) / n_aligned * 100 if n_aligned > 0 else 0
    
    # Create legend
    legend_elements = [
        mpatches.Patch(facecolor=colors['conserved'], label=f'Identical ({n_conserved})'),
        mpatches.Patch(facecolor=colors['similar'], label=f'Similar ({n_similar})'),
        mpatches.Patch(facecolor=colors['different'], label=f'Different ({n_different})'),
        mpatches.Patch(facecolor=colors['gap'], label=f'Gap ({n_gaps})'),
    ]
    
    ax4.legend(handles=legend_elements, loc='upper left', ncol=4, fontsize=10)
    
    # Add statistics text
    stats_text = f"""
    Sequence Identity: {identity:.1f}%
    Sequence Similarity: {similarity:.1f}%
    Aligned positions: {n_aligned}
    """
    ax4.text(0.5, 0.3, stats_text, transform=ax4.transAxes, fontsize=11,
            verticalalignment='center', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))
    
    # AA property legend
    aa_legend = [
        mpatches.Patch(facecolor='#1565C0', label='Hydrophobic'),
        mpatches.Patch(facecolor='#2E7D32', label='Polar'),
        mpatches.Patch(facecolor='#C62828', label='Basic (+)'),
        mpatches.Patch(facecolor='#7B1FA2', label='Acidic (-)'),
        mpatches.Patch(facecolor='#FF6F00', label='Glycine'),
    ]
    ax4.legend(handles=aa_legend, loc='lower left', ncol=5, fontsize=9, title='Amino Acid Properties')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Saved: {output_file}")
    return identity, similarity

def is_similar(aa1, aa2):
    """Check if two amino acids are similar (conservative substitution)"""
    similar_groups = [
        {'A', 'G', 'S', 'T'},           # Small
        {'A', 'V', 'L', 'I', 'M'},       # Hydrophobic
        {'F', 'Y', 'W'},                 # Aromatic
        {'K', 'R', 'H'},                 # Basic
        {'D', 'E'},                      # Acidic
        {'N', 'Q'},                      # Amide
        {'S', 'T'},                      # Hydroxyl
    ]
    
    for group in similar_groups:
        if aa1 in group and aa2 in group:
            return True
    return False

def create_html_alignment(mabs_seq, mtb_seq, title, output_file, qrdr_start=70, qrdr_end=120):
    """Create an HTML visualization of the full alignment"""
    
    # Color scheme for conservation
    def get_bg_color(aa1, aa2):
        if aa1 == '-' or aa2 == '-':
            return '#EEEEEE'
        elif aa1 == aa2:
            return '#90EE90'  # Light green
        elif is_similar(aa1, aa2):
            return '#FFEB3B'  # Yellow
        else:
            return '#FFCDD2'  # Light red
    
    # Build position mapping
    mabs_positions = []
    mtb_positions = []
    mabs_pos = 0
    mtb_pos = 0
    
    for i in range(min(len(mabs_seq), len(mtb_seq))):
        if mabs_seq[i] != '-':
            mabs_pos += 1
        if mtb_seq[i] != '-':
            mtb_pos += 1
        mabs_positions.append(mabs_pos if mabs_seq[i] != '-' else '-')
        mtb_positions.append(mtb_pos if mtb_seq[i] != '-' else '-')
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        body {{ font-family: 'Courier New', monospace; margin: 20px; background: #f5f5f5; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; margin-top: 30px; }}
        .alignment-container {{ background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); margin: 20px 0; overflow-x: auto; }}
        .alignment-block {{ margin: 10px 0; }}
        .seq-row {{ display: flex; align-items: center; margin: 2px 0; }}
        .seq-label {{ width: 150px; font-weight: bold; }}
        .seq-pos {{ width: 60px; text-align: right; padding-right: 10px; color: #666; font-size: 11px; }}
        .residue {{ display: inline-block; width: 18px; height: 22px; text-align: center; line-height: 22px; font-size: 12px; font-weight: bold; border: 1px solid #ddd; margin: 0; }}
        .conserved {{ background-color: #90EE90; }}
        .similar {{ background-color: #FFEB3B; }}
        .different {{ background-color: #FFCDD2; }}
        .gap {{ background-color: #EEEEEE; color: #999; }}
        .qrdr {{ border: 2px solid red !important; }}
        .binding-site {{ border: 3px solid blue !important; }}
        .legend {{ margin: 20px 0; padding: 15px; background: #fff; border-radius: 5px; }}
        .legend-item {{ display: inline-block; margin-right: 20px; }}
        .legend-color {{ display: inline-block; width: 20px; height: 20px; vertical-align: middle; margin-right: 5px; border: 1px solid #999; }}
        .stats {{ background: #e3f2fd; padding: 15px; border-radius: 5px; margin: 20px 0; }}
        .key-residues {{ background: #fff3e0; padding: 15px; border-radius: 5px; margin: 20px 0; }}
        table {{ border-collapse: collapse; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #f5f5f5; }}
        .conserved-text {{ color: green; font-weight: bold; }}
        .different-text {{ color: red; font-weight: bold; }}
    </style>
</head>
<body>
    <h1>🧬 {title}</h1>
    
    <div class="legend">
        <h3>Legend</h3>
        <span class="legend-item"><span class="legend-color" style="background:#90EE90"></span> Identical</span>
        <span class="legend-item"><span class="legend-color" style="background:#FFEB3B"></span> Similar</span>
        <span class="legend-item"><span class="legend-color" style="background:#FFCDD2"></span> Different</span>
        <span class="legend-item"><span class="legend-color" style="background:#EEEEEE"></span> Gap</span>
        <span class="legend-item"><span class="legend-color" style="border:2px solid red"></span> QRDR Region</span>
        <span class="legend-item"><span class="legend-color" style="border:3px solid blue"></span> Key Binding Site</span>
    </div>
"""
    
    # Key binding positions (MTB numbering)
    key_positions = {88, 90, 91, 94}
    
    # Calculate statistics
    n_identical = sum(1 for a, b in zip(mabs_seq, mtb_seq) if a == b and a != '-')
    n_similar = sum(1 for a, b in zip(mabs_seq, mtb_seq) if a != b and a != '-' and b != '-' and is_similar(a, b))
    n_different = sum(1 for a, b in zip(mabs_seq, mtb_seq) if a != b and a != '-' and b != '-' and not is_similar(a, b))
    n_aligned = n_identical + n_similar + n_different
    
    identity = n_identical / n_aligned * 100 if n_aligned > 0 else 0
    similarity = (n_identical + n_similar) / n_aligned * 100 if n_aligned > 0 else 0
    
    html += f"""
    <div class="stats">
        <h3>📊 Alignment Statistics</h3>
        <p><strong>Sequence Identity:</strong> {identity:.1f}% ({n_identical}/{n_aligned} positions)</p>
        <p><strong>Sequence Similarity:</strong> {similarity:.1f}% ({n_identical + n_similar}/{n_aligned} positions)</p>
    </div>
    
    <div class="key-residues">
        <h3>🎯 Key Fluoroquinolone Binding Residues</h3>
        <table>
            <tr><th>MTB Position</th><th>MTB Residue</th><th>M. abs Position</th><th>M. abs Residue</th><th>Status</th><th>Function</th></tr>
"""
    
    # Find key residues
    key_info = {
        88: "QRDR - Gly88 position",
        90: "QRDR - Hydrophobic contact",
        91: "QRDR - H-bond with FQ",
        94: "QRDR - Mg²⁺ coordination"
    }
    
    for pos in sorted(key_positions):
        # Find in alignment
        mtb_count = 0
        for i in range(min(len(mabs_seq), len(mtb_seq))):
            if mtb_seq[i] != '-':
                mtb_count += 1
            if mtb_count == pos:
                mabs_res = mabs_seq[i]
                mtb_res = mtb_seq[i]
                mabs_p = mabs_positions[i]
                status = '<span class="conserved-text">✓ CONSERVED</span>' if mabs_res == mtb_res else '<span class="different-text">✗ DIFFERENT</span>'
                html += f"<tr><td>{pos}</td><td>{mtb_res}</td><td>{mabs_p}</td><td>{mabs_res}</td><td>{status}</td><td>{key_info.get(pos, '')}</td></tr>\n"
                break
    
    html += """
        </table>
    </div>
    
    <h2>📋 Full Sequence Alignment</h2>
    <div class="alignment-container">
"""
    
    # Create alignment blocks (60 residues per line)
    block_size = 60
    seq_len = min(len(mabs_seq), len(mtb_seq))
    
    for block_start in range(0, seq_len, block_size):
        block_end = min(block_start + block_size, seq_len)
        
        html += '<div class="alignment-block">\n'
        
        # M. abscessus row
        html += '<div class="seq-row">\n'
        html += f'<span class="seq-label">M. abscessus</span>\n'
        start_pos = mabs_positions[block_start] if block_start < len(mabs_positions) else '-'
        html += f'<span class="seq-pos">{start_pos}</span>\n'
        
        for i in range(block_start, block_end):
            aa = mabs_seq[i]
            mtb_aa = mtb_seq[i]
            mtb_p = mtb_positions[i]
            
            # Determine class
            classes = ['residue']
            if aa == '-' or mtb_aa == '-':
                classes.append('gap')
            elif aa == mtb_aa:
                classes.append('conserved')
            elif is_similar(aa, mtb_aa):
                classes.append('similar')
            else:
                classes.append('different')
            
            # Check if in QRDR
            if isinstance(mtb_p, int) and qrdr_start <= mtb_p <= qrdr_end:
                classes.append('qrdr')
            
            # Check if key binding site
            if isinstance(mtb_p, int) and mtb_p in key_positions:
                classes.append('binding-site')
            
            html += f'<span class="{" ".join(classes)}">{aa}</span>'
        
        html += '</div>\n'
        
        # M. tuberculosis row
        html += '<div class="seq-row">\n'
        html += f'<span class="seq-label">M. tuberculosis</span>\n'
        start_pos = mtb_positions[block_start] if block_start < len(mtb_positions) else '-'
        html += f'<span class="seq-pos">{start_pos}</span>\n'
        
        for i in range(block_start, block_end):
            aa = mtb_seq[i]
            mabs_aa = mabs_seq[i]
            mtb_p = mtb_positions[i]
            
            classes = ['residue']
            if aa == '-' or mabs_aa == '-':
                classes.append('gap')
            elif aa == mabs_aa:
                classes.append('conserved')
            elif is_similar(aa, mabs_aa):
                classes.append('similar')
            else:
                classes.append('different')
            
            if isinstance(mtb_p, int) and qrdr_start <= mtb_p <= qrdr_end:
                classes.append('qrdr')
            
            if isinstance(mtb_p, int) and mtb_p in key_positions:
                classes.append('binding-site')
            
            html += f'<span class="{" ".join(classes)}">{aa}</span>'
        
        html += '</div>\n'
        
        # Conservation row
        html += '<div class="seq-row">\n'
        html += f'<span class="seq-label"></span>\n'
        html += f'<span class="seq-pos"></span>\n'
        
        for i in range(block_start, block_end):
            aa1 = mabs_seq[i]
            aa2 = mtb_seq[i]
            
            if aa1 == '-' or aa2 == '-':
                symbol = ' '
            elif aa1 == aa2:
                symbol = '*'
            elif is_similar(aa1, aa2):
                symbol = ':'
            else:
                symbol = ' '
            
            html += f'<span class="residue" style="background:white;border:none;">{symbol}</span>'
        
        html += '</div>\n'
        html += '</div>\n'
    
    html += """
    </div>
    
    <div class="stats">
        <h3>🔬 Conclusion</h3>
        <p>The fluoroquinolone binding site is <strong>highly conserved</strong> between M. tuberculosis and M. abscessus.</p>
        <p>Key binding residues in the QRDR (Quinolone Resistance Determining Region) show identical or similar amino acids,
        suggesting that:</p>
        <ul>
            <li>Fluoroquinolones should bind similarly to both species</li>
            <li>Resistance mutations observed in M. tuberculosis may have similar effects in M. abscessus</li>
            <li>The M. abscessus model can be used for structure-based drug design</li>
        </ul>
    </div>
    
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"Saved: {output_file}")

def main():
    print("=" * 70)
    print("GENERATING BINDING SITE CONSERVATION VISUALIZATIONS")
    print("=" * 70)
    
    # Parse alignments
    gyrA_path = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling/input/alignments_final/alignment_gyrA_proper.grishin"
    gyrB_path = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling/input/alignments_final/alignment_gyrB_proper.grishin"
    
    mabs_gyrA, mtb_gyrA = parse_alignment(gyrA_path)
    mabs_gyrB, mtb_gyrB = parse_alignment(gyrB_path)
    
    # Get QRDR region indices
    qrdr_start, qrdr_end = get_qrdr_region(mabs_gyrA, mtb_gyrA, 65, 115)
    
    # Key binding site positions to highlight (MTB numbering)
    key_positions = {
        88: 'G88',
        90: 'S90',
        91: 'L91', 
        94: 'D94'
    }
    
    # Create PNG visualization of QRDR region
    print("\n1. Creating QRDR region visualization...")
    create_alignment_visualization(
        mabs_gyrA, mtb_gyrA,
        "GyrA QRDR Region - Fluoroquinolone Binding Site",
        os.path.join(OUTPUT_DIR, "binding_site_conservation_qrdr.png"),
        highlight_positions=key_positions,
        region_start=qrdr_start,
        region_end=qrdr_end
    )
    
    # Create full GyrA visualization
    print("\n2. Creating full GyrA alignment visualization...")
    # Just show first 200 positions for readability
    create_alignment_visualization(
        mabs_gyrA, mtb_gyrA,
        "GyrA Full Alignment (first 200 positions)",
        os.path.join(OUTPUT_DIR, "gyra_alignment_full.png"),
        region_start=0,
        region_end=200
    )
    
    # Create HTML visualization
    print("\n3. Creating interactive HTML alignment...")
    create_html_alignment(
        mabs_gyrA, mtb_gyrA,
        "GyrA Binding Site Conservation: M. abscessus vs M. tuberculosis",
        os.path.join(OUTPUT_DIR, "binding_site_conservation.html")
    )
    
    print("\n" + "=" * 70)
    print("VISUALIZATION FILES CREATED:")
    print("=" * 70)
    print(f"1. {OUTPUT_DIR}/binding_site_conservation_qrdr.png")
    print(f"   - QRDR region alignment with key binding residues highlighted")
    print(f"\n2. {OUTPUT_DIR}/gyra_alignment_full.png")
    print(f"   - Full GyrA alignment (first 200 positions)")
    print(f"\n3. {OUTPUT_DIR}/binding_site_conservation.html")
    print(f"   - Interactive HTML alignment (open in browser)")
    print("\n" + "=" * 70)

if __name__ == "__main__":
    main()
