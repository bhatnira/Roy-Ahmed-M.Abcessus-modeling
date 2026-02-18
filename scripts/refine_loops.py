#!/usr/bin/env python3
"""
Loop Refinement Script for M. abscessus Gyrase Model
Focuses on improving disjoint loop regions identified in structural alignment

Options:
1. Rosetta loop modeling (KIC protocol)
2. Energy minimization of loop regions
3. AlphaFold2 prediction for specific regions
"""

import os
import sys
import subprocess

# Disjoint regions identified from structural alignment
DISJOINT_REGIONS = {
    'GyrA': [(17, 19)],  # Only 3 residues - minor
    'GyrB': [
        (425, 430),   # 6 residues
        (446, 454),   # 9 residues  
        (472, 475),   # 4 residues
        (520, 538),   # 19 residues
        (559, 620),   # 62 residues - MAJOR LOOP
        (629, 630),   # 2 residues
        (655, 657),   # 3 residues
        (675, 675),   # 1 residue (C-terminal)
    ]
}

def create_loop_file(output_file, chain, regions):
    """Create a Rosetta loop file for KIC refinement"""
    with open(output_file, 'w') as f:
        f.write("LOOP  START  END  CUT  SKIP_RATE  EXTENDED\n")
        for start, end in regions:
            cut = (start + end) // 2  # Cut point in middle
            f.write(f"LOOP  {start}  {end}  {cut}  0.0  0\n")
    print(f"Created loop file: {output_file}")

def create_rosetta_loop_xml(output_file, loop_regions):
    """Create Rosetta XML for loop modeling"""
    
    # Build residue selector for loop regions
    loop_residues = ",".join([f"{s}-{e}" for s, e in loop_regions])
    
    xml_content = f'''<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref2015" weights="ref2015"/>
        <ScoreFunction name="ref2015_cart" weights="ref2015_cart"/>
    </SCOREFXNS>
    
    <RESIDUE_SELECTORS>
        <!-- Select the disjoint loop regions for refinement -->
        <Index name="loop_regions" resnums="{loop_residues}"/>
        <Not name="not_loops" selector="loop_regions"/>
    </RESIDUE_SELECTORS>
    
    <MOVE_MAP_FACTORIES>
        <!-- Only allow loop regions to move -->
        <MoveMapFactory name="loop_mm" bb="false" chi="false">
            <Backbone residue_selector="loop_regions" enable="true"/>
            <Chi residue_selector="loop_regions" enable="true"/>
        </MoveMapFactory>
        
        <!-- Allow everything to move for final refinement -->
        <MoveMapFactory name="full_mm" bb="true" chi="true"/>
    </MOVE_MAP_FACTORIES>
    
    <MOVERS>
        <!-- Loop modeling with KIC -->
        <LoopModeler name="loop_model" 
            scorefxn="ref2015"
            config="kic"
            loops_file="loops.txt"
            max_kic_build_attempts="1000"/>
        
        <!-- Minimize loop regions -->
        <MinMover name="min_loops" 
            scorefxn="ref2015_cart" 
            movemap_factory="loop_mm"
            type="lbfgs_armijo_nonmonotone"
            tolerance="0.01"
            cartesian="true"/>
        
        <!-- Final relaxation of loop regions -->
        <FastRelax name="relax_loops" 
            scorefxn="ref2015" 
            repeats="3"
            movemap_factory="loop_mm"/>
        
        <!-- Optional: Full relax at the end -->
        <FastRelax name="full_relax" 
            scorefxn="ref2015" 
            repeats="1"
            movemap_factory="full_mm"/>
    </MOVERS>
    
    <PROTOCOLS>
        <!-- Option 1: Just minimize loops (fast) -->
        <!-- <Add mover="min_loops"/> -->
        
        <!-- Option 2: Relax loops (moderate) -->
        <Add mover="relax_loops"/>
        
        <!-- Option 3: Full protocol (slow but best) -->
        <!-- <Add mover="loop_model"/> -->
        <!-- <Add mover="relax_loops"/> -->
        <!-- <Add mover="full_relax"/> -->
    </PROTOCOLS>
    
    <OUTPUT scorefxn="ref2015"/>
</ROSETTASCRIPTS>
'''
    
    with open(output_file, 'w') as f:
        f.write(xml_content)
    print(f"Created Rosetta XML: {output_file}")

def create_pymol_visualization(output_file):
    """Create PyMOL script to visualize before/after refinement"""
    
    pml_content = '''# PyMOL script to compare original vs refined loops
# Usage: pymol compare_loops.pml

# Load structures
load output/tetramer/mabs_gyrase_tetramer_protein_only.pdb, original
# load output/refined/mabs_gyrase_refined.pdb, refined  # Uncomment after refinement

# Setup visualization
hide everything
show cartoon

# Color scheme
color cyan, original
# color green, refined  # Uncomment after refinement

# Highlight loop regions (GyrB disjoint loops)
select loops, chain B and resi 425-430+446-454+472-475+520-538+559-620+629-630+655-657+675
color red, loops and original
# color lime, loops and refined  # Uncomment after refinement

# Show loops as sticks for detail
show sticks, loops

# Zoom to major loop
zoom chain B and resi 559-620

# Set nice rendering
bg_color white
set cartoon_fancy_helices, 1
set ray_shadows, 0

print("Loop regions highlighted in red (original)")
print("Major disjoint loop: residues 559-620 in GyrB")
'''
    
    with open(output_file, 'w') as f:
        f.write(pml_content)
    print(f"Created PyMOL script: {output_file}")

def print_refinement_options():
    """Print available refinement strategies"""
    
    print("""
╔═══════════════════════════════════════════════════════════════════════════════╗
║                     LOOP REFINEMENT OPTIONS                                    ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║                                                                               ║
║  OPTION 1: Rosetta FastRelax (Recommended - Quick)                            ║
║  ─────────────────────────────────────────────────                            ║
║  • Relaxes loop conformations while keeping core fixed                        ║
║  • Runtime: ~30 min - 2 hours                                                 ║
║  • Command:                                                                   ║
║    rosetta_scripts -s input.pdb -parser:protocol loop_relax.xml               ║
║                                                                               ║
║  OPTION 2: Rosetta KIC Loop Modeling (Best Quality - Slow)                    ║
║  ─────────────────────────────────────────────────────────                    ║
║  • Rebuilds loops from scratch using kinematic closure                        ║
║  • Runtime: Several hours to days                                             ║
║  • Best for the major 559-620 loop                                            ║
║                                                                               ║
║  OPTION 3: AlphaFold2/ColabFold (Best for Long Loops)                         ║
║  ─────────────────────────────────────────────────────                        ║
║  • Use AlphaFold to predict just the GyrB subunit                             ║
║  • Then align and graft the better loop onto your model                       ║
║  • Especially good for the 62-residue loop (559-620)                          ║
║                                                                               ║
║  OPTION 4: Molecular Dynamics (Most Rigorous)                                 ║
║  ────────────────────────────────────────────                                 ║
║  • Run short MD simulation to relax loop conformations                        ║
║  • Use GROMACS or OpenMM                                                      ║
║  • Runtime: Hours to days                                                     ║
║                                                                               ║
║  OPTION 5: Accept Current Model ✓                                             ║
║  ────────────────────────────────                                             ║
║  • Core RMSD is excellent (0.75-1.22 Å)                                       ║
║  • Disjoint loops are often flexible/disordered in nature                     ║
║  • If loops are not in drug binding site, refinement may not be needed        ║
║                                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝
""")

def main():
    base_dir = "/home/nbhatta1/Desktop/Roy-Ahmed-M.Abcessus-modeling"
    os.makedirs(f"{base_dir}/output/refined", exist_ok=True)
    
    print_refinement_options()
    
    # Create loop file for GyrB (the problematic chain)
    gyrb_regions = DISJOINT_REGIONS['GyrB']
    create_loop_file(f"{base_dir}/output/refined/gyrB_loops.txt", 'B', gyrb_regions)
    
    # Create Rosetta XML
    create_rosetta_loop_xml(f"{base_dir}/xml/loop_refinement.xml", gyrb_regions)
    
    # Create visualization script
    create_pymol_visualization(f"{base_dir}/compare_loops.pml")
    
    print("\n" + "="*70)
    print("FILES CREATED:")
    print("="*70)
    print(f"  1. {base_dir}/output/refined/gyrB_loops.txt")
    print(f"  2. {base_dir}/xml/loop_refinement.xml")
    print(f"  3. {base_dir}/compare_loops.pml")
    print("\n" + "="*70)
    print("TO RUN LOOP REFINEMENT:")
    print("="*70)
    print("""
  # Using Rosetta (if installed):
  rosetta_scripts.default.linuxgccrelease \\
    -s output/tetramer/mabs_gyrase_tetramer_protein_only.pdb \\
    -parser:protocol xml/loop_refinement.xml \\
    -out:prefix refined_ \\
    -nstruct 5
  
  # Or use the Docker version:
  ./scripts/run_docker_rosettacm.sh xml/loop_refinement.xml \\
    output/tetramer/mabs_gyrase_tetramer_protein_only.pdb
""")

if __name__ == "__main__":
    main()
