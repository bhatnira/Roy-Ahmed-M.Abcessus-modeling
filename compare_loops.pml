# PyMOL script to compare original vs refined loops
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
