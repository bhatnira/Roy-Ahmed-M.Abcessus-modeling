# PyMOL script to visualize fluoroquinolone binding site
# Run: pymol visualize_binding_site.pml

# Load structures
load ../output/tetramer/mabs_gyrase_tetramer_protein_only.pdb, protein
load reference_gyrase_fq.pdb, reference

# Load ligands
load ligands/moxifloxacin.sdf, moxifloxacin
load ligands/ciprofloxacin.sdf, ciprofloxacin

# Hide everything initially
hide everything

# Show protein as cartoon
show cartoon, protein
color cyan, protein and chain A+C
color lightblue, protein and chain B+D

# Superimpose reference to our model (GyrA chains)
super reference and chain A, protein and chain A

# Show reference ligand (MFX from 5CDQ)
show sticks, reference and resn MFX
show spheres, reference and resn MFX and name MG
color yellow, reference and resn MFX
set sphere_scale, 0.3

# QRDR region (Quinolone Resistance Determining Region)
# Residues ~83-87 in E. coli numbering, adjust for M. abscessus
select qrdr, protein and chain A and resi 80-95
color red, qrdr
show sticks, qrdr

# Catalytic tyrosine (Y129 equivalent)
select cat_tyr, protein and chain A and resi 120-135 and resn TYR
color magenta, cat_tyr
show sticks, cat_tyr

# Show DNA if present
show cartoon, reference and chain E+F
color orange, reference and chain E+F

# Binding site sphere
# pseudoatom binding_center, pos=[BINDING_CENTER_X, BINDING_CENTER_Y, BINDING_CENTER_Z]
# show spheres, binding_center
# set sphere_scale, 1.0, binding_center
# color green, binding_center

# Zoom to binding site
zoom reference and resn MFX, 15

# Set nice rendering
bg_color white
set ray_shadows, 0
set cartoon_fancy_helices, 1

# Labels
set label_size, 14
label qrdr and name CA, "QRDR"

print("Binding site visualization loaded")
print("Yellow sticks: Moxifloxacin from reference (PDB 5CDQ)")
print("Red: QRDR region (fluoroquinolone resistance mutations)")
print("Magenta: Catalytic tyrosine")
