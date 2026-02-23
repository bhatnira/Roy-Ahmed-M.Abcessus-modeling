# PyMOL script to generate all presentation figures
# Run: pymol -cq generate_all_figures.pml

# Output directory
import os
os.makedirs("slides/figures", exist_ok=True)

# High resolution settings
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set ray_opaque_background, 1

# ============================================================
# FIGURE 1: Template Structure (5BS8 - M. tuberculosis)
# ============================================================
reinitialize
bg_color white

# Load template
load input/templates/5bs8.pdb, template

hide everything
show cartoon, template

# Color by chain
color cyan, template and chain A
color lightblue, template and chain C
color orange, template and chain B
color lightorange, template and chain D
color red, template and chain E
color salmon, template and chain F

# Show DNA
set cartoon_ring_mode, 3
set cartoon_ring_finder, 2

# Add labels
set label_size, 20
set label_color, black

# Zoom and orient
orient template
zoom template

# Ray trace and save
ray 1600, 1200
png slides/figures/fig1_template_5bs8.png, dpi=300

print "Figure 1: Template 5BS8 saved"

# ============================================================
# FIGURE 2: Template Protein Only
# ============================================================
hide cartoon, template and chain E+F
orient template and (chain A or chain B or chain C or chain D)
zoom template and (chain A or chain B or chain C or chain D)

ray 1600, 1200
png slides/figures/fig2_template_protein_only.png, dpi=300

print "Figure 2: Template protein only saved"

# ============================================================
# FIGURE 3: GyrA Monomer Comparison
# ============================================================
reinitialize
bg_color white

load input/templates/5bs8_chainA.pdb, mtb_gyrA
load input/templates/mabs_gyrA_threaded.pdb, mabs_gyrA

hide everything
show cartoon, all

# Superimpose
super mabs_gyrA, mtb_gyrA

# Color
color cyan, mtb_gyrA
color forest, mabs_gyrA

# Transparency for overlay visibility
set cartoon_transparency, 0.4, mabs_gyrA

orient
zoom all

ray 1600, 1200
png slides/figures/fig3_gyrA_comparison.png, dpi=300

print "Figure 3: GyrA monomer comparison saved"

# ============================================================
# FIGURE 4: GyrB Monomer Comparison
# ============================================================
reinitialize
bg_color white

load input/templates/5bs8_chainB.pdb, mtb_gyrB
load input/templates/mabs_gyrB_threaded.pdb, mabs_gyrB

hide everything
show cartoon, all

super mabs_gyrB, mtb_gyrB

color lightblue, mtb_gyrB
color limon, mabs_gyrB

set cartoon_transparency, 0.4, mabs_gyrB

orient
zoom all

ray 1600, 1200
png slides/figures/fig4_gyrB_comparison.png, dpi=300

print "Figure 4: GyrB monomer comparison saved"

# ============================================================
# FIGURE 5: Full Tetramer Model (if exists)
# ============================================================
reinitialize
bg_color white

# Try to load tetramer
python
import os
tetramer_path = "output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"
if os.path.exists(tetramer_path):
    cmd.load(tetramer_path, "mabs_tetramer")
    print("Loaded tetramer model")
else:
    # Build from monomers
    cmd.load("input/templates/mabs_gyrA_threaded.pdb", "gyrA1")
    cmd.load("input/templates/mabs_gyrB_threaded.pdb", "gyrB1")
    print("Loaded monomer models")
python end

hide everything
show cartoon, all

# Color chains
color marine, chain A or gyrA1
color lightblue, chain C
color orange, chain B or gyrB1
color lightorange, chain D

orient
zoom all

ray 1600, 1200
png slides/figures/fig5_mabs_model.png, dpi=300

print "Figure 5: M. abscessus model saved"

# ============================================================
# FIGURE 6: Binding Site Conservation (Close-up)
# ============================================================
reinitialize
bg_color white

load input/templates/5bs8_chainA.pdb, mtb_gyrA
load input/templates/mabs_gyrA_threaded.pdb, mabs_gyrA

hide everything
show cartoon, all

super mabs_gyrA, mtb_gyrA

color cyan, mtb_gyrA
color forest, mabs_gyrA

# QRDR region
select mtb_qrdr, mtb_gyrA and resi 70-110
select mabs_qrdr, mabs_gyrA and resi 86-126
color red, mtb_qrdr
color magenta, mabs_qrdr

# Key binding residues as sticks
select mtb_key, mtb_gyrA and resi 88+90+91+94
select mabs_key, mabs_gyrA and resi 104+106+107+110
show sticks, mtb_key or mabs_key
color yellow, mtb_key
color orange, mabs_key

# Zoom to binding site
zoom mtb_qrdr, 12
center mtb_qrdr

ray 1600, 1200
png slides/figures/fig6_binding_site_closeup.png, dpi=300

print "Figure 6: Binding site close-up saved"

# ============================================================
# FIGURE 7: Sequence Alignment Visualization (Structural)
# ============================================================
# Same as fig6 but different angle
turn y, 90
ray 1600, 1200
png slides/figures/fig7_binding_site_side.png, dpi=300

turn y, -90
turn x, 45
ray 1600, 1200
png slides/figures/fig7b_binding_site_top.png, dpi=300

print "Figure 7: Binding site alternate views saved"

# ============================================================
# FIGURE 8: Docking - Moxifloxacin
# ============================================================
reinitialize
bg_color white

# Load receptor and docked ligand
python
import os
if os.path.exists("docking/docked_moxifloxacin.pdbqt"):
    cmd.load("docking/docked_moxifloxacin.pdbqt", "moxifloxacin")
if os.path.exists("output/tetramer/mabs_gyrase_tetramer_protein_only.pdb"):
    cmd.load("output/tetramer/mabs_gyrase_tetramer_protein_only.pdb", "receptor")
elif os.path.exists("input/templates/mabs_gyrA_threaded.pdb"):
    cmd.load("input/templates/mabs_gyrA_threaded.pdb", "receptor")
python end

hide everything
show cartoon, receptor
show sticks, moxifloxacin
show spheres, moxifloxacin

# Color
color palegreen, receptor
color yellow, moxifloxacin
set stick_radius, 0.2, moxifloxacin
set sphere_scale, 0.25, moxifloxacin

# Highlight binding site
select binding_region, receptor and (resi 86-126 or resi 104-110)
color forest, binding_region

# Zoom to ligand
zoom moxifloxacin, 10

ray 1600, 1200
png slides/figures/fig8_docking_moxifloxacin.png, dpi=300

print "Figure 8: Moxifloxacin docking saved"

# ============================================================
# FIGURE 9: Docking - All Ligands Comparison
# ============================================================
reinitialize
bg_color white

# Load all docked ligands
python
import os
ligands = ["moxifloxacin", "levofloxacin", "ciprofloxacin"]
colors = ["yellow", "magenta", "cyan"]
for lig, col in zip(ligands, colors):
    path = f"docking/docked_{lig}.pdbqt"
    if os.path.exists(path):
        cmd.load(path, lig)
        cmd.color(col, lig)
        print(f"Loaded {lig}")

# Load receptor
if os.path.exists("input/templates/mabs_gyrA_threaded.pdb"):
    cmd.load("input/templates/mabs_gyrA_threaded.pdb", "receptor")
python end

hide everything
show cartoon, receptor
color gray80, receptor

# Show ligands
show sticks, moxifloxacin or levofloxacin or ciprofloxacin
set stick_radius, 0.15

# Binding site
select qrdr, receptor and resi 86-126
color palegreen, qrdr

# Zoom
zoom moxifloxacin or levofloxacin or ciprofloxacin, 8

ray 1600, 1200
png slides/figures/fig9_docking_comparison.png, dpi=300

print "Figure 9: Docking comparison saved"

# ============================================================
# FIGURE 10: Template vs Model Comparison (M.tb vs M.abs)
# ============================================================
reinitialize
bg_color white

# Load both template and model GyrA
load input/templates/5bs8_chainA.pdb, mtb_template
load input/templates/mabs_gyrA_threaded.pdb, mabs_model

hide everything
show cartoon, all

# Align
super mabs_model, mtb_template

# Side by side coloring
color deepteal, mtb_template
color tv_orange, mabs_model

# Show as surface for one
#show surface, mtb_template
#set surface_color, cyan, mtb_template
#set transparency, 0.7, mtb_template

orient
zoom all, 5

ray 1600, 1200
png slides/figures/fig10_template_vs_model.png, dpi=300

print "Figure 10: Template vs Model comparison saved"

# ============================================================
# FIGURE 11: DNA Interaction Site
# ============================================================
reinitialize
bg_color white

# Load template with DNA
load input/templates/5bs8.pdb, complex

hide everything
show cartoon, complex and (chain A or chain C)
show cartoon, complex and (chain E or chain F)

# Color
color cyan, complex and chain A
color lightblue, complex and chain C
color red, complex and chain E
color salmon, complex and chain F

# DNA rendering
set cartoon_ring_mode, 3
set cartoon_ring_finder, 2

# Highlight DNA-binding residues
select dna_binding, complex and chain A and resi 120-145
show sticks, dna_binding
color yellow, dna_binding

orient complex and chain E
zoom complex and chain E, 15

ray 1600, 1200
png slides/figures/fig11_dna_interaction.png, dpi=300

print "Figure 11: DNA interaction site saved"

# ============================================================
# FIGURE 12: Catalytic Site with Mg2+
# ============================================================
reinitialize
bg_color white

load input/templates/5bs8.pdb, complex

hide everything
show cartoon, complex and chain A
color cyan, complex and chain A

# Show Mg if present
select mg_ions, complex and resn MG
show spheres, mg_ions
color green, mg_ions
set sphere_scale, 0.6, mg_ions

# Catalytic tyrosine region
select cat_site, complex and chain A and resi 125-135
show sticks, cat_site
color magenta, cat_site

zoom cat_site, 12

ray 1600, 1200
png slides/figures/fig12_catalytic_site.png, dpi=300

print "Figure 12: Catalytic site saved"

# ============================================================
# Print summary
# ============================================================
print ""
print "=============================================="
print "ALL FIGURES GENERATED SUCCESSFULLY"
print "=============================================="
print "Location: slides/figures/"
print ""
print "Figures created:"
print "  fig1_template_5bs8.png - Full template with DNA"
print "  fig2_template_protein_only.png - Template protein"
print "  fig3_gyrA_comparison.png - GyrA MTB vs M.abs"
print "  fig4_gyrB_comparison.png - GyrB MTB vs M.abs"
print "  fig5_mabs_model.png - M. abscessus model"
print "  fig6_binding_site_closeup.png - Binding site detail"
print "  fig7_binding_site_side.png - Binding site side view"
print "  fig8_docking_moxifloxacin.png - Moxifloxacin docked"
print "  fig9_docking_comparison.png - All ligands"
print "  fig10_template_vs_model.png - Overlay comparison"
print "  fig11_dna_interaction.png - DNA binding site"
print "  fig12_catalytic_site.png - Catalytic site + Mg"
print "=============================================="

quit
