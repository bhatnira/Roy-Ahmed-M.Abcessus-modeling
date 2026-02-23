# PyMOL script: Sequence alignment-based binding site conservation
# Simpler version using template structures only
# Run: pymol compare_binding_sites.pml

reinitialize

# ============================================================
# LOAD STRUCTURES
# ============================================================
# Load M. tuberculosis reference (5BS8)
load input/templates/5bs8_chainA.pdb, mtb_gyrA
load input/templates/5bs8_chainB.pdb, mtb_gyrB

# Load M. abscessus threaded models
load input/templates/mabs_gyrA_threaded.pdb, mabs_gyrA
load input/templates/mabs_gyrB_threaded.pdb, mabs_gyrB

# ============================================================
# INITIAL SETUP
# ============================================================
hide everything
bg_color white

# Show cartoon representation
show cartoon, all

# Superimpose M. abscessus onto M. tuberculosis
super mabs_gyrA, mtb_gyrA
super mabs_gyrB, mtb_gyrB

# ============================================================
# COLOR BY CONSERVATION
# ============================================================
# MTB = cyan, M.abs = green
color cyan, mtb_gyrA
color lightblue, mtb_gyrB
color forest, mabs_gyrA
color limon, mabs_gyrB

# ============================================================
# QRDR REGION (Quinolone Resistance Determining Region)
# ============================================================
# This is the main fluoroquinolone binding site

# M. tuberculosis QRDR (residues 70-110)
select mtb_qrdr_full, mtb_gyrA and resi 70-110
color red, mtb_qrdr_full

# Key binding residues in MTB
select mtb_binding, mtb_gyrA and resi 88+90+91+94
show sticks, mtb_binding
color yellow, mtb_binding
set stick_radius, 0.3, mtb_binding

# M. abscessus QRDR equivalent 
# Based on alignment: MTB 70 ~ Mabs 86, so offset is +16
select mabs_qrdr_full, mabs_gyrA and resi 86-126
color magenta, mabs_qrdr_full

# Key binding residues in M. abscessus
select mabs_binding, mabs_gyrA and resi 104+106+107+110
show sticks, mabs_binding
color orange, mabs_binding
set stick_radius, 0.3, mabs_binding

# ============================================================
# SEQUENCE ALIGNMENT VISUALIZATION
# Based on pairwise alignment from Grishin file
# ============================================================
# 
# MTB:     ...YHPHGDASIYDSLVRMAQPWSLRYPLVDGQGNFGSPGNDPP...
# M.abs:   ...YHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGSPGNDPA...
#                         *                          *
# 
# Key positions (MTB numbering -> M.abs numbering):
#   88 -> 104 : P -> P (CONSERVED)
#   90 -> 106 : S -> S (CONSERVED)  
#   91 -> 107 : L -> L (CONSERVED)
#   94 -> 110 : P -> P (CONSERVED)
#
# Differences in QRDR:
#   81 -> 97  : S -> T (similar - both polar)
#   110 -> 126: P -> A (different)

# ============================================================
# CONSERVED RESIDUES - show as spheres
# ============================================================
# Position 88/104
select cons_88, (mtb_gyrA and resi 88) or (mabs_gyrA and resi 104)
show spheres, cons_88
color tv_blue, cons_88
set sphere_scale, 0.5, cons_88

# Position 90/106 - SERINE (H-bonding)
select cons_90, (mtb_gyrA and resi 90) or (mabs_gyrA and resi 106)
show spheres, cons_90
color tv_green, cons_90
set sphere_scale, 0.5, cons_90

# Position 91/107 - LEUCINE (hydrophobic)
select cons_91, (mtb_gyrA and resi 91) or (mabs_gyrA and resi 107)
show spheres, cons_91
color tv_yellow, cons_91
set sphere_scale, 0.5, cons_91

# Position 94/110 - critical for Mg2+ coordination
select cons_94, (mtb_gyrA and resi 94) or (mabs_gyrA and resi 110)
show spheres, cons_94
color hotpink, cons_94
set sphere_scale, 0.5, cons_94

# Group all conserved binding residues
group conserved_binding, cons_88 cons_90 cons_91 cons_94

# ============================================================
# DIFFERENT RESIDUES
# ============================================================
# Position 81/97 - S->T substitution
select diff_81, (mtb_gyrA and resi 81) or (mabs_gyrA and resi 97)
show sticks, diff_81
color white, diff_81

# ============================================================
# LABELS
# ============================================================
set label_size, 12
set label_color, black
set label_position, (0, 0, 3)

# Label conserved residues
label mtb_gyrA and resi 88 and name CA, "P88 (MTB)"
label mtb_gyrA and resi 90 and name CA, "S90 (MTB)"
label mtb_gyrA and resi 91 and name CA, "L91 (MTB)"
label mtb_gyrA and resi 94 and name CA, "P94 (MTB)"

# ============================================================
# VIEWS
# ============================================================
# Create binding site view
zoom mtb_binding, 12

# Store views
view binding_site_close, store
zoom mtb_gyrA, 30
view full_view, store

# ============================================================
# DISPLAY SETTINGS
# ============================================================
set ray_shadows, off
set cartoon_fancy_helices, 1
set cartoon_highlight_color, grey50
set depth_cue, 0

# Transparency for better visualization
set cartoon_transparency, 0.3, mabs_gyrA
set cartoon_transparency, 0.3, mabs_gyrB

# ============================================================
# CALCULATE AND DISPLAY RMSD
# ============================================================
# RMSD of QRDR region
rms_cur mtb_qrdr_full and name CA, mabs_qrdr_full and name CA

# ============================================================
# OUTPUT
# ============================================================
print ""
print "============================================================"
print "  BINDING SITE CONSERVATION: M. tuberculosis vs M. abscessus"
print "============================================================"
print ""
print "Structure Colors:"
print "  Cyan/Light Blue : M. tuberculosis (5BS8)"
print "  Green/Lime      : M. abscessus (threaded model)"
print "  Red/Magenta     : QRDR regions"
print ""
print "Conserved Binding Residues (spheres):"
print "  Blue   : Position 88/104 (Pro)"
print "  Green  : Position 90/106 (Ser) - H-bond donor"
print "  Yellow : Position 91/107 (Leu) - Hydrophobic"
print "  Pink   : Position 94/110 (Pro/Asp equiv)"
print ""
print "Alignment Summary:"
print "  QRDR Identity: 95.1%"
print "  Key binding residues: 100% CONSERVED"
print ""
print "Use these commands:"
print "  view binding_site_close  - zoom to binding site"
print "  view full_view           - see full structure"
print "  set cartoon_transparency, 0, all  - remove transparency"
print "============================================================"
