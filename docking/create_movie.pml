# PyMOL Movie Script: M. abscessus DNA Gyrase Modeling & Docking
# Creates a multi-scene movie showcasing template, model, alignment, and docking

# Initialize
reinitialize
set ray_trace_mode, 1
set antialias, 2
bg_color white

# Label settings for on-screen text
set label_font_id, 7
set label_size, 24
set label_color, black
set label_outline_color, white
set label_position, [0, 0, 0]

# ============================================
# LOAD ALL STRUCTURES
# ============================================

# Template (M. tuberculosis 5BS8)
load ../input/templates/5bs8.pdb, mtb_template

# Generated Model (M. abscessus)
load ../output/tetramer/mabs_gyrase_tetramer_dna_mg.pdb, mabs_model

# Docked poses - M. abscessus
load docked_moxifloxacin.pdbqt, mabs_mox
load docked_ciprofloxacin.pdbqt, mabs_cip
load docked_levofloxacin.pdbqt, mabs_lev

# Docked poses - Template
load template_docked_moxifloxacin.pdbqt, mtb_mox
load template_docked_ciprofloxacin.pdbqt, mtb_cip
load template_docked_levofloxacin.pdbqt, mtb_lev

# ============================================
# INITIAL SETUP - HIDE ALL
# ============================================
hide everything

# ============================================
# SCENE 1: M. TUBERCULOSIS TEMPLATE
# ============================================
# Show template protein
show cartoon, mtb_template and polymer.protein
color orange, mtb_template and chain A+C
color lightorange, mtb_template and chain B+D

# Show template DNA
show cartoon, mtb_template and polymer.nucleic
show sticks, mtb_template and polymer.nucleic
color forest, mtb_template and polymer.nucleic

# Show Mg ions
show spheres, mtb_template and resn MG
color green, mtb_template and resn MG
set sphere_scale, 0.6

# Add scene label
pseudoatom label_scene1, pos=[-80, 80, 0]
label label_scene1, "SCENE 1: M. tuberculosis Template (PDB 5BS8)"
pseudoatom label_desc1a, pos=[-80, 70, 0]
label label_desc1a, "Orange: GyrA/GyrB protein"
pseudoatom label_desc1b, pos=[-80, 60, 0]
label label_desc1b, "Green: DNA + Mg2+ ions"

# Orient and save scene
orient mtb_template
scene 001_template, store, Template: M. tuberculosis (PDB 5BS8)

# ============================================
# SCENE 2: M. ABSCESSUS MODEL
# ============================================
hide everything
delete label_scene1
delete label_desc1a
delete label_desc1b

# Show model protein
show cartoon, mabs_model and polymer.protein
color cyan, mabs_model and chain A+C
color lightblue, mabs_model and chain B+D

# Show model DNA
show cartoon, mabs_model and polymer.nucleic
show sticks, mabs_model and polymer.nucleic
color marine, mabs_model and polymer.nucleic

# Show Mg ions
show spheres, mabs_model and resn MG
color green, mabs_model and resn MG

# Add scene label
pseudoatom label_scene2, pos=[-80, 80, 0]
label label_scene2, "SCENE 2: M. abscessus Model (Generated)"
pseudoatom label_desc2a, pos=[-80, 70, 0]
label label_desc2a, "Cyan: GyrA/GyrB protein"
pseudoatom label_desc2b, pos=[-80, 60, 0]
label label_desc2b, "Marine blue: DNA + Mg2+ ions"

orient mabs_model
scene 002_model, store, Generated Model: M. abscessus Gyrase

# ============================================
# SCENE 3: STRUCTURAL ALIGNMENT
# ============================================
hide everything
delete label_scene2
delete label_desc2a
delete label_desc2b

# Align template to model
super mtb_template and chain A, mabs_model and chain A

# Show both proteins
show cartoon, mabs_model and polymer.protein
show cartoon, mtb_template and polymer.protein
color cyan, mabs_model and polymer.protein
color orange, mtb_template and polymer.protein

# Show DNA from model only (cleaner view)
show cartoon, mabs_model and polymer.nucleic
color marine, mabs_model and polymer.nucleic

# Add scene label
pseudoatom label_scene3, pos=[-80, 80, 0]
label label_scene3, "SCENE 3: Structural Alignment"
pseudoatom label_desc3a, pos=[-80, 70, 0]
label label_desc3a, "Cyan: M. abscessus model"
pseudoatom label_desc3b, pos=[-80, 60, 0]
label label_desc3b, "Orange: M. tuberculosis template"
pseudoatom label_desc3c, pos=[-80, 50, 0]
label label_desc3c, "RMSD: GyrA 1.22A, GyrB 0.75A"

orient
scene 003_alignment, store, Structural Alignment: Model (cyan) vs Template (orange)

# ============================================
# SCENE 4: QRDR BINDING SITE FOCUS
# ============================================
delete label_scene3
delete label_desc3a
delete label_desc3b
delete label_desc3c

# Highlight QRDR region
select qrdr_mabs, mabs_model and chain A and resi 80-95
select qrdr_mtb, mtb_template and chain A and resi 83-98
color red, qrdr_mabs
color salmon, qrdr_mtb
show sticks, qrdr_mabs or qrdr_mtb

# Add scene label
pseudoatom label_scene4, pos=[-30, 40, 0]
label label_scene4, "SCENE 4: QRDR Binding Site"
pseudoatom label_desc4a, pos=[-30, 35, 0]
label label_desc4a, "Red: Quinolone Resistance Region"
pseudoatom label_desc4b, pos=[-30, 30, 0]
label label_desc4b, "Residues 80-95 (GyrA)"
pseudoatom label_desc4c, pos=[-30, 25, 0]
label label_desc4c, "Mutations here cause resistance"

zoom qrdr_mabs, 15
scene 004_qrdr, store, QRDR Region: Fluoroquinolone Binding Site

# ============================================
# SCENE 5: M. ABSCESSUS DOCKING
# ============================================
hide everything
delete label_scene4
delete label_desc4a
delete label_desc4b
delete label_desc4c

# Show M. abscessus model
show cartoon, mabs_model and polymer.protein
color gray80, mabs_model and polymer.protein

# Show DNA
show cartoon, mabs_model and polymer.nucleic
show sticks, mabs_model and polymer.nucleic
color marine, mabs_model and polymer.nucleic

# Show QRDR
color salmon, qrdr_mabs
show sticks, qrdr_mabs

# Show docked ligands
show sticks, mabs_mox
show sticks, mabs_cip
show sticks, mabs_lev
color yellow, mabs_mox
color green, mabs_cip
color magenta, mabs_lev

# Add scene label
pseudoatom label_scene5, pos=[-30, 50, 0]
label label_scene5, "SCENE 5: M. abscessus Docking"
pseudoatom label_desc5a, pos=[-30, 45, 0]
label label_desc5a, "Yellow: Moxifloxacin -7.91 kcal/mol"
pseudoatom label_desc5b, pos=[-30, 40, 0]
label label_desc5b, "Green: Ciprofloxacin -7.36 kcal/mol"
pseudoatom label_desc5c, pos=[-30, 35, 0]
label label_desc5c, "Magenta: Levofloxacin -7.65 kcal/mol"

zoom mabs_mox, 15
scene 005_mabs_docking, store, M. abscessus Docking: Mox(-7.91) Cip(-7.36) Lev(-7.65)

# ============================================
# SCENE 6: M. TUBERCULOSIS DOCKING
# ============================================
hide everything
delete label_scene5
delete label_desc5a
delete label_desc5b
delete label_desc5c

# Show template
show cartoon, mtb_template and polymer.protein
color gray60, mtb_template and polymer.protein

# Show DNA
show cartoon, mtb_template and polymer.nucleic
show sticks, mtb_template and polymer.nucleic
color forest, mtb_template and polymer.nucleic

# Show QRDR
color salmon, qrdr_mtb
show sticks, qrdr_mtb

# Show docked ligands
show sticks, mtb_mox
show sticks, mtb_cip
show sticks, mtb_lev
color wheat, mtb_mox
color palegreen, mtb_cip
color lightpink, mtb_lev

# Add scene label
pseudoatom label_scene6, pos=[-30, 50, 0]
label label_scene6, "SCENE 6: M. tuberculosis Docking"
pseudoatom label_desc6a, pos=[-30, 45, 0]
label label_desc6a, "Wheat: Moxifloxacin -7.24 kcal/mol"
pseudoatom label_desc6b, pos=[-30, 40, 0]
label label_desc6b, "Pale green: Ciprofloxacin -7.27 kcal/mol"
pseudoatom label_desc6c, pos=[-30, 35, 0]
label label_desc6c, "Light pink: Levofloxacin -7.31 kcal/mol"

zoom mtb_mox, 15
scene 006_mtb_docking, store, M. tuberculosis Docking: Mox(-7.24) Cip(-7.27) Lev(-7.31)

# ============================================
# SCENE 7: SIDE-BY-SIDE DOCKING COMPARISON
# ============================================
hide everything
delete label_scene6
delete label_desc6a
delete label_desc6b
delete label_desc6c

# Show both proteins
show cartoon, mabs_model and polymer.protein
show cartoon, mtb_template and polymer.protein
color cyan, mabs_model and polymer.protein
color orange, mtb_template and polymer.protein

# Show all ligands
show sticks, mabs_mox or mabs_cip or mabs_lev
show sticks, mtb_mox or mtb_cip or mtb_lev

# Color M. abscessus ligands (bold)
color yellow, mabs_mox
color green, mabs_cip
color magenta, mabs_lev

# Color template ligands (pastel)
color wheat, mtb_mox
color palegreen, mtb_cip
color lightpink, mtb_lev

# Show DNA
show cartoon, mabs_model and polymer.nucleic
color marine, mabs_model and polymer.nucleic

# Add scene label
pseudoatom label_scene7, pos=[-30, 55, 0]
label label_scene7, "SCENE 7: Docking Comparison"
pseudoatom label_desc7a, pos=[-30, 50, 0]
label label_desc7a, "M. abscessus (cyan) vs M. tuberculosis (orange)"
pseudoatom label_desc7b, pos=[-30, 45, 0]
label label_desc7b, "Bold colors: M. abscessus ligands"
pseudoatom label_desc7c, pos=[-30, 40, 0]
label label_desc7c, "Pastel colors: M. tuberculosis ligands"

zoom mabs_mox, 18
scene 007_comparison, store, Docking Comparison: M. abscessus vs M. tuberculosis

# ============================================
# SCENE 8: FINAL OVERVIEW
# ============================================
hide everything
delete label_scene7
delete label_desc7a
delete label_desc7b
delete label_desc7c

# Show M. abscessus model with all components
show cartoon, mabs_model
color cyan, mabs_model and polymer.protein and chain A+C
color lightblue, mabs_model and polymer.protein and chain B+D
color marine, mabs_model and polymer.nucleic
show sticks, mabs_model and polymer.nucleic
show spheres, mabs_model and resn MG
color green, mabs_model and resn MG

# Show best docked ligand (moxifloxacin)
show sticks, mabs_mox
show spheres, mabs_mox
set sphere_scale, 0.3, mabs_mox
color yellow, mabs_mox

# Add scene label
pseudoatom label_scene8, pos=[-80, 80, 0]
label label_scene8, "SCENE 8: Final Model"
pseudoatom label_desc8a, pos=[-80, 70, 0]
label label_desc8a, "M. abscessus DNA Gyrase + Moxifloxacin"
pseudoatom label_desc8b, pos=[-80, 60, 0]
label label_desc8b, "Docking score: -7.91 kcal/mol"
pseudoatom label_desc8c, pos=[-80, 50, 0]
label label_desc8c, "Complete complex ready for analysis"

orient mabs_model
scene 008_final, store, Final: M. abscessus Gyrase with Moxifloxacin

# ============================================
# CREATE MOVIE
# ============================================

# Movie settings
mset 1 x720
# 30 fps, so 720 frames = 24 seconds

# Scene transitions with rotation
# Scene 1: Template (frames 1-90, 3 sec)
frame 1
scene 001_template, animate=0
mview store, object=label_scene1
mview store, object=label_desc1a
mview store, object=label_desc1b
mview store

frame 45
turn y, 180
mview store

frame 90
turn y, 180
mview store

# Scene 2: Model (frames 91-180, 3 sec)
frame 91
scene 002_model, animate=0
mview store

frame 135
turn y, 180
mview store

frame 180
turn y, 180
mview store

# Scene 3: Alignment (frames 181-270, 3 sec)
frame 181
scene 003_alignment, animate=0
mview store

frame 225
turn y, 180
mview store

frame 270
turn y, 180
mview store

# Scene 4: QRDR focus (frames 271-360, 3 sec)
frame 271
scene 004_qrdr, animate=0
mview store

frame 315
turn y, 180
mview store

frame 360
turn y, 180
mview store

# Scene 5: M. abscessus docking (frames 361-450, 3 sec)
frame 361
scene 005_mabs_docking, animate=0
mview store

frame 405
turn y, 180
mview store

frame 450
turn y, 180
mview store

# Scene 6: Template docking (frames 451-540, 3 sec)
frame 451
scene 006_mtb_docking, animate=0
mview store

frame 495
turn y, 180
mview store

frame 540
turn y, 180
mview store

# Scene 7: Comparison (frames 541-630, 3 sec) 
frame 541
scene 007_comparison, animate=0
mview store

frame 585
turn y, 180
mview store

frame 630
turn y, 180
mview store

# Scene 8: Final (frames 631-720, 3 sec)
frame 631
scene 008_final, animate=0
mview store

frame 675
turn y, 180
mview store

frame 720
turn y, 180
mview store

# Interpolate all keyframes
mview interpolate

# ============================================
# DISPLAY INFO
# ============================================
print ""
print "============================================"
print "MOVIE CREATED: 8 scenes, 24 seconds total"
print "============================================"
print ""
print "SCENES:"
print "  1. M. tuberculosis Template (5BS8)"
print "  2. M. abscessus Generated Model"
print "  3. Structural Alignment (cyan vs orange)"
print "  4. QRDR Binding Site Focus"
print "  5. M. abscessus Docking Results"
print "  6. M. tuberculosis Docking Results"
print "  7. Side-by-Side Docking Comparison"
print "  8. Final Overview with Moxifloxacin"
print ""
print "CONTROLS:"
print "  - Press PLAY button or type 'mplay' to play"
print "  - Press STOP or type 'mstop' to stop"
print "  - Use slider to navigate frames"
print ""
print "TO EXPORT MOVIE:"
print "  File > Export Movie As > MPEG..."
print "  Or: mpng movie_frames/frame (saves PNG sequence)"
print ""
