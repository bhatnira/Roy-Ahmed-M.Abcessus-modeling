# PyMOL Visualization Script for M. abscessus Gyrase Model vs Template Alignment
# Run with: pymol visualize_structural_alignment.pml
# Or load in PyMOL: @visualize_structural_alignment.pml

# Load the structural alignment file
load output/tetramer/structural_alignment_model_vs_template.pdb, alignment

# Create separate objects for model and template
create mabs_model, chain A+B+C+D
create mtb_template, chain W+X+Y+Z

# Create per-chain objects for detailed analysis
create gyrA_model, chain A
create gyrA_template, chain W
create gyrB_model, chain B
create gyrB_template, chain X

# Hide original combined object
disable alignment

# === VISUALIZATION SETUP ===
hide everything

# Show cartoon representation
show cartoon, mabs_model
show cartoon, mtb_template

# Color scheme: Model in cyan/blue, Template in orange/red
color cyan, mabs_model
color orange, mtb_template

# Alternative: color by chain
# util.cbc mabs_model
# util.cbc mtb_template

# === STRUCTURAL ALIGNMENT ANALYSIS ===
# Superimpose and get RMSD
super mabs_model, mtb_template

# Per-chain superposition for detailed RMSD
super gyrA_model, gyrA_template
super gyrB_model, gyrB_template

# === HIGHLIGHT STRUCTURAL DIFFERENCES ===
# Create selections for conserved core vs variable regions
# (adjust cutoff as needed - lower = stricter)
select conserved_core, mabs_model within 2.0 of mtb_template
select variable_regions, mabs_model and not conserved_core

# Color variable regions differently
color red, variable_regions

# === VISUALIZATION OPTIONS ===

# Option 1: Side-by-side comparison
# translate [50, 0, 0], mtb_template

# Option 2: Overlay (default - structures superimposed)

# Option 3: Show only backbone for cleaner view
# hide cartoon
# show ribbon

# === ACTIVE SITE VISUALIZATION ===
# Highlight catalytic residues (typical for DNA gyrase)
# GyrA: catalytic tyrosine (around residue 122 in E. coli numbering)
# GyrB: ATP binding site

# Select potential active site residues
select active_site_gyrA, (gyrA_model or gyrA_template) and resi 90-130
select active_site_gyrB, (gyrB_model or gyrB_template) and resi 1-50

# Show active site as sticks
show sticks, active_site_gyrA
show sticks, active_site_gyrB

# === FINAL SETUP ===
# Orient view
orient mabs_model or mtb_template

# Set background
bg_color white

# Set rendering quality
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# Enable depth cueing for 3D perception
set depth_cue, 1
set fog_start, 0.4

# === PRINT ALIGNMENT STATISTICS ===
python
from pymol import cmd

print("\n" + "="*60)
print("STRUCTURAL ALIGNMENT QUALITY REPORT")
print("M. abscessus Gyrase Model vs M. tuberculosis Template (5BS8)")
print("="*60 + "\n")

# Overall alignment
result = cmd.super('mabs_model', 'mtb_template')
print(f"Overall Tetramer Alignment:")
print(f"  RMSD: {result[0]:.2f} Angstroms")
print(f"  Aligned atoms: {result[1]}")
print(f"  Rejection cycles: {result[2]}")
print()

# Per-chain alignment
result_A = cmd.super('gyrA_model', 'gyrA_template')
print(f"GyrA (Chain A vs W):")
print(f"  RMSD: {result_A[0]:.2f} Angstroms")
print(f"  Aligned atoms: {result_A[1]}")

result_B = cmd.super('gyrB_model', 'gyrB_template')
print(f"\nGyrB (Chain B vs X):")
print(f"  RMSD: {result_B[0]:.2f} Angstroms")
print(f"  Aligned atoms: {result_B[1]}")

print("\n" + "-"*60)
print("QUALITY INTERPRETATION:")
print("-"*60)
print("  < 1.0 A: Excellent (near identical)")
print("  1.0-2.0 A: Good (very similar)")
print("  2.0-4.0 A: Moderate (same fold, local differences)")
print("  > 4.0 A: Fair/Poor (significant divergence)")
print("="*60 + "\n")
python end

# === INTERACTIVE TIPS ===
# Uncomment commands below as needed:

# Rock animation to see 3D differences
# rock

# Toggle between model and template
# To see model only: disable mtb_template
# To see template only: disable mabs_model
# To see both: enable mabs_model; enable mtb_template

# Save a high-quality image
# ray 2400, 1800
# png structural_alignment_figure.png, dpi=300

# Save session for later
# save structural_alignment_session.pse

print("\nVisualization loaded. Use object panel to toggle visibility.")
print("Commands: rock (animate), ray (render), png filename.png (save image)")
