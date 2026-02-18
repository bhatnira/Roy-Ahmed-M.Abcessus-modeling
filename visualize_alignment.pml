# PyMOL visualization script for M. abscessus gyrase structural alignment
# Run in PyMOL: @visualize_alignment.pml

# Load the aligned structure
load output/tetramer/structural_alignment_model_vs_template.pdb, alignment

# Set up display
bg_color white
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set antialias, 2

# Color the model (M. abscessus) - Blue shades
select model_gyrA, chain A or chain C
select model_gyrB, chain B or chain D
color marine, model_gyrA
color cyan, model_gyrB

# Color the template (M. tuberculosis) - Orange shades
select template_gyrA, chain W or chain Y
select template_gyrB, chain X or chain Z
color orange, template_gyrA
color salmon, template_gyrB

# Show as cartoon
hide everything
show cartoon, all

# Center view
center alignment
zoom alignment

# Create selections for easy comparison
select mabs_model, chain A or chain B or chain C or chain D
select mtb_template, chain W or chain X or chain Y or chain Z

# Print info
print "============================================"
print "M. abscessus DNA Gyrase Structural Alignment"
print "============================================"
print "MODEL (M. abscessus): Chains A, B, C, D (blue/cyan)"
print "TEMPLATE (5BS8):      Chains W, X, Y, Z (orange/salmon)"
print ""
print "Commands:"
print "  show mabs_model      - Show only model"
print "  show mtb_template    - Show only template"
print "  show all             - Show both"
print "============================================"
