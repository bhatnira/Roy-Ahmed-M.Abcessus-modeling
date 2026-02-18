#!/bin/bash
#
# Quick-start script for RosettaCM dimeric homology modeling
# 
# This is an all-in-one example that demonstrates the complete workflow.
#

echo "=============================================="
echo "RosettaCM Dimeric Homology Modeling Example"
echo "=============================================="
echo ""

# Set up directory structure
mkdir -p input/templates
mkdir -p input/prepared_templates
mkdir -p input/alignments
mkdir -p output/threaded
mkdir -p output/models

# Check for Rosetta
if [ -z "${ROSETTA_PATH}" ]; then
    echo "WARNING: ROSETTA_PATH not set"
    echo "Please set it: export ROSETTA_PATH=/path/to/rosetta"
    echo ""
fi

echo "STEP-BY-STEP WORKFLOW:"
echo "======================"
echo ""
echo "1. PREPARE YOUR INPUT FILES:"
echo "   - Place your target sequence in: input/target.fasta"
echo "   - Place template PDB files in: input/templates/"
echo ""
echo "   Target FASTA format (for homodimer):"
echo "   >ChainA"
echo "   MKXXXXXX..."
echo "   >ChainB"
echo "   MKXXXXXX..."
echo ""
echo "2. PREPARE TEMPLATES:"
echo "   python scripts/prepare_templates.py"
echo ""
echo "3. GENERATE ALIGNMENTS:"
echo "   python scripts/generate_alignment.py"
echo ""
echo "   Or manually create GRISHIN format alignments in:"
echo "   input/alignments/"
echo ""
echo "4. (OPTIONAL) GENERATE CONSTRAINTS:"
echo "   python scripts/generate_constraints.py --template input/templates/template.pdb"
echo ""
echo "5. RUN ROSETTACM:"
echo "   ./scripts/run_rosettacm.sh"
echo ""
echo "   Or manually:"
echo "   \$ROSETTA/partial_thread ... # Threading"
echo "   \$ROSETTA/rosetta_scripts -parser:protocol xml/rosettacm_dimer.xml ..."
echo ""
echo "6. ANALYZE RESULTS:"
echo "   python scripts/analyze_models.py"
echo ""
echo "=============================================="
echo "Example Rosetta Commands"
echo "=============================================="
echo ""
echo "# Threading command:"
echo 'partial_thread.default.linuxgccrelease \'
echo '  -database $ROSETTA_DB \'
echo '  -in:file:fasta input/target.fasta \'
echo '  -in:file:alignment input/alignments/alignment.grishin \'
echo '  -in:file:template_pdb input/templates/template.pdb'
echo ""
echo "# RosettaCM command:"
echo 'rosetta_scripts.default.linuxgccrelease \'
echo '  -database $ROSETTA_DB \'
echo '  -parser:protocol xml/rosettacm_dimer.xml \'
echo '  -in:file:fasta input/target.fasta \'
echo '  -in:file:silent output/threaded/template.silent \'
echo '  -nstruct 100 \'
echo '  -out:path:pdb output/models/'
echo ""
echo "=============================================="
echo ""
echo "For symmetric homodimers, use:"
echo "  xml/rosettacm_symmetric.xml"
echo "  with: input/symmetry_c2.symm"
echo ""
