#!/bin/bash
#
# RosettaCM workflow for M. abscessus DNA Gyrase
# Template: 5BS8 (M. tuberculosis)
#
# This script runs the complete homology modeling pipeline
#

set -e

# ============================================
# Configuration
# ============================================

# Rosetta paths - MODIFY THESE
ROSETTA_PATH="${ROSETTA_PATH:-/path/to/rosetta}"
ROSETTA_BIN="${ROSETTA_PATH}/main/source/bin"
ROSETTA_DB="${ROSETTA_PATH}/main/database"

# Detect OS for executable suffix
if [[ "$OSTYPE" == "darwin"* ]]; then
    SUFFIX="macosclangrelease"
else
    SUFFIX="linuxgccrelease"
fi

# Number of models
NSTRUCT=${NSTRUCT:-10}

# Working directory
WORKDIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=============================================="
echo "M. abscessus DNA Gyrase Homology Modeling"
echo "Using RosettaCM with 5BS8 template"
echo "=============================================="
echo ""
echo "Working directory: $WORKDIR"
echo "Rosetta path: $ROSETTA_PATH"
echo "Number of models: $NSTRUCT"
echo ""

# ============================================
# Check prerequisites
# ============================================

check_files() {
    echo "Checking input files..."
    
    if [ ! -f "$WORKDIR/input/target_for_5bs8.fasta" ]; then
        echo "ERROR: Target FASTA not found"
        echo "Run: python3 scripts/setup_gyrase_model.py"
        exit 1
    fi
    
    if [ ! -f "$WORKDIR/input/templates/5bs8_protein.pdb" ]; then
        echo "ERROR: Template PDB not found"
        echo "Creating protein-only template..."
        grep -v "^ATOM.* [EFGH] " "$WORKDIR/input/templates/5bs8.pdb" | \
        grep -v "^HETATM.* [EFGH] " > "$WORKDIR/input/templates/5bs8_protein.pdb"
    fi
    
    if [ ! -f "$WORKDIR/input/alignments/gyrase_5bs8.grishin" ]; then
        echo "ERROR: Alignment file not found"
        echo "Run: python3 scripts/setup_gyrase_model.py"
        exit 1
    fi
    
    echo "All input files present."
}

check_rosetta() {
    if [ ! -d "$ROSETTA_PATH" ]; then
        echo ""
        echo "WARNING: Rosetta not found at $ROSETTA_PATH"
        echo ""
        echo "To run this workflow, you need Rosetta installed."
        echo "Set the ROSETTA_PATH environment variable:"
        echo "  export ROSETTA_PATH=/path/to/rosetta"
        echo ""
        echo "The following commands would be run:"
        echo ""
        print_commands
        exit 0
    fi
}

print_commands() {
    cat << EOF
# Step 1: Threading (partial_thread)
${ROSETTA_BIN}/partial_thread.default.${SUFFIX} \\
    -database ${ROSETTA_DB} \\
    -in:file:fasta ${WORKDIR}/input/target_for_5bs8.fasta \\
    -in:file:alignment ${WORKDIR}/input/alignments/gyrase_5bs8.grishin \\
    -in:file:template_pdb ${WORKDIR}/input/templates/5bs8_protein.pdb \\
    -out:file:silent ${WORKDIR}/output/threaded/gyrase_threaded.silent \\
    -ignore_unrecognized_res true

# Step 2: RosettaCM Hybridization
${ROSETTA_BIN}/rosetta_scripts.default.${SUFFIX} \\
    -database ${ROSETTA_DB} \\
    -parser:protocol ${WORKDIR}/xml/rosettacm_gyrase.xml \\
    -in:file:fasta ${WORKDIR}/input/target_for_5bs8.fasta \\
    -in:file:silent ${WORKDIR}/output/threaded/gyrase_threaded.silent \\
    -nstruct ${NSTRUCT} \\
    -out:path:pdb ${WORKDIR}/output/models/ \\
    -out:file:scorefile ${WORKDIR}/output/scores.sc \\
    -out:prefix gyrase_mabs_ \\
    -ignore_unrecognized_res true \\
    -relax:default_repeats 1

EOF
}

# ============================================
# Run pipeline
# ============================================

run_threading() {
    echo ""
    echo "Step 1: Threading target onto template..."
    echo "----------------------------------------"
    
    mkdir -p "$WORKDIR/output/threaded"
    
    ${ROSETTA_BIN}/partial_thread.default.${SUFFIX} \
        -database ${ROSETTA_DB} \
        -in:file:fasta ${WORKDIR}/input/target_for_5bs8.fasta \
        -in:file:alignment ${WORKDIR}/input/alignments/gyrase_5bs8.grishin \
        -in:file:template_pdb ${WORKDIR}/input/templates/5bs8_protein.pdb \
        -out:file:silent ${WORKDIR}/output/threaded/gyrase_threaded.silent \
        -ignore_unrecognized_res true
    
    echo "Threading complete."
}

run_rosettacm() {
    echo ""
    echo "Step 2: Running RosettaCM..."
    echo "----------------------------"
    
    mkdir -p "$WORKDIR/output/models"
    
    ${ROSETTA_BIN}/rosetta_scripts.default.${SUFFIX} \
        -database ${ROSETTA_DB} \
        -parser:protocol ${WORKDIR}/xml/rosettacm_gyrase.xml \
        -in:file:fasta ${WORKDIR}/input/target_for_5bs8.fasta \
        -in:file:silent ${WORKDIR}/output/threaded/gyrase_threaded.silent \
        -nstruct ${NSTRUCT} \
        -out:path:pdb ${WORKDIR}/output/models/ \
        -out:file:scorefile ${WORKDIR}/output/scores.sc \
        -out:prefix gyrase_mabs_ \
        -ignore_unrecognized_res true \
        -relax:default_repeats 1
    
    echo "RosettaCM complete."
}

# ============================================
# Main
# ============================================

cd "$WORKDIR"

check_files
check_rosetta

echo ""
echo "Starting RosettaCM pipeline..."
echo ""

run_threading
run_rosettacm

echo ""
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  Models: output/models/"
echo "  Scores: output/scores.sc"
echo ""
echo "Analyze results:"
echo "  python3 scripts/analyze_models.py"
echo ""
