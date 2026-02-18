#!/bin/bash
# ==============================================
# RosettaCM with Threading: M. abscessus DNA Gyrase
# Step 1: Thread sequences onto template
# Step 2: Run RosettaCM hybridize
# Step 3: Merge DNA and Mg2+
# ==============================================

set -e

ROSETTA_BIN="/Users/nb/Desktop/rosetta/source/bin"
ROSETTA_DB="/Users/nb/Desktop/rosetta/database"
WORK_DIR="/Users/nb/Desktop/rosetta-cm"
NSTRUCT=${1:-2}

cd ${WORK_DIR}

echo "=============================================="
echo "RosettaCM: M. abscessus DNA Gyrase Complex"
echo "=============================================="

# Step 1: Use partial_thread to create threaded models
echo ""
echo "Step 1: Threading sequences onto template..."
echo "----------------------------------------"

# Thread GyrA onto chain A
${ROSETTA_BIN}/partial_thread.default.macosclangrelease \
    -database ${ROSETTA_DB} \
    -in:file:fasta input/gyrA_single.fasta \
    -in:file:template_pdb input/templates/5bs8_chainA.pdb \
    -out:file:o output/gyrA_threaded.pdb \
    -ignore_unrecognized_res \
    -ignore_zero_occupancy false 2>&1 | tail -20

echo ""
echo "Threaded model created: output/gyrA_threaded.pdb"
echo ""
echo "Step 2: Running RosettaCM hybridize..."
echo "----------------------------------------"

# Now run RosettaCM with the threaded template
${ROSETTA_BIN}/rosetta_scripts.default.macosclangrelease \
    -database ${ROSETTA_DB} \
    -parser:protocol xml/rosettacm_threaded.xml \
    -in:file:fasta input/gyrA_single.fasta \
    -nstruct ${NSTRUCT} \
    -out:path:pdb output \
    -out:prefix mabs_gyrA_ \
    -out:file:scorefile output/scores.sc \
    -overwrite \
    -relax:default_repeats 2

echo ""
echo "=============================================="
echo "COMPLETE!"
echo "=============================================="
ls -la output/mabs_gyrA_*.pdb 2>/dev/null || echo "Check output directory"
