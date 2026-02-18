#!/bin/bash
# =============================================================================
# RosettaCM for M. abscessus DNA Gyrase Heterotetramer (GyrA2-GyrB2)
# Using locally compiled Rosetta from GitHub source
# =============================================================================

set -e

# Paths
ROSETTA_DIR="/Users/nb/Desktop/rosetta"
ROSETTA_BIN="${ROSETTA_DIR}/source/bin"
ROSETTA_DB="${ROSETTA_DIR}/database"
WORK_DIR="/Users/nb/Desktop/rosetta-cm"

# Find the exact binary name (varies by platform)
ROSETTA_SCRIPTS=$(ls ${ROSETTA_BIN}/rosetta_scripts.* 2>/dev/null | head -1)
PARTIAL_THREAD=$(ls ${ROSETTA_BIN}/partial_thread.* 2>/dev/null | head -1)

if [[ -z "$ROSETTA_SCRIPTS" ]]; then
    echo "ERROR: rosetta_scripts binary not found in ${ROSETTA_BIN}"
    echo "Build may still be in progress. Check with:"
    echo "  ls ${ROSETTA_BIN}/"
    exit 1
fi

if [[ -z "$PARTIAL_THREAD" ]]; then
    echo "ERROR: partial_thread binary not found in ${ROSETTA_BIN}"
    exit 1
fi

echo "=============================================="
echo "RosettaCM: M. abscessus DNA Gyrase Modeling"
echo "=============================================="
echo "Rosetta binaries:"
echo "  rosetta_scripts: $ROSETTA_SCRIPTS"
echo "  partial_thread:  $PARTIAL_THREAD"
echo "Database: $ROSETTA_DB"
echo ""

# Create output directory
OUTPUT_DIR="${WORK_DIR}/output"
mkdir -p ${OUTPUT_DIR}

# Change to work directory
cd ${WORK_DIR}

echo "Step 1: Threading target sequence onto template..."
echo "----------------------------------------"

# Thread each chain separately
# Chain A (GyrA_chain1) -> template chain A
${PARTIAL_THREAD} \
    -database ${ROSETTA_DB} \
    -in:file:fasta input/target_for_5bs8.fasta \
    -in:file:template_pdb input/templates/5bs8_protein.pdb \
    -in:file:alignment input/alignments/gyrase_5bs8.grishin \
    -out:file:silent ${OUTPUT_DIR}/threaded.silent \
    -out:silent_struct_type binary \
    -overwrite

echo ""
echo "Step 2: Running RosettaCM hybridization..."
echo "----------------------------------------"

# Number of output structures
NSTRUCT=5  # Start with 5 for testing, increase to 50-100 for production

${ROSETTA_SCRIPTS} \
    -database ${ROSETTA_DB} \
    -parser:protocol xml/rosettacm_gyrase.xml \
    -in:file:fasta input/target_for_5bs8.fasta \
    -in:file:template_pdb input/templates/5bs8_protein.pdb \
    -in:file:alignment input/alignments/gyrase_5bs8.grishin \
    -nstruct ${NSTRUCT} \
    -out:path:pdb ${OUTPUT_DIR} \
    -out:prefix gyrase_model_ \
    -overwrite \
    -relax:default_repeats 2 \
    -default_max_cycles 200

echo ""
echo "=============================================="
echo "RosettaCM modeling complete!"
echo "Output structures in: ${OUTPUT_DIR}"
echo "=============================================="

# List output files
echo ""
echo "Generated files:"
ls -la ${OUTPUT_DIR}/*.pdb 2>/dev/null || echo "No PDB files generated yet"
