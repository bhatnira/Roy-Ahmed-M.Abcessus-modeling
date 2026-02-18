#!/bin/bash
# =============================================================================
# RosettaCM for M. abscessus DNA Gyrase with DNA, Mg2+, and Ligand Sites
# Complete Complex: GyrA2-GyrB2 + G-segment DNA + Mg ions
# =============================================================================

set -e

# Paths
ROSETTA_DIR="/Users/nb/Desktop/rosetta"
ROSETTA_BIN="${ROSETTA_DIR}/source/bin"
ROSETTA_DB="${ROSETTA_DIR}/database"
WORK_DIR="/Users/nb/Desktop/rosetta-cm"

# Find the exact binary name (varies by platform)
ROSETTA_SCRIPTS=$(ls ${ROSETTA_BIN}/rosetta_scripts.* 2>/dev/null | grep -v "mpi" | head -1)
PARTIAL_THREAD=$(ls ${ROSETTA_BIN}/partial_thread.* 2>/dev/null | grep -v "mpi" | head -1)

if [[ -z "$ROSETTA_SCRIPTS" ]]; then
    echo "ERROR: rosetta_scripts binary not found in ${ROSETTA_BIN}"
    echo "Build may still be in progress. Check with:"
    echo "  ls ${ROSETTA_BIN}/"
    exit 1
fi

echo "=============================================="
echo "RosettaCM: M. abscessus DNA Gyrase Complex"
echo "=============================================="
echo "Components:"
echo "  - Protein: GyrA2-GyrB2 heterotetramer"
echo "  - DNA: G-segment (24bp duplex)"
echo "  - Metals: Mg2+ catalytic ions"
echo ""
echo "Rosetta binary: $ROSETTA_SCRIPTS"
echo "Database: $ROSETTA_DB"
echo ""

# Create output directory
OUTPUT_DIR="${WORK_DIR}/output"
mkdir -p ${OUTPUT_DIR}

# Change to work directory
cd ${WORK_DIR}

# Number of output structures
NSTRUCT=${1:-5}  # Default 5, can be overridden from command line

echo "Generating ${NSTRUCT} models..."
echo "----------------------------------------"

${ROSETTA_SCRIPTS} \
    -database ${ROSETTA_DB} \
    -parser:protocol xml/rosettacm_gyrase_dna.xml \
    -in:file:fasta input/target_with_dna.fasta \
    -nstruct ${NSTRUCT} \
    -out:path:pdb ${OUTPUT_DIR} \
    -out:prefix mabs_gyrase_dna_ \
    -out:file:scorefile ${OUTPUT_DIR}/scores.sc \
    -overwrite \
    -relax:default_repeats 2 \
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -ignore_zero_occupancy false \
    -load_PDB_components true

echo ""
echo "=============================================="
echo "RosettaCM modeling complete!"
echo "=============================================="
echo ""
echo "Output location: ${OUTPUT_DIR}"
echo ""
echo "Generated files:"
ls -la ${OUTPUT_DIR}/mabs_gyrase_dna_*.pdb 2>/dev/null || echo "  (PDB files pending)"
echo ""
echo "Score file: ${OUTPUT_DIR}/scores.sc"
echo ""
echo "To view the best model (lowest score):"
echo "  sort -k2 -n ${OUTPUT_DIR}/scores.sc | head -5"
