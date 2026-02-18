#!/bin/bash
# ==============================================
# RosettaCM: M. abscessus DNA Gyrase with DNA
# Two-stage approach:
#   1. Model protein heterotetramer with RosettaCM
#   2. Merge DNA and Mg2+ from template
# ==============================================

set -e

# Configuration
ROSETTA_SCRIPTS="/Users/nb/Desktop/rosetta/source/bin/rosetta_scripts.default.macosclangrelease"
ROSETTA_DB="/Users/nb/Desktop/rosetta/database"
INPUT_DIR="/Users/nb/Desktop/rosetta-cm/input"
OUTPUT_DIR="/Users/nb/Desktop/rosetta-cm/output"
XML_FILE="/Users/nb/Desktop/rosetta-cm/xml/rosettacm_protein.xml"
NSTRUCT=${1:-2}

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=============================================="
echo "RosettaCM: M. abscessus DNA Gyrase Complex"
echo "=============================================="
echo "Stage 1: Model protein heterotetramer (GyrA2-GyrB2)"
echo "Stage 2: Add DNA and Mg2+ from template"
echo ""
echo "Rosetta binary: ${ROSETTA_SCRIPTS}"
echo "Database: ${ROSETTA_DB}"
echo ""
echo "Generating ${NSTRUCT} models..."
echo "----------------------------------------"

# Stage 1: Run RosettaCM for protein only
${ROSETTA_SCRIPTS} \
    -database ${ROSETTA_DB} \
    -parser:protocol ${XML_FILE} \
    -in:file:fasta ${INPUT_DIR}/target_protein_only.fasta \
    -nstruct ${NSTRUCT} \
    -out:path:pdb ${OUTPUT_DIR} \
    -out:prefix mabs_gyrase_protein_ \
    -out:file:scorefile ${OUTPUT_DIR}/scores_protein.sc \
    -overwrite \
    -relax:default_repeats 2 \
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -ignore_zero_occupancy false

echo ""
echo "----------------------------------------"
echo "Stage 1 complete. Protein models generated."
echo ""
echo "Stage 2: Merging DNA and Mg2+ from template..."

# Stage 2: Merge DNA and metals into each model
DNA_MG_FILE="${INPUT_DIR}/templates/5bs8_dna_mg.pdb"

for i in $(seq -w 1 ${NSTRUCT}); do
    PROTEIN_PDB="${OUTPUT_DIR}/mabs_gyrase_protein_${i}.pdb"
    FINAL_PDB="${OUTPUT_DIR}/mabs_gyrase_dna_${i}.pdb"
    
    if [ -f "${PROTEIN_PDB}" ]; then
        echo "  Processing model ${i}..."
        
        # Remove END line, append DNA/Mg, add END
        grep -v "^END" "${PROTEIN_PDB}" > "${FINAL_PDB}"
        cat "${DNA_MG_FILE}" >> "${FINAL_PDB}"
        echo "END" >> "${FINAL_PDB}"
        
        echo "    Created: ${FINAL_PDB}"
    fi
done

echo ""
echo "=============================================="
echo "COMPLETE!"
echo "=============================================="
echo "Output files:"
ls -la ${OUTPUT_DIR}/mabs_gyrase_dna_*.pdb 2>/dev/null || echo "  (check output directory)"
echo ""
echo "Each model contains:"
echo "  - GyrA chains (A, C) - modeled"
echo "  - GyrB chains (B, D) - modeled"
echo "  - DNA G-segment (E, F, G, H) - from template"
echo "  - Mg2+ ions - from template"
