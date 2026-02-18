#!/bin/bash
#
# Run RosettaCM for dimeric homology modeling
#
# This script performs:
# 1. Threading target sequence onto templates
# 2. Running RosettaCM hybridization to generate models
#
# Prerequisites:
# - Rosetta compiled with rosetta_scripts, partial_thread
# - Template PDBs prepared
# - Alignments generated
#

# ============================================
# Configuration - MODIFY THESE PATHS
# ============================================

# Path to Rosetta installation
ROSETTA_PATH="${ROSETTA_PATH:-/path/to/rosetta}"
ROSETTA_BIN="${ROSETTA_PATH}/main/source/bin"
ROSETTA_DB="${ROSETTA_PATH}/main/database"

# Number of models to generate
NSTRUCT=100

# Number of CPU cores to use
NCPU=4

# Input/Output paths
WORKING_DIR="$(pwd)"
INPUT_DIR="${WORKING_DIR}/input"
OUTPUT_DIR="${WORKING_DIR}/output"
TEMPLATE_DIR="${INPUT_DIR}/prepared_templates"
ALIGNMENT_DIR="${INPUT_DIR}/alignments"

# Target sequence
TARGET_FASTA="${INPUT_DIR}/target.fasta"

# ============================================
# Check prerequisites
# ============================================

check_rosetta() {
    if [ ! -d "${ROSETTA_PATH}" ]; then
        echo "ERROR: Rosetta not found at ${ROSETTA_PATH}"
        echo "Please set ROSETTA_PATH environment variable"
        echo "Example: export ROSETTA_PATH=/path/to/rosetta"
        exit 1
    fi
    
    if [ ! -f "${ROSETTA_BIN}/partial_thread.default.linuxgccrelease" ] && \
       [ ! -f "${ROSETTA_BIN}/partial_thread.default.macosclangrelease" ]; then
        echo "WARNING: partial_thread executable not found"
        echo "Please check your Rosetta installation"
    fi
}

# ============================================
# Step 1: Create threaded models
# ============================================

run_threading() {
    echo "================================================"
    echo "Step 1: Threading target onto templates"
    echo "================================================"
    
    mkdir -p "${OUTPUT_DIR}/threaded"
    
    # Find the correct executable suffix
    if [ -f "${ROSETTA_BIN}/partial_thread.default.macosclangrelease" ]; then
        SUFFIX="macosclangrelease"
    else
        SUFFIX="linuxgccrelease"
    fi
    
    PARTIAL_THREAD="${ROSETTA_BIN}/partial_thread.default.${SUFFIX}"
    
    # Thread each template
    for alignment in ${ALIGNMENT_DIR}/*.grishin; do
        if [ ! -f "$alignment" ]; then
            continue
        fi
        
        base=$(basename "$alignment" .grishin)
        template_name=$(echo "$base" | sed 's/alignment_//')
        template_pdb="${TEMPLATE_DIR}/${template_name}_clean.pdb"
        
        if [ ! -f "$template_pdb" ]; then
            echo "WARNING: Template PDB not found: $template_pdb"
            continue
        fi
        
        echo "Threading: $template_name"
        
        ${PARTIAL_THREAD} \
            -database ${ROSETTA_DB} \
            -in:file:fasta ${TARGET_FASTA} \
            -in:file:alignment ${alignment} \
            -in:file:template_pdb ${template_pdb} \
            -out:file:silent "${OUTPUT_DIR}/threaded/${template_name}.silent" \
            -ignore_unrecognized_res true \
            -ignore_zero_occupancy false
            
    done
    
    echo "Threading complete."
}

# ============================================
# Step 2: Generate RosettaCM XML protocol
# ============================================

generate_xml() {
    echo "================================================"
    echo "Step 2: Generating RosettaCM XML protocol"
    echo "================================================"
    
    # Find all threaded models
    TEMPLATE_LIST=""
    for silent in ${OUTPUT_DIR}/threaded/*.silent; do
        if [ -f "$silent" ]; then
            TEMPLATE_LIST="${TEMPLATE_LIST}$(basename $silent .silent),"
        fi
    done
    TEMPLATE_LIST=${TEMPLATE_LIST%,}  # Remove trailing comma
    
    cat > "${OUTPUT_DIR}/rosettacm_dimer.xml" << 'XMLEOF'
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="stage1" weights="score3" symmetric="0">
            <Reweight scoretype="hbond_sr_bb" weight="1.0"/>
            <Reweight scoretype="hbond_lr_bb" weight="1.0"/>
            <Reweight scoretype="atom_pair_constraint" weight="1.0"/>
        </ScoreFunction>
        <ScoreFunction name="stage2" weights="score4_smooth_cart" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
        <ScoreFunction name="fullatom" weights="ref2015_cart" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
    </SCOREFXNS>
    
    <RESIDUE_SELECTORS>
        <Chain name="chainA" chains="A"/>
        <Chain name="chainB" chains="B"/>
    </RESIDUE_SELECTORS>
    
    <MOVERS>
        <!-- Stage 1: Fragment insertion and centroid-level minimization -->
        <Hybridize name="hybridize" 
            stage1_scorefxn="stage1"
            stage2_scorefxn="stage2"
            fa_scorefxn="fullatom"
            batch="1"
            stage1_increase_cycles="1.0"
            stage2_increase_cycles="1.0"
            linmin_only="0">
            <Template pdb="%%template1%%" cst_file="AUTO" weight="1.0"/>
        </Hybridize>
        
        <!-- Final relaxation -->
        <FastRelax name="relax" scorefxn="fullatom" repeats="1"/>
    </MOVERS>
    
    <PROTOCOLS>
        <Add mover="hybridize"/>
        <Add mover="relax"/>
    </PROTOCOLS>
    
    <OUTPUT scorefxn="fullatom"/>
</ROSETTASCRIPTS>
XMLEOF

    echo "XML protocol written to: ${OUTPUT_DIR}/rosettacm_dimer.xml"
}

# ============================================
# Step 3: Run RosettaCM hybridization
# ============================================

run_rosettacm() {
    echo "================================================"
    echo "Step 3: Running RosettaCM hybridization"
    echo "================================================"
    
    mkdir -p "${OUTPUT_DIR}/models"
    
    # Find the correct executable suffix
    if [ -f "${ROSETTA_BIN}/rosetta_scripts.default.macosclangrelease" ]; then
        SUFFIX="macosclangrelease"
    else
        SUFFIX="linuxgccrelease"
    fi
    
    ROSETTA_SCRIPTS="${ROSETTA_BIN}/rosetta_scripts.default.${SUFFIX}"
    
    # Collect all threaded PDBs
    THREADED_PDBS=""
    for silent in ${OUTPUT_DIR}/threaded/*.silent; do
        if [ -f "$silent" ]; then
            THREADED_PDBS="${THREADED_PDBS} -in:file:silent ${silent}"
        fi
    done
    
    # Run RosettaCM
    ${ROSETTA_SCRIPTS} \
        -database ${ROSETTA_DB} \
        -parser:protocol "${OUTPUT_DIR}/rosettacm_dimer.xml" \
        -in:file:fasta ${TARGET_FASTA} \
        ${THREADED_PDBS} \
        -nstruct ${NSTRUCT} \
        -out:path:pdb "${OUTPUT_DIR}/models" \
        -out:file:scorefile "${OUTPUT_DIR}/scores.sc" \
        -relax:default_repeats 1 \
        -default_max_cycles 200 \
        -ignore_unrecognized_res true \
        -ignore_zero_occupancy false \
        -out:suffix "_model"
        
    echo "RosettaCM complete. Models in: ${OUTPUT_DIR}/models"
}

# ============================================
# Alternative: Run with MPI for parallel execution
# ============================================

run_rosettacm_mpi() {
    echo "================================================"
    echo "Running RosettaCM with MPI (${NCPU} cores)"
    echo "================================================"
    
    mkdir -p "${OUTPUT_DIR}/models"
    
    # Find the MPI executable
    if [ -f "${ROSETTA_BIN}/rosetta_scripts.mpi.macosclangrelease" ]; then
        SUFFIX="macosclangrelease"
    else
        SUFFIX="linuxgccrelease"
    fi
    
    ROSETTA_SCRIPTS_MPI="${ROSETTA_BIN}/rosetta_scripts.mpi.${SUFFIX}"
    
    if [ ! -f "${ROSETTA_SCRIPTS_MPI}" ]; then
        echo "MPI executable not found, falling back to serial execution"
        run_rosettacm
        return
    fi
    
    # Collect threaded PDBs
    THREADED_PDBS=""
    for silent in ${OUTPUT_DIR}/threaded/*.silent; do
        if [ -f "$silent" ]; then
            THREADED_PDBS="${THREADED_PDBS} -in:file:silent ${silent}"
        fi
    done
    
    mpirun -np ${NCPU} ${ROSETTA_SCRIPTS_MPI} \
        -database ${ROSETTA_DB} \
        -parser:protocol "${OUTPUT_DIR}/rosettacm_dimer.xml" \
        -in:file:fasta ${TARGET_FASTA} \
        ${THREADED_PDBS} \
        -nstruct ${NSTRUCT} \
        -out:path:pdb "${OUTPUT_DIR}/models" \
        -out:file:scorefile "${OUTPUT_DIR}/scores.sc" \
        -relax:default_repeats 1 \
        -ignore_unrecognized_res true \
        -out:suffix "_model"
}

# ============================================
# Main execution
# ============================================

main() {
    echo "================================================"
    echo "RosettaCM Dimeric Homology Modeling Pipeline"
    echo "================================================"
    echo ""
    echo "Working directory: ${WORKING_DIR}"
    echo "Rosetta path: ${ROSETTA_PATH}"
    echo "Number of models: ${NSTRUCT}"
    echo ""
    
    # Check Rosetta installation
    check_rosetta
    
    # Run pipeline
    run_threading
    generate_xml
    
    # Use MPI if available, otherwise serial
    if command -v mpirun &> /dev/null; then
        run_rosettacm_mpi
    else
        run_rosettacm
    fi
    
    echo ""
    echo "================================================"
    echo "Pipeline complete!"
    echo "================================================"
    echo ""
    echo "Output files:"
    echo "  - Models: ${OUTPUT_DIR}/models/"
    echo "  - Scores: ${OUTPUT_DIR}/scores.sc"
    echo ""
    echo "Next steps:"
    echo "  1. Analyze scores: python scripts/analyze_models.py"
    echo "  2. Select best models based on total_score"
    echo "  3. Validate with molecular dynamics or experimental data"
}

# Run main function
main "$@"
