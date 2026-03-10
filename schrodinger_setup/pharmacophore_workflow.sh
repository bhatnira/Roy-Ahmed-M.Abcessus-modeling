#!/bin/bash
# Pharmacophore Modeling Workflow for M. abscessus GyrA/B
# Using Schrödinger Phase

SCHRODINGER=/opt/schrodinger/schrodinger2026-1
WORKDIR=$(pwd)

echo "=============================================="
echo "  Pharmacophore Modeling Workflow"
echo "  Schrödinger Phase 2026-1"
echo "=============================================="

# Step 1: Prepare ligands (if needed)
prepare_ligands() {
    echo "[Step 1] Preparing ligands..."
    $SCHRODINGER/ligprep \
        -ismi "$1" \
        -omae "${1%.smi}_prep.mae" \
        -g \
        -s 1 \
        -i 0 \
        -WAIT
    echo "Prepared ligands saved to: ${1%.smi}_prep.mae"
}

# Step 2: Create Phase database
create_phase_db() {
    echo "[Step 2] Creating Phase database..."
    $SCHRODINGER/phase_database "$1" import \
        -structure "$2" \
        -sites auto \
        -confs auto \
        -WAIT
    echo "Database created: $1"
}

# Step 3: Find common pharmacophore
find_common_pharmacophore() {
    echo "[Step 3] Finding common pharmacophore features..."
    $SCHRODINGER/phase_find_common \
        -actives "$1" \
        -o "${1%.mae}_hypothesis" \
        -sites 4-6 \
        -match 50 \
        -WAIT
    echo "Hypothesis saved to: ${1%.mae}_hypothesis.phypo"
}

# Step 4: Screen database
screen_database() {
    echo "[Step 4] Screening database with hypothesis..."
    $SCHRODINGER/phase_screen \
        "$1" \
        "$2" \
        "$3" \
        -flex \
        -match 4 \
        -WAIT
    echo "Screen results saved to: ${3}_phasedb_hits.mae"
}

# Step 5: Build 3D-QSAR model
build_qsar() {
    echo "[Step 5] Building 3D-QSAR model..."
    $SCHRODINGER/phase_build_qsar \
        "$1" \
        "$2" \
        -activity "$3" \
        -o "${2%.phypo}_qsar" \
        -WAIT
    echo "QSAR model saved to: ${2%.phypo}_qsar"
}

# Main workflow
print_usage() {
    echo ""
    echo "Usage:"
    echo "  ./pharmacophore_workflow.sh prepare <ligands.smi>"
    echo "  ./pharmacophore_workflow.sh create_db <database.phdb> <structures.mae>"
    echo "  ./pharmacophore_workflow.sh find_common <actives.mae>"
    echo "  ./pharmacophore_workflow.sh screen <database.phdb> <hypothesis.phypo> <jobname>"
    echo "  ./pharmacophore_workflow.sh qsar <training.mae> <hypothesis.phypo> <activity_property>"
    echo ""
}

case "$1" in
    prepare)
        prepare_ligands "$2"
        ;;
    create_db)
        create_phase_db "$2" "$3"
        ;;
    find_common)
        find_common_pharmacophore "$2"
        ;;
    screen)
        screen_database "$2" "$3" "$4"
        ;;
    qsar)
        build_qsar "$2" "$3" "$4"
        ;;
    *)
        print_usage
        ;;
esac
