#!/bin/bash
# Glide Docking Workflow for M. abscessus GyrA/B
# Using Schrödinger Glide

SCHRODINGER=/opt/schrodinger/schrodinger2026-1
WORKDIR=$(pwd)

echo "=============================================="
echo "  Glide Docking Workflow"
echo "  Schrödinger Glide 2026-1"
echo "=============================================="

# Step 1: Prepare protein
prepare_protein() {
    echo "[Step 1] Preparing protein with PrepWizard..."
    $SCHRODINGER/utilities/prepwizard \
        "$1" \
        "${1%.pdb}_prep.mae" \
        -fillsidechains \
        -fillloops \
        -disulfides \
        -mse \
        -propka_pH 7.0 \
        -WAIT
    echo "Prepared protein saved to: ${1%.pdb}_prep.mae"
}

# Step 2: Generate receptor grid (requires grid input file)
generate_grid() {
    echo "[Step 2] Generating receptor grid..."
    cat > grid_temp.in << EOF
GRIDFILE    ${1%.mae}_grid.zip
INNERBOX    10,10,10
OUTERBOX    30,30,30
RECEP_FILE  $1
EOF
    $SCHRODINGER/glide grid_temp.in -WAIT
    rm grid_temp.in
    echo "Grid file saved to: ${1%.mae}_grid.zip"
}

# Step 3: Prepare ligands
prepare_ligands() {
    echo "[Step 3] Preparing ligands with LigPrep..."
    $SCHRODINGER/ligprep \
        -ismi "$1" \
        -omae "${1%.smi}_prep.mae" \
        -g \
        -s 1 \
        -i 2 -t 0 \
        -W e,-ph,7.0,-pht,2.0 \
        -WAIT
    echo "Prepared ligands saved to: ${1%.smi}_prep.mae"
}

# Step 4: Run Glide SP docking
dock_sp() {
    echo "[Step 4] Running Glide SP docking..."
    cat > dock_temp.in << EOF
GRIDFILE    $1
LIGANDFILE  $2
PRECISION   SP
POSTDOCK_NPOSE  5
EOF
    $SCHRODINGER/glide dock_temp.in -WAIT -NJOBS 4
    rm dock_temp.in
    echo "Docking complete!"
}

# Step 5: Run Glide XP docking (for rescoring)
dock_xp() {
    echo "[Step 5] Running Glide XP docking..."
    cat > dock_temp.in << EOF
GRIDFILE    $1
LIGANDFILE  $2
PRECISION   XP
WRITE_XP_DESC  TRUE
POSTDOCK_NPOSE  3
EOF
    $SCHRODINGER/glide dock_temp.in -WAIT -NJOBS 4
    rm dock_temp.in
    echo "XP Docking complete!"
}

# Print usage
print_usage() {
    echo ""
    echo "Usage:"
    echo "  ./glide_docking_workflow.sh prepare_protein <protein.pdb>"
    echo "  ./glide_docking_workflow.sh generate_grid <receptor.mae>"
    echo "  ./glide_docking_workflow.sh prepare_ligands <ligands.smi>"
    echo "  ./glide_docking_workflow.sh dock_sp <grid.zip> <ligands.mae>"
    echo "  ./glide_docking_workflow.sh dock_xp <grid.zip> <ligands.mae>"
    echo ""
}

case "$1" in
    prepare_protein)
        prepare_protein "$2"
        ;;
    generate_grid)
        generate_grid "$2"
        ;;
    prepare_ligands)
        prepare_ligands "$2"
        ;;
    dock_sp)
        dock_sp "$2" "$3"
        ;;
    dock_xp)
        dock_xp "$2" "$3"
        ;;
    *)
        print_usage
        ;;
esac
