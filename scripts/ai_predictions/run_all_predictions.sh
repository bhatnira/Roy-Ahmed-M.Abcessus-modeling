#!/bin/bash
# =============================================================================
# Master Script - Run All AI Structure Predictions
# Utilizes: 24 CPU cores, 124GB RAM, 48GB GPU (RTX PRO 5000 Blackwell)
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
LOG_DIR="$PROJECT_DIR/logs/ai_predictions"

# Create log directory
mkdir -p "$LOG_DIR"

echo "============================================================"
echo "AI Structure Prediction Suite"
echo "============================================================"
echo "Date: $(date)"
echo "Project: M. abscessus DNA Gyrase"
echo ""
echo "Resources:"
echo "  CPU: $(nproc) cores"
echo "  RAM: $(free -h | awk '/^Mem:/ {print $2}')"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "  GPU: None detected"
echo ""
echo "Methods to run:"
echo "  1. ESMFold (Meta AI)"
echo "  2. ColabFold (AlphaFold2 + MMseqs2)"
echo "  3. AlphaFold2 (DeepMind) - if databases available"
echo "  4. RoseTTAFold (Baker Lab) - if installed"
echo ""
echo "============================================================"

# Function to check if a method is available
check_method() {
    local method=$1
    case $method in
        esmfold)
            python -c "import esm" 2>/dev/null && return 0 || return 1
            ;;
        colabfold)
            command -v colabfold_batch >/dev/null 2>&1 && return 0 || return 1
            ;;
        alphafold)
            [ -d "$HOME/alphafold" ] && [ -d "$HOME/alphafold_data" ] && return 0 || return 1
            ;;
        rosettafold)
            [ -d "$HOME/RoseTTAFold" ] && return 0 || return 1
            ;;
    esac
}

# Run ESMFold (fastest, single GPU)
run_esmfold() {
    echo "[$(date '+%H:%M:%S')] Starting ESMFold..."
    bash "$SCRIPT_DIR/run_esmfold_local.sh" > "$LOG_DIR/esmfold.log" 2>&1
    echo "[$(date '+%H:%M:%S')] ESMFold complete."
}

# Run ColabFold (best balance of speed/accuracy)
run_colabfold() {
    echo "[$(date '+%H:%M:%S')] Starting ColabFold..."
    bash "$SCRIPT_DIR/run_colabfold_local.sh" > "$LOG_DIR/colabfold.log" 2>&1
    echo "[$(date '+%H:%M:%S')] ColabFold complete."
}

# Run AlphaFold2 (highest accuracy, requires databases)
run_alphafold() {
    echo "[$(date '+%H:%M:%S')] Starting AlphaFold2..."
    bash "$SCRIPT_DIR/run_alphafold2_local.sh" > "$LOG_DIR/alphafold2.log" 2>&1
    echo "[$(date '+%H:%M:%S')] AlphaFold2 complete."
}

# Run RoseTTAFold
run_rosettafold() {
    echo "[$(date '+%H:%M:%S')] Starting RoseTTAFold..."
    bash "$SCRIPT_DIR/run_rosettafold_local.sh" > "$LOG_DIR/rosettafold.log" 2>&1
    echo "[$(date '+%H:%M:%S')] RoseTTAFold complete."
}

# Main execution
echo "Starting predictions at $(date)"
echo ""

# Since GPU can't be shared efficiently, run sequentially but use all CPU for each
# For methods that support it, enable multi-threading

export OMP_NUM_THREADS=24
export MKL_NUM_THREADS=24
export OPENBLAS_NUM_THREADS=24

# Run predictions
# Note: These run sequentially because they share the GPU
# Each uses all available CPU cores for data processing

echo "============================================================"
echo "PHASE 1: ESMFold (fastest prediction)"
echo "============================================================"
run_esmfold &
ESMFOLD_PID=$!

# Wait for ESMFold to finish before starting next GPU task
wait $ESMFOLD_PID

echo ""
echo "============================================================"
echo "PHASE 2: ColabFold (AlphaFold2 with MMseqs2)"
echo "============================================================"
run_colabfold &
COLABFOLD_PID=$!

wait $COLABFOLD_PID

# Check if AlphaFold2 databases are available
if [ -d "$HOME/alphafold_data" ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 3: AlphaFold2 (full database search)"
    echo "============================================================"
    run_alphafold &
    ALPHAFOLD_PID=$!
    wait $ALPHAFOLD_PID
else
    echo ""
    echo "SKIPPED: AlphaFold2 (databases not found at ~/alphafold_data)"
fi

# Check if RoseTTAFold is available
if [ -d "$HOME/RoseTTAFold" ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 4: RoseTTAFold"
    echo "============================================================"
    run_rosettafold &
    ROSETTAFOLD_PID=$!
    wait $ROSETTAFOLD_PID
else
    echo ""
    echo "SKIPPED: RoseTTAFold (not installed at ~/RoseTTAFold)"
fi

echo ""
echo "============================================================"
echo "ALL PREDICTIONS COMPLETE"
echo "============================================================"
echo "Finished at: $(date)"
echo ""
echo "Results directories:"
echo "  ESMFold:     $PROJECT_DIR/output/ai_predictions/esmfold/"
echo "  ColabFold:   $PROJECT_DIR/output/ai_predictions/colabfold/"
echo "  AlphaFold2:  $PROJECT_DIR/output/ai_predictions/alphafold2/"
echo "  RoseTTAFold: $PROJECT_DIR/output/ai_predictions/rosettafold/"
echo ""
echo "Logs: $LOG_DIR/"
echo ""
echo "Next step: Run comparison analysis"
echo "  python scripts/ai_predictions/compare_predictions.py"
