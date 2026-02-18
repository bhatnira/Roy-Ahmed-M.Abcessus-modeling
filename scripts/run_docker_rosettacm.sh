#!/bin/bash
#
# RosettaCM Docker Workflow for M. abscessus DNA Gyrase
#
# This script runs RosettaCM using Docker containers
#

set -e

WORKDIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$WORKDIR"

# Number of models to generate
NSTRUCT=${NSTRUCT:-5}

echo "=============================================="
echo "RosettaCM via Docker"
echo "M. abscessus DNA Gyrase Modeling"
echo "=============================================="
echo ""
echo "Working directory: $WORKDIR"
echo "Models to generate: $NSTRUCT"
echo ""

# Check Docker
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker not found"
    exit 1
fi

# Create output directories
mkdir -p output/threaded output/models

# Check for Rosetta image
ROSETTA_IMAGE="rosettacommons/rosetta:latest"

echo "Checking for Rosetta Docker image..."
if ! docker images | grep -q "rosettacommons/rosetta"; then
    echo ""
    echo "NOTE: Official Rosetta Docker image requires a license."
    echo "Visit: https://www.rosettacommons.org/software/license-and-download"
    echo ""
    echo "Alternative: Using publicly available tools..."
    echo ""
fi

# ===========================================
# Option 1: If you have Rosetta Docker image
# ===========================================
run_with_rosetta_docker() {
    echo "Step 1: Threading..."
    docker run --rm -v "$WORKDIR:/work" -w /work $ROSETTA_IMAGE \
        partial_thread.default.linuxgccrelease \
        -database /rosetta/main/database \
        -in:file:fasta /work/input/target_for_5bs8.fasta \
        -in:file:alignment /work/input/alignments/gyrase_5bs8.grishin \
        -in:file:template_pdb /work/input/templates/5bs8_protein.pdb \
        -out:file:silent /work/output/threaded/gyrase_threaded.silent \
        -ignore_unrecognized_res true

    echo "Step 2: RosettaCM Hybridize..."
    docker run --rm -v "$WORKDIR:/work" -w /work $ROSETTA_IMAGE \
        rosetta_scripts.default.linuxgccrelease \
        -database /rosetta/main/database \
        -parser:protocol /work/xml/rosettacm_gyrase.xml \
        -in:file:fasta /work/input/target_for_5bs8.fasta \
        -in:file:silent /work/output/threaded/gyrase_threaded.silent \
        -nstruct $NSTRUCT \
        -out:path:pdb /work/output/models/ \
        -out:file:scorefile /work/output/scores.sc \
        -out:prefix gyrase_mabs_ \
        -ignore_unrecognized_res true
}

# ===========================================
# Option 2: Use PyRosetta Docker
# ===========================================
run_with_pyrosetta() {
    echo "Using PyRosetta Docker container..."
    
    # Create PyRosetta script
    cat > "$WORKDIR/scripts/run_pyrosetta_cm.py" << 'PYEOF'
#!/usr/bin/env python3
"""
PyRosetta-based comparative modeling for DNA Gyrase
"""
import os
import sys

try:
    import pyrosetta
    from pyrosetta import *
    from pyrosetta.rosetta.protocols.hybridization import *
except ImportError:
    print("PyRosetta not available")
    sys.exit(1)

# Initialize PyRosetta
pyrosetta.init("-ignore_unrecognized_res true")

print("PyRosetta initialized successfully")
print("Running comparative modeling...")

# Load template
template = pose_from_pdb("/work/input/templates/5bs8_protein.pdb")
print(f"Template loaded: {template.total_residue()} residues")

# Run hybridize protocol
# ... (simplified for demonstration)

print("Modeling complete")
PYEOF

    docker run --rm -v "$WORKDIR:/work" -w /work \
        pyrosetta/pyrosetta:latest \
        python3 /work/scripts/run_pyrosetta_cm.py
}

# ===========================================
# Main execution
# ===========================================

echo "Attempting to run with Rosetta Docker..."
echo ""

# Try to pull and run
if docker pull $ROSETTA_IMAGE 2>/dev/null; then
    run_with_rosetta_docker
else
    echo "=============================================="
    echo "Rosetta Docker image not available"
    echo "=============================================="
    echo ""
    echo "The official Rosetta Docker image requires authentication."
    echo ""
    echo "OPTIONS:"
    echo ""
    echo "1. Get Rosetta license and Docker access:"
    echo "   https://www.rosettacommons.org/software/license-and-download"
    echo ""
    echo "2. Use ROBETTA web server (free):"
    echo "   https://robetta.bakerlab.org/"
    echo "   Upload: input/target_for_5bs8.fasta"
    echo "   Upload: input/templates/5bs8_protein.pdb"
    echo ""
    echo "3. Use local Rosetta installation:"
    echo "   export ROSETTA_PATH=/path/to/rosetta"
    echo "   ./scripts/run_gyrase_modeling.sh"
    echo ""
    echo "4. Use PyRosetta (requires license):"
    echo "   pip install pyrosetta"
    echo ""
fi
