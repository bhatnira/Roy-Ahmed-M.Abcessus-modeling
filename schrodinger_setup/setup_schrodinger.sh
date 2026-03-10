#!/bin/bash
# Schrödinger Suite 2026-1 Environment Setup
# Source this file: source setup_schrodinger.sh

export SCHRODINGER=/opt/schrodinger/schrodinger2026-1
export PATH=$SCHRODINGER:$SCHRODINGER/utilities:$PATH

# Verify installation
echo "Schrödinger Suite 2026-1 configured at: $SCHRODINGER"
echo "Available commands: glide, ligprep, epik, confgen, phase_*, sitemap, ifd, fep_plus"
echo "GUI: maestro"
echo ""
echo "Quick test: $SCHRODINGER/glide -h | head -1"
