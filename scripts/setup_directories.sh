#!/bin/bash
# Setup directory structure for RosettaCM dimeric modeling

mkdir -p input/templates
mkdir -p input/alignments
mkdir -p output/models
mkdir -p scripts

echo "Directory structure created successfully!"
echo ""
echo "Next steps:"
echo "1. Place your target sequence in input/target.fasta"
echo "2. Place template PDB files in input/templates/"
echo "3. Run prepare_templates.py"
