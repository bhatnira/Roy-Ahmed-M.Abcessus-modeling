# Dimeric Homology Modeling with Rosetta CM

This guide walks through the process of building a homology model for a **dimeric protein** using Rosetta Comparative Modeling (RosettaCM).

## Overview

RosettaCM is a powerful tool for building protein structure models using one or more template structures. For dimeric proteins, we need to:

1. Prepare the target sequence (both chains)
2. Identify and prepare template structures
3. Generate sequence alignments
4. Create threading models
5. Run RosettaCM hybridization
6. Analyze and select the best models

## Prerequisites

- Rosetta software suite installed
- Python 3.x with BioPython
- Template PDB structures
- Target sequence in FASTA format

## Directory Structure

```
rosetta-cm/
├── input/
│   ├── target.fasta          # Target sequence (both chains)
│   ├── templates/            # Template PDB files
│   └── alignments/           # Sequence alignments
├── scripts/
│   ├── prepare_templates.py  # Template preparation script
│   ├── generate_alignment.py # Alignment generation
│   └── run_rosettacm.sh      # Main execution script
├── output/
│   └── models/               # Output models
└── README.md
```

## Quick Start

```bash
# 1. Set up directories
./scripts/setup_directories.sh

# 2. Prepare templates
python scripts/prepare_templates.py

# 3. Generate alignments
python scripts/generate_alignment.py

# 4. Run RosettaCM
./scripts/run_rosettacm.sh
```

## Detailed Instructions

See the individual script files for detailed documentation.
