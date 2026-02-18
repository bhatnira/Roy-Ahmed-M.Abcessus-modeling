# M. abscessus DNA Gyrase Homology Model

Structural model of **Mycobacterium abscessus DNA Gyrase** built using Rosetta Comparative Modeling.

## Project Summary

| Item | Details |
|------|---------|
| **Target Organism** | *Mycobacterium abscessus* |
| **Target Proteins** | DNA Gyrase subunit A (GyrA) and subunit B (GyrB) |
| **UniProt IDs** | **B1ME58** (GyrA, 839 aa) and **B1ME45** (GyrB, 677 aa) |
| **Template** | PDB **5BS8** - *M. tuberculosis* DNA gyrase-DNA complex (X-ray, 2.56 Å) |
| **Complex** | Heterotetramer (GyrA₂-GyrB₂) with DNA and Mg²⁺ ions |
| **Method** | Rosetta partial_thread + FastRelax |

## Overview

RosettaCM is a powerful tool for building protein structure models using one or more template structures. For this heterotetrameric complex, we:

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

## Model Details

### What was modeled:
| Chain | Protein | Residues Modeled | Position Range | Sequence Match |
|-------|---------|------------------|----------------|----------------|
| A | GyrA | 486 | 17-502 | 100% |
| B | GyrB | 245 | 425-669 | 100% |
| C | GyrA | 486 | 17-502 | 100% |
| D | GyrB | 245 | 425-669 | 100% |
| E-H | DNA | 4 strands | - | From template |

### Rosetta Scores:
- Individual GyrA chain: **-913 REU** (excellent)
- Individual GyrB chain: Relaxed and optimized

## Output Files

| File | Description |
|------|-------------|
| `output/mabs_gyrA_relaxed_*.pdb` | Relaxed GyrA monomer |
| `output/mabs_gyrB_relaxed_*.pdb` | Relaxed GyrB monomer |
| `output/models/mabs_gyrase_dimer_AB_DNA.pdb` | Heterodimer with DNA |
| `output/models/mabs_gyrase_tetramer_A2B2_DNA.pdb` | Full heterotetramer |

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

## References

- Template: Blower TR et al. (2016) Crystal structure of DNA gyrase. PDB: 5BS8
- UniProt: B1ME58 (GyrA), B1ME45 (GyrB) - *Mycobacterium abscessus*
- Rosetta: Song Y et al. (2013) High-resolution comparative modeling with RosettaCM

## Detailed Instructions

See the individual script files for detailed documentation.
