# M. abscessus DNA Gyrase Homology Model

Structural model of **Mycobacterium abscessus DNA Gyrase** built using Rosetta Comparative Modeling.

## Project Summary

| Item | Details |
|------|---------|
| **Target Organism** | *Mycobacterium abscessus* |
| **Target Proteins** | DNA Gyrase subunit A (GyrA) and subunit B (GyrB) |
| **UniProt IDs** | **B1ME58** (GyrA, 839 aa) and **B1ME45** (GyrB, 675 aa) |
| **Template** | PDB **5BS8** - *M. tuberculosis* DNA gyrase-DNA complex (X-ray, 2.8 Å) |
| **Complex** | Heterotetramer (GyrA₂-GyrB₂) with DNA and Mg²⁺ ions |
| **Method** | Rosetta partial_thread + FastRelax |

## Modeling Results

### Sequence Identity
- **GyrA alignment**: 90.3% identity (439/486 aligned residues)
  - Coverage: residues 17-503 of target
- **GyrB alignment**: 95.5% identity (234/245 aligned residues)
  - Coverage: residues 425-675 of target

### Rosetta Scores (ref2015)
- **GyrA relaxed model**: -1551.42 REU
- **GyrB relaxed model**: -717.83 REU

### Output Files

| File | Description |
|------|-------------|
| `output/mabs_gyrA_new_relaxed_*.pdb` | Relaxed GyrA model (486 residues) |
| `output/mabs_gyrB_new_relaxed_*.pdb` | Relaxed GyrB model (245 residues) |
| `output/mabs_gyrase_tetramer_dna_mg.pdb` | Full tetramer with DNA and Mg²⁺ |
| `output/mabs_gyrase_tetramer_protein_only.pdb` | Protein-only tetramer |

### Tetramer Chain Assignment
- **Chain A**: GyrA subunit 1
- **Chain B**: GyrB subunit 1  
- **Chain C**: GyrA subunit 2 (symmetry copy)
- **Chain D**: GyrB subunit 2 (symmetry copy)
- **Chains E-H**: DNA from template
- **Mg²⁺ ions**: 4 ions from template

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
