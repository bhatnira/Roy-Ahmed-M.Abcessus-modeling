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
| `output/threaded/` | Threaded models before relaxation |
| `output/relaxed/` | FastRelax-refined individual chains |
| `output/tetramer/` | Assembled heterotetramer models |
| `output/scores/` | Rosetta score files |

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

## Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| **Rosetta** | 2021.16+ (weekly release) | Comparative modeling, FastRelax |
| **Python** | 3.8+ | Scripts and analysis |
| **BioPython** | 1.79+ | Sequence/structure parsing |
| **AutoDock Vina** | 1.2.3 | Molecular docking |
| **PyMOL** | 2.5+ | Visualization |
| **Open Babel** | 3.1.1 | Ligand preparation |

### Installation

```bash
# Install Python dependencies
pip install -r requirements.txt

# Rosetta installation (requires license)
# See: https://www.rosettacommons.org/software/license-and-download
```

## Model Validation Metrics

### Structural Quality (from MolProbity/Rosetta analysis)
| Metric | GyrA | GyrB | Tetramer |
|--------|------|------|----------|
| Ramachandran favored | 96.2% | 95.8% | 95.4% |
| Ramachandran outliers | 0.4% | 0.6% | 0.5% |
| Rotamer outliers | 1.2% | 1.8% | 1.5% |
| Clashscore | 3.2 | 4.1 | 4.8 |
| Cα RMSD to template | 0.82 Å | 0.75 Å | 1.22 Å |

### Rosetta Energy Breakdown (REU)
| Term | GyrA | GyrB |
|------|------|------|
| total_score | -1551.42 | -717.83 |
| fa_atr | -2845.6 | -1421.2 |
| fa_rep | +421.3 | +198.7 |
| fa_sol | +721.4 | +398.2 |
| hbond_total | -198.5 | -112.4 |

## Molecular Docking Results

### Fluoroquinolone Binding Affinities (AutoDock Vina)

| Ligand | *M. abscessus* | *M. tuberculosis* | ΔΔG |
|--------|----------------|-------------------|-----|
| **Moxifloxacin** | -7.91 kcal/mol | -7.24 kcal/mol | **-0.67** |
| Levofloxacin | -7.65 kcal/mol | -7.31 kcal/mol | -0.34 |
| Ciprofloxacin | -7.36 kcal/mol | -7.27 kcal/mol | -0.09 |

**Key Finding**: All fluoroquinolones show stronger predicted binding to *M. abscessus* DNA gyrase compared to *M. tuberculosis*, with moxifloxacin showing the largest improvement.

## Directory Structure

```
rosetta-cm/
├── input/
│   ├── sequences/            # Target FASTA sequences
│   │   ├── B1ME58_GyrA.fasta
│   │   └── B1ME45_GyrB.fasta
│   ├── alignments_final/     # Final Grishin alignments
│   │   ├── alignment_gyrA_proper.grishin
│   │   └── alignment_gyrB_proper.grishin
│   └── templates/            # Template PDB files
│       └── 5bs8.pdb
├── output/
│   ├── threaded/             # Threaded models
│   ├── relaxed/              # FastRelax-refined chains
│   ├── tetramer/             # Assembled complexes
│   │   ├── mabs_gyrase_tetramer_dna_mg.pdb
│   │   └── mabs_gyrase_tetramer_protein_only.pdb
│   └── scores/               # Rosetta score files
├── scripts/
│   └── assemble_tetramer.py  # Tetramer assembly script
├── xml/
│   └── fastrelax.xml         # Rosetta XML protocols
├── logs/                     # Rosetta logs and crash files
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
