# Schrödinger Suite 2026-1 - Module Status Report
Generated: March 10, 2026

## Installation Path
`/opt/schrodinger/schrodinger2026-1`

## Module Status (Verified Working)

### Ligand Preparation
| Module | Command | Status |
|--------|---------|--------|
| LigPrep | `$SCHRODINGER/ligprep` | ✅ Working |
| Epik | `$SCHRODINGER/epik` | ✅ Working |
| ConfGen | `$SCHRODINGER/confgen` | ✅ Working |

### Docking & Scoring
| Module | Command | Status |
|--------|---------|--------|
| Glide | `$SCHRODINGER/glide` | ✅ Working |
| IFD (Induced Fit) | `$SCHRODINGER/ifd` | ✅ Working |
| CovDock | `$SCHRODINGER/covalent_docking` | ✅ Available |

### Pharmacophore Modeling (Phase)
| Module | Command | Status |
|--------|---------|--------|
| Phase Database | `$SCHRODINGER/phase_database` | ✅ Working |
| Phase Screen | `$SCHRODINGER/phase_screen` | ✅ Working |
| Phase Find Common | `$SCHRODINGER/phase_find_common` | ✅ Working |
| Phase Hypo Refine | `$SCHRODINGER/phase_hypo_refine` | ✅ Working |
| Phase Build QSAR | `$SCHRODINGER/phase_build_qsar` | ✅ Available |
| Phase QSAR | `$SCHRODINGER/phase_qsar` | ✅ Available |
| Field QSAR | `$SCHRODINGER/phase_fqsar` | ✅ Available |

### Structure Analysis
| Module | Command | Status |
|--------|---------|--------|
| Protein Prep Wizard | `$SCHRODINGER/utilities/prepwizard` | ✅ Working |
| SiteMap | `$SCHRODINGER/sitemap` | ✅ Working |
| BioLuminate | `$SCHRODINGER/bioluminate` | ✅ Available |

### Free Energy Calculations
| Module | Command | Status |
|--------|---------|--------|
| FEP+ | `$SCHRODINGER/fep_plus` | ✅ Working |
| FEP Absolute Binding | `$SCHRODINGER/fep_absolute_binding` | ✅ Available |
| FEP Residue Scanning | `$SCHRODINGER/fep_residue_scanning` | ✅ Available |
| FEP Solubility | `$SCHRODINGER/fep_solubility` | ✅ Available |

### Molecular Dynamics
| Module | Command | Status |
|--------|---------|--------|
| Desmond | `$SCHRODINGER/desmond` | ✅ Available |

### Cheminformatics
| Module | Command | Status |
|--------|---------|--------|
| Canvas | `$SCHRODINGER/canvas-v6.7/` | ✅ Available |

### GUI
| Module | Command | Status |
|--------|---------|--------|
| Maestro | `$SCHRODINGER/maestro` | ✅ Working |

## Quick Start Commands

```bash
# Source environment
source setup_schrodinger.sh

# Launch Maestro GUI
$SCHRODINGER/maestro

# Prepare protein structure
$SCHRODINGER/utilities/prepwizard input.pdb output.mae

# Prepare ligands
$SCHRODINGER/ligprep -ismi ligands.smi -omae ligands_prep.mae

# Grid generation (requires Maestro GUI or grid file)
$SCHRODINGER/glide grid.in

# Docking
$SCHRODINGER/glide dock.in

# Pharmacophore screening
$SCHRODINGER/phase_screen ligands.mae hypothesis.phypo jobname
```

## Total Modules Available: 97
