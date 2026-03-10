# Schrödinger Setup for M. abscessus Gyrase Project

## Quick Start

```bash
# 1. Source the environment
source setup_schrodinger.sh

# 2. Launch Maestro GUI
$SCHRODINGER/maestro
```

## Available Workflow Scripts

### 1. Pharmacophore Modeling (`pharmacophore_workflow.sh`)
```bash
# Prepare ligands from SMILES
./pharmacophore_workflow.sh prepare ligands.smi

# Create Phase database
./pharmacophore_workflow.sh create_db my_database.phdb structures.mae

# Find common pharmacophore from actives
./pharmacophore_workflow.sh find_common actives.mae

# Screen database with hypothesis
./pharmacophore_workflow.sh screen database.phdb hypothesis.phypo results

# Build 3D-QSAR model
./pharmacophore_workflow.sh qsar training.mae hypothesis.phypo activity_prop
```

### 2. Glide Docking (`glide_docking_workflow.sh`)
```bash
# Prepare protein structure
./glide_docking_workflow.sh prepare_protein receptor.pdb

# Generate receptor grid
./glide_docking_workflow.sh generate_grid receptor_prep.mae

# Prepare ligands
./glide_docking_workflow.sh prepare_ligands ligands.smi

# Run SP docking
./glide_docking_workflow.sh dock_sp grid.zip ligands_prep.mae

# Run XP docking for rescoring
./glide_docking_workflow.sh dock_xp grid.zip sp_hits.mae
```

## Key Modules

| Category | Module | Command |
|----------|--------|---------|
| GUI | Maestro | `$SCHRODINGER/maestro` |
| Protein Prep | PrepWizard | `$SCHRODINGER/utilities/prepwizard` |
| Ligand Prep | LigPrep | `$SCHRODINGER/ligprep` |
| pKa/Tautomers | Epik | `$SCHRODINGER/epik` |
| Conformers | ConfGen | `$SCHRODINGER/confgen` |
| Docking | Glide | `$SCHRODINGER/glide` |
| Induced Fit | IFD | `$SCHRODINGER/ifd` |
| Binding Sites | SiteMap | `$SCHRODINGER/sitemap` |
| FEP | FEP+ | `$SCHRODINGER/fep_plus` |

### Phase (Pharmacophore) Tools
- `phase_database` - Create/manage 3D databases
- `phase_screen` - Screen with pharmacophore hypothesis
- `phase_find_common` - Find common pharmacophore features
- `phase_hypo_refine` - Refine pharmacophore hypothesis
- `phase_build_qsar` - Build 3D-QSAR model
- `phase_qsar` - Apply QSAR model

## For the Gyrase Project

### Suggested Workflow
1. Prepare protein model: `prepwizard gyrase_model.pdb gyrase_prep.mae`
2. Identify binding site: `sitemap -prot gyrase_prep.mae -j sitemap_job`
3. Generate grid centered on QRDR region
4. Prepare fluoroquinolones with LigPrep
5. Dock with Glide SP/XP
6. Build pharmacophore from docked poses
7. Screen compound library

### Integration with Existing Docking Data
Your project already has AutoDock Vina docking results in `docking/`. You can:
1. Import docked poses into Maestro
2. Rescore with Glide
3. Generate pharmacophore from best poses
4. Use Phase for virtual screening

## Files
- `setup_schrodinger.sh` - Environment setup script
- `pharmacophore_workflow.sh` - Pharmacophore modeling workflow
- `glide_docking_workflow.sh` - Glide docking workflow
- `MODULE_STATUS.md` - Detailed module availability report
