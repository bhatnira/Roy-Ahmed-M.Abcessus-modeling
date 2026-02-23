# Detailed Methods

## 1. Sequence Retrieval and Preparation

### Target Sequences
- **GyrA**: UniProt B1ME58 (*Mycobacterium abscessus* ATCC 19977)
  - Full length: 839 amino acids
  - Modeled region: residues 17-503 (DNA breakage-reunion domain)
  
- **GyrB**: UniProt B1ME45 (*Mycobacterium abscessus* ATCC 19977)
  - Full length: 675 amino acids
  - Modeled region: residues 425-675 (TOPRIM domain)

### Template Structure
- **PDB ID**: 5BS8
- **Organism**: *Mycobacterium tuberculosis*
- **Resolution**: 2.8 Å
- **Method**: X-ray crystallography
- **Citation**: Blower TR et al. (2016) Crystal structure and stability of gyrase–fluoroquinolone cleaved complexes from Mycobacterium tuberculosis. PNAS 113(7):1706-13.

## 2. Sequence Alignment

### Method
Pairwise sequence alignment using Needleman-Wunsch algorithm (global alignment).

### Parameters
- Substitution matrix: BLOSUM62
- Gap open penalty: -10
- Gap extend penalty: -0.5

### Output Format
Grishin alignment format required by RosettaCM:
```
## target_name template.pdb
#
scores_from_program: 0
0 TARGET_SEQUENCE...
0 TEMPLATE_SEQUENCE...
```

### Alignment Statistics
| Subunit | Identity | Similarity | Gaps | Coverage |
|---------|----------|------------|------|----------|
| GyrA | 90.3% | 94.2% | 2.1% | 486/839 |
| GyrB | 95.5% | 97.8% | 0.8% | 245/675 |

## 3. Template Preparation

### Steps
1. Download PDB 5BS8 from RCSB
2. Extract individual chains:
   - Chain A: GyrA (template for target GyrA)
   - Chain B: GyrB (template for target GyrB)
   - Chains E-H: DNA
   - MG: Magnesium ions
3. Fix missing atoms using REDUCE
4. Renumber residues to match alignment
5. Clean HETATM records

### Files Generated
- `5bs8_chainA_fixed.pdb`
- `5bs8_chainB_fixed.pdb`
- `5bs8_dna_mg.pdb`
- `5bs8_protein_only.pdb`

## 4. RosettaCM Modeling Protocol

### Software
- Rosetta 2021.16 (weekly release)
- RosettaScripts XML protocol

### Hybridize Protocol Parameters
```xml
<Hybridize name="hybridize"
    stage1_scorefxn="score3"
    stage2_scorefxn="score4_smooth"
    fa_scorefxn="ref2015"
    stage1_increase_cycles="1.0"
    stage2_increase_cycles="1.0"
    dualspace="true">
```

### Three-Stage Protocol
1. **Stage 1 - Centroid Assembly**
   - Score function: score3
   - Fragment insertion for gaps
   - Monte Carlo optimization
   
2. **Stage 2 - Full-Atom Conversion**
   - Score function: score4_smooth
   - Rotamer packing (Dunbrack library)
   - Side-chain optimization
   
3. **Stage 3 - Cartesian Refinement**
   - Score function: ref2015
   - All-atom minimization
   - Backbone and side-chain relaxation

### Command Line
```bash
rosetta_scripts.default.linuxgccrelease \
    -parser:protocol rosettacm_gyrase.xml \
    -in:file:fasta target.fasta \
    -nstruct 10 \
    -out:prefix mabs_
```

## 5. FastRelax Refinement

### Protocol
Iterative all-atom refinement with fa_rep ramping.

### Cycles
1. fa_rep = 0.02 → Pack + Minimize
2. fa_rep = 0.25 → Pack + Minimize
3. fa_rep = 0.55 → Pack + Minimize
4. fa_rep = 1.0 → Final minimization

### MoveMap
- Backbone: movable for all residues
- Side chains: movable for all residues
- Jumps: fixed (preserve relative chain positions)

### Command Line
```bash
rosetta_scripts.default.linuxgccrelease \
    -parser:protocol fastrelax.xml \
    -s input.pdb \
    -nstruct 5 \
    -relax:default_repeats 5
```

## 6. Tetramer Assembly

### Method
1. Thread GyrA model onto chains A and C positions from template
2. Thread GyrB model onto chains B and D positions
3. Merge coordinates into single PDB file
4. Add DNA coordinates from template (conserved)
5. Add Mg²⁺ ions at catalytic sites

### Script
`scripts/assemble_tetramer.py`

## 7. Molecular Docking

### Software
- AutoDock Vina 1.2.3
- Open Babel 3.1.1 (ligand preparation)
- MGLTools 1.5.7 (receptor preparation)

### Ligand Preparation
1. Download 2D structures from PubChem
2. Generate 3D conformers
3. Add hydrogens at pH 7.4
4. Minimize with MMFF94 force field
5. Convert to PDBQT format

### Receptor Preparation
1. Extract binding site region
2. Remove water molecules
3. Add polar hydrogens
4. Assign Gasteiger charges
5. Convert to PDBQT format

### Docking Parameters
```
center_x = 12.5
center_y = -8.3
center_z = 45.2
size_x = 25
size_y = 25
size_z = 25
exhaustiveness = 32
num_modes = 9
energy_range = 3
```

### Binding Site Definition
- Centered on QRDR region (residues 80-95)
- Includes Mg²⁺ coordination site
- Encompasses DNA intercalation pocket

## 8. Model Validation

### Metrics Evaluated
1. **Ramachandran analysis**: φ/ψ backbone angles
2. **Rotamer analysis**: Side-chain conformations
3. **Clashscore**: Steric clashes per 1000 atoms
4. **RMSD**: Comparison to template structure
5. **MolProbity score**: Combined quality metric

### Software
- Rosetta score functions
- MolProbity server (backup validation)

## References

1. Song Y et al. (2013) High-resolution comparative modeling with RosettaCM. Structure 21:1735-1742.
2. Alford RF et al. (2017) The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. J Chem Theory Comput 13:3031-3048.
3. Trott O, Olson AJ (2010) AutoDock Vina: improving the speed and accuracy of docking with a new scoring function. J Comput Chem 31:455-461.
4. Blower TR et al. (2016) Crystal structure and stability of gyrase–fluoroquinolone cleaved complexes. PNAS 113:1706-1713.
