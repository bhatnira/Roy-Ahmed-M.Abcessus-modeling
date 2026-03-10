# Google Colab Notebooks for M. abscessus Gyrase Modeling

This folder contains Google Colab-compatible notebooks for homology modeling and drug-target binding prediction.

## 📚 Notebooks

### Part 1: Homology Modeling
**[01_Homology_Modeling.ipynb](01_Homology_Modeling.ipynb)**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/YOUR_REPO/blob/main/colab_notebooks/01_Homology_Modeling.ipynb)

Contents:
- Target sequence analysis
- Template-based sequence alignment (vs M. tuberculosis 5BS8)
- ESMFold structure prediction
- Structure quality assessment (pLDDT confidence scores)
- RosettaCM reference for local modeling

### Part 2: DrugCLIP Analysis
**[02_DrugCLIP_Analysis.ipynb](02_DrugCLIP_Analysis.ipynb)**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/YOUR_REPO/blob/main/colab_notebooks/02_DrugCLIP_Analysis.ipynb)

Contents:
- DrugCLIP contrastive learning framework
- Drug-target binding predictions
- Fluoroquinolone analysis (Ciprofloxacin, Levofloxacin, Moxifloxacin)
- QRDR (binding site) analysis
- Molecular docking validation

## 🚀 Quick Start

### Option 1: Open in Google Colab
1. Click the "Open in Colab" badge for the notebook you want
2. Connect to a GPU runtime (Runtime → Change runtime type → GPU)
3. Run Part 1 first, then Part 2

### Option 2: Local Jupyter
```bash
# Install dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook
```

## 📋 Requirements

```
biopython>=1.79
numpy>=1.21.0
pandas>=1.3.0
matplotlib>=3.5.0
seaborn>=0.11.0
torch>=1.10.0
transformers>=4.20.0
rdkit-pypi>=2022.03.1
py3Dmol>=2.0.0
fair-esm>=2.0.0
vina>=1.2.3
meeko>=0.4.0
prody>=2.0
```

## 📁 Data Files

All required data files are included in the `data/` folder:

```
data/
├── GyrA_B1ME58.fasta          # GyrA target sequence (839 aa)
├── GyrB_B1ME45.fasta          # GyrB target sequence (675 aa)
├── drugs.csv                   # Drug SMILES database
├── templates/                  # Template PDB structures
│   ├── 5bs8_chainA_fixed.pdb  # M. tuberculosis GyrA
│   ├── 5bs8_chainB_fixed.pdb  # M. tuberculosis GyrB
│   └── 5bs8_protein_only.pdb  # Full complex
├── alignments/                 # Sequence alignments
│   ├── alignment_gyrA_proper.grishin
│   └── alignment_gyrB_proper.grishin
├── ligands/                    # Ligand structures
│   ├── ciprofloxacin.sdf
│   ├── levofloxacin.sdf
│   └── moxifloxacin.sdf
└── models/                     # Pre-built relaxed models
    ├── mabs_gyrA_relaxed.pdb
    └── mabs_gyrB_relaxed.pdb
```

## 🧬 Target Proteins

| Subunit | UniProt ID | Length | Function | Template Identity |
|---------|------------|--------|----------|-------------------|
| **GyrA** | [B1ME58](https://www.uniprot.org/uniprot/B1ME58) | 839 aa | DNA breakage-reunion domain | 90.3% to M.tb |
| **GyrB** | [B1ME45](https://www.uniprot.org/uniprot/B1ME45) | 675 aa | ATPase/TOPRIM domain | 95.5% to M.tb |

### Template Structure
- **PDB**: 5BS8 - *M. tuberculosis* DNA Gyrase-DNA-Fluoroquinolone complex
- **Resolution**: 2.8 Å
- **Method**: X-ray Crystallography

## 💊 DrugCLIP: Drug-Target Binding Prediction

### What is DrugCLIP?

DrugCLIP is a contrastive learning framework that predicts drug-target interactions by learning joint representations of molecules and proteins.

### How it Works:

```
Drug (SMILES) → Drug Encoder → Drug Embedding ─┐
                                                ├→ Cosine Similarity → Binding Score
Protein (Seq) → Protein Encoder → Protein Embedding ─┘
```

### Components:
- **Drug Encoder**: ChemBERTa (transformer trained on SMILES)
- **Protein Encoder**: ESM-2 (transformer trained on protein sequences)
- **Contrastive Loss**: Learns to maximize similarity for binding pairs

### Drugs Tested:

| Drug | Mechanism | Expected Binding |
|------|-----------|------------------|
| Ciprofloxacin | DNA gyrase inhibitor | ✅ Yes |
| Levofloxacin | DNA gyrase inhibitor | ✅ Yes |
| Moxifloxacin | DNA gyrase inhibitor | ✅ Yes |
| Rifampicin | RNA polymerase inhibitor | ❌ No (control) |

## 🔬 Homology Modeling Pipeline

### ESMFold (Cloud-Compatible)

```python
# ESMFold API usage
response = requests.post(
    "https://api.esmatlas.com/foldSequence/v1/pdb/",
    data=sequence,
    headers={"Content-Type": "text/plain"}
)
```

**Advantages:**
- No local GPU required
- Fast predictions (seconds)
- No MSA computation needed

**Limitations:**
- Sequence length limit (~400 aa for free tier)
- Requires internet connection

### Alternative: RosettaCM (Local)

For higher accuracy with full-length sequences, see the main project's RosettaCM pipeline in `/scripts/`.

## 📊 Output Files

### Part 1 - Homology Modeling:
```
output/
├── esmfold/           # ESMFold predicted structures
│   ├── GyrA_esmfold.pdb
│   └── GyrB_esmfold.pdb
├── analysis/          # Quality assessment plots
│   ├── aa_composition.png
│   ├── alignment_identity.png
│   └── plddt_distribution.png
└── alignments/        # Generated Grishin alignments
    ├── alignment_gyrA.grishin
    └── alignment_gyrB.grishin
```

### Part 2 - DrugCLIP Analysis:
```
output/
├── drugclip/          # Binding predictions
│   ├── binding_predictions.csv
│   ├── drugclip_results.png
│   ├── qrdr_binding.png
│   └── final_summary.csv
└── docking/           # Prepared ligand files
    ├── Ciprofloxacin.sdf
    ├── Levofloxacin.sdf
    └── Moxifloxacin.sdf
```

## 📖 Methods Overview

### 1. Sequence Analysis
- Amino acid composition
- Hydrophobicity analysis
- Functional domain identification

### 2. Structure Prediction (ESMFold)
- Single-sequence protein folding
- pLDDT confidence scores
- 3D visualization with py3Dmol

### 3. DrugCLIP Binding Prediction
- ChemBERTa drug encoding
- ESM-2 protein encoding
- Contrastive similarity scoring

### 4. Molecular Docking (AutoDock Vina)
- QRDR region targeting
- Binding affinity estimation
- Pose visualization

## 🔗 References

1. **ESMFold**: Lin Z, et al. (2023) Evolutionary-scale prediction of atomic-level protein structure. *Science* 379:1123-1130

2. **DrugCLIP**: Gao Z, et al. (2024) DrugCLIP: Contrastive learning for drug-target interaction prediction. *Nature Communications* 15:2043

3. **ChemBERTa**: Chithrananda S, et al. (2020) ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction. *arXiv*:2010.09885

4. **ESM-2**: Lin Z, et al. (2022) Language models of protein sequences at the scale of evolution. *bioRxiv*

5. **AutoDock Vina**: Trott O, Olson AJ (2010) AutoDock Vina: improving the speed and accuracy of docking. *J Comput Chem* 31:455-461

## 📧 Contact

For questions about this analysis, please open an issue in the main repository.

---

*Last updated: February 2026*
