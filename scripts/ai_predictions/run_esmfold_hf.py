#!/usr/bin/env python3
"""
ESMFold Structure Prediction using Hugging Face Transformers
M. abscessus DNA Gyrase (GyrA & GyrB)
"""

import os
import sys
import torch
import gc
from pathlib import Path
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

# Setup paths
PROJECT_DIR = Path(__file__).parent.parent.parent
SEQ_DIR = PROJECT_DIR / "input" / "sequences"
OUTPUT_DIR = PROJECT_DIR / "output" / "ai_predictions" / "esmfold"

def load_sequence(fasta_file):
    """Load sequence from FASTA file"""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq), record.id

def convert_outputs_to_pdb(outputs):
    """Convert ESMFold outputs to PDB format"""
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs

def predict_structure(model, tokenizer, sequence, output_path, name="structure"):
    """Run ESMFold prediction and save PDB"""
    print(f"\nPredicting structure for {name}")
    print(f"  Sequence length: {len(sequence)} residues")
    
    # Check GPU memory
    if torch.cuda.is_available():
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        used_mem = torch.cuda.memory_allocated() / 1e9
        print(f"  GPU memory: {gpu_mem:.1f} GB total, {used_mem:.1f} GB used")
    
    try:
        # Tokenize
        inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
        
        if torch.cuda.is_available():
            inputs = {k: v.cuda() for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
        
        # Convert to PDB
        pdbs = convert_outputs_to_pdb(outputs)
        
        # Save PDB
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(pdbs[0])
        print(f"  Saved: {output_path}")
        
        # Get pLDDT score
        plddt = outputs["plddt"].mean().item()
        print(f"  Mean pLDDT: {plddt:.1f}")
        
        # Clear GPU memory
        del outputs, inputs
        torch.cuda.empty_cache()
        gc.collect()
        
        return True, plddt
        
    except RuntimeError as e:
        if "out of memory" in str(e).lower():
            print(f"  ERROR: Out of GPU memory for {name}")
            torch.cuda.empty_cache()
            gc.collect()
            return False, 0
        raise

def main():
    print("=" * 70)
    print("ESMFold Structure Prediction (Hugging Face)")
    print("M. abscessus DNA Gyrase (GyrA & GyrB)")
    print("=" * 70)
    
    # Check GPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nUsing device: {device}")
    if device == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    # Load ESMFold model from Hugging Face
    print("\nLoading ESMFold model from Hugging Face...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    
    if device == "cuda":
        model = model.cuda()
    model.eval()
    
    # Enable memory efficient attention
    model.esm = model.esm.half()  # Use FP16 for ESM trunk
    model.trunk.set_chunk_size(64)  # Trade speed for memory
    
    # Load sequences
    print("\nLoading sequences...")
    gyrA_seq, gyrA_id = load_sequence(SEQ_DIR / "B1ME58_GyrA.fasta")
    gyrB_seq, gyrB_id = load_sequence(SEQ_DIR / "B1ME45_GyrB.fasta")
    
    print(f"  GyrA: {len(gyrA_seq)} residues")
    print(f"  GyrB: {len(gyrB_seq)} residues")
    
    results = {}
    
    # 1. Predict GyrA monomer
    print("\n" + "-" * 70)
    print("1. GyrA Monomer Prediction")
    print("-" * 70)
    success, plddt = predict_structure(
        model, tokenizer, gyrA_seq, 
        OUTPUT_DIR / "GyrA" / "gyrA_esmfold.pdb",
        "GyrA monomer"
    )
    results["GyrA"] = {"status": "Success" if success else "Failed", "plddt": plddt}
    
    # 2. Predict GyrB monomer
    print("\n" + "-" * 70)
    print("2. GyrB Monomer Prediction")
    print("-" * 70)
    success, plddt = predict_structure(
        model, tokenizer, gyrB_seq,
        OUTPUT_DIR / "GyrB" / "gyrB_esmfold.pdb", 
        "GyrB monomer"
    )
    results["GyrB"] = {"status": "Success" if success else "Failed", "plddt": plddt}
    
    # 3. Predict GyrA-GyrB heterodimer
    # For multimer, we need to process differently or use ESM-MSA
    print("\n" + "-" * 70)
    print("3. GyrA-GyrB Heterodimer Prediction")
    print("-" * 70)
    heterodimer_seq = gyrA_seq + gyrB_seq  # Concatenate for heterodimer
    print(f"  Total residues: {len(heterodimer_seq)}")
    print(f"  NOTE: Large complex, may run out of memory")
    
    success, plddt = predict_structure(
        model, tokenizer, heterodimer_seq,
        OUTPUT_DIR / "complex" / "gyrAB_heterodimer_esmfold.pdb",
        "GyrA-GyrB heterodimer"
    )
    results["Heterodimer"] = {"status": "Success" if success else "Failed (OOM)", "plddt": plddt}
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Structure':<20} {'Status':<15} {'pLDDT':<10}")
    print("-" * 45)
    for name, data in results.items():
        plddt_str = f"{data['plddt']:.1f}" if data['plddt'] > 0 else "N/A"
        print(f"{name:<20} {data['status']:<15} {plddt_str:<10}")
    
    print("\nOutput directory:", OUTPUT_DIR)
    print("=" * 70)
    
    # Save summary
    summary_path = OUTPUT_DIR / "prediction_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("ESMFold Prediction Summary\n")
        f.write("=" * 50 + "\n\n")
        for name, data in results.items():
            f.write(f"{name}: {data['status']}, pLDDT={data['plddt']:.1f}\n")
    print(f"Summary saved to: {summary_path}")

if __name__ == "__main__":
    main()
