#!/usr/bin/env python3
"""
ESMFold Structure Prediction for M. abscessus DNA Gyrase
Predicts GyrA, GyrB monomers and GyrA-GyrB heterodimer
"""

import os
import sys
import torch
import esm
from pathlib import Path
from Bio import SeqIO
import gc

# Setup paths
PROJECT_DIR = Path(__file__).parent.parent.parent
SEQ_DIR = PROJECT_DIR / "input" / "sequences"
OUTPUT_DIR = PROJECT_DIR / "output" / "ai_predictions" / "esmfold"

def load_sequence(fasta_file):
    """Load sequence from FASTA file"""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq), record.id

def predict_structure(model, sequence, output_path, name="structure"):
    """Run ESMFold prediction and save PDB"""
    print(f"\nPredicting structure for {name}")
    print(f"  Sequence length: {len(sequence)} residues")
    
    # Check GPU memory
    if torch.cuda.is_available():
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        print(f"  GPU memory: {gpu_mem:.1f} GB")
        
        # Estimate memory needed (~0.1 GB per 100 residues for ESMFold)
        est_mem = len(sequence) * 0.0015  # Conservative estimate
        print(f"  Estimated memory needed: {est_mem:.1f} GB")
        
        if est_mem > gpu_mem * 0.9:
            print(f"  WARNING: May run out of GPU memory!")
    
    with torch.no_grad():
        try:
            output = model.infer_pdb(sequence)
            
            # Save PDB
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(output)
            print(f"  Saved: {output_path}")
            
            # Clear GPU memory
            torch.cuda.empty_cache()
            gc.collect()
            
            return True
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                print(f"  ERROR: Out of GPU memory for {name}")
                torch.cuda.empty_cache()
                gc.collect()
                return False
            raise

def main():
    print("=" * 70)
    print("ESMFold Structure Prediction")
    print("M. abscessus DNA Gyrase (GyrA & GyrB)")
    print("=" * 70)
    
    # Check GPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nUsing device: {device}")
    if device == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    # Load ESMFold model
    print("\nLoading ESMFold model...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    
    if device == "cuda":
        model = model.cuda()
        # Use half precision to save memory
        model.set_chunk_size(128)  # Trade speed for memory
    
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
    success = predict_structure(
        model, gyrA_seq, 
        OUTPUT_DIR / "GyrA" / "gyrA_esmfold.pdb",
        "GyrA monomer"
    )
    results["GyrA"] = "Success" if success else "Failed"
    
    # 2. Predict GyrB monomer
    print("\n" + "-" * 70)
    print("2. GyrB Monomer Prediction")
    print("-" * 70)
    success = predict_structure(
        model, gyrB_seq,
        OUTPUT_DIR / "GyrB" / "gyrB_esmfold.pdb", 
        "GyrB monomer"
    )
    results["GyrB"] = "Success" if success else "Failed"
    
    # 3. Predict GyrA-GyrB heterodimer (chains separated by ":")
    print("\n" + "-" * 70)
    print("3. GyrA-GyrB Heterodimer Prediction")
    print("-" * 70)
    heterodimer_seq = f"{gyrA_seq}:{gyrB_seq}"
    print(f"  Total residues: {len(gyrA_seq) + len(gyrB_seq)}")
    
    # This may fail due to memory - it's a large complex
    success = predict_structure(
        model, heterodimer_seq,
        OUTPUT_DIR / "complex" / "gyrAB_heterodimer_esmfold.pdb",
        "GyrA-GyrB heterodimer"
    )
    results["Heterodimer"] = "Success" if success else "Failed (OOM expected)"
    
    # 4. Optionally try A2B2 tetramer (likely to fail on most GPUs)
    print("\n" + "-" * 70)
    print("4. A2B2 Tetramer Prediction (experimental)")
    print("-" * 70)
    tetramer_seq = f"{gyrA_seq}:{gyrA_seq}:{gyrB_seq}:{gyrB_seq}"
    total_res = 2 * len(gyrA_seq) + 2 * len(gyrB_seq)
    print(f"  Total residues: {total_res}")
    print(f"  NOTE: This is very large and will likely run out of memory")
    
    # Only attempt if we have > 40GB GPU memory
    if torch.cuda.is_available():
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        if gpu_mem >= 40:
            success = predict_structure(
                model, tetramer_seq,
                OUTPUT_DIR / "complex" / "gyrA2B2_tetramer_esmfold.pdb",
                "A2B2 Tetramer"
            )
            results["Tetramer"] = "Success" if success else "Failed (OOM)"
        else:
            print(f"  Skipping - need ~40GB GPU memory, have {gpu_mem:.1f} GB")
            results["Tetramer"] = "Skipped (insufficient GPU memory)"
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, status in results.items():
        print(f"  {name}: {status}")
    
    print("\nOutput directory:", OUTPUT_DIR)
    print("=" * 70)

if __name__ == "__main__":
    main()
