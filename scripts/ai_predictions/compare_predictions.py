#!/usr/bin/env python3
"""
Compare AI Structure Predictions
Analyzes models from ESMFold, ColabFold, AlphaFold2, RoseTTAFold, and RosettaCM
"""

import os
import sys
import glob
import json
from pathlib import Path

# Try to import analysis libraries
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("Warning: numpy not installed, some analyses will be skipped")

PROJECT_DIR = Path(__file__).parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "output" / "ai_predictions"
ROSETTACM_DIR = PROJECT_DIR / "output" / "relaxed"

def parse_pdb_bfactor(pdb_file):
    """Extract B-factors (often pLDDT scores) from PDB file"""
    bfactors = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                try:
                    bfactor = float(line[60:66].strip())
                    bfactors.append(bfactor)
                except:
                    pass
    return bfactors

def get_pdb_files(directory, pattern="*.pdb"):
    """Find all PDB files in a directory"""
    return list(Path(directory).glob(f"**/{pattern}"))

def calculate_stats(values):
    """Calculate basic statistics"""
    if not HAS_NUMPY or not values:
        return {"mean": 0, "std": 0, "min": 0, "max": 0, "n": 0}
    arr = np.array(values)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "n": len(arr)
    }

def analyze_method(method_name, method_dir):
    """Analyze predictions from a single method"""
    results = {"method": method_name, "proteins": {}}
    
    pdbs = get_pdb_files(method_dir)
    
    for pdb_file in pdbs:
        # Determine protein name from path or filename
        if "GyrA" in str(pdb_file) or "gyrA" in str(pdb_file):
            protein = "GyrA"
        elif "GyrB" in str(pdb_file) or "gyrB" in str(pdb_file):
            protein = "GyrB"
        else:
            protein = pdb_file.stem
        
        # Get confidence scores (pLDDT stored in B-factor column)
        bfactors = parse_pdb_bfactor(pdb_file)
        
        results["proteins"][protein] = {
            "pdb_file": str(pdb_file),
            "n_residues": len(bfactors),
            "confidence": calculate_stats(bfactors)
        }
    
    return results

def main():
    print("=" * 70)
    print("AI Structure Prediction Comparison")
    print("M. abscessus DNA Gyrase (GyrA & GyrB)")
    print("=" * 70)
    print()
    
    # Define methods and their directories
    methods = {
        "RosettaCM": ROSETTACM_DIR,
        "ESMFold": OUTPUT_DIR / "esmfold",
        "ColabFold": OUTPUT_DIR / "colabfold",
        "AlphaFold2": OUTPUT_DIR / "alphafold2",
        "RoseTTAFold": OUTPUT_DIR / "rosettafold"
    }
    
    all_results = {}
    
    for method_name, method_dir in methods.items():
        print(f"\n{'-' * 70}")
        print(f"Analyzing: {method_name}")
        print(f"Directory: {method_dir}")
        print(f"{'-' * 70}")
        
        if not method_dir.exists():
            print(f"  Directory not found, skipping...")
            continue
            
        pdbs = get_pdb_files(method_dir)
        if not pdbs:
            print(f"  No PDB files found, skipping...")
            continue
        
        print(f"  Found {len(pdbs)} PDB file(s)")
        
        results = analyze_method(method_name, method_dir)
        all_results[method_name] = results
        
        for protein, data in results["proteins"].items():
            conf = data["confidence"]
            print(f"\n  {protein}:")
            print(f"    File: {Path(data['pdb_file']).name}")
            print(f"    Residues: {data['n_residues']}")
            if conf["n"] > 0:
                print(f"    Confidence (pLDDT/B-factor):")
                print(f"      Mean: {conf['mean']:.1f}")
                print(f"      Std:  {conf['std']:.1f}")
                print(f"      Range: {conf['min']:.1f} - {conf['max']:.1f}")
    
    # Summary comparison table
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON")
    print("=" * 70)
    print()
    print(f"{'Method':<15} {'GyrA pLDDT':<15} {'GyrB pLDDT':<15} {'Files':<10}")
    print("-" * 55)
    
    for method_name, results in all_results.items():
        proteins = results.get("proteins", {})
        gyrA_conf = proteins.get("GyrA", {}).get("confidence", {}).get("mean", "N/A")
        gyrB_conf = proteins.get("GyrB", {}).get("confidence", {}).get("mean", "N/A")
        n_files = len(proteins)
        
        gyrA_str = f"{gyrA_conf:.1f}" if isinstance(gyrA_conf, (int, float)) else "N/A"
        gyrB_str = f"{gyrB_conf:.1f}" if isinstance(gyrB_conf, (int, float)) else "N/A"
        
        print(f"{method_name:<15} {gyrA_str:<15} {gyrB_str:<15} {n_files:<10}")
    
    print()
    print("=" * 70)
    print("Note: pLDDT scores range from 0-100")
    print("  >90: Very high confidence")
    print("  70-90: Confident")
    print("  50-70: Low confidence")
    print("  <50: Very low confidence (likely disordered)")
    print("=" * 70)
    
    # Save results to JSON
    output_file = OUTPUT_DIR / "comparison_results.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {output_file}")

if __name__ == "__main__":
    main()
