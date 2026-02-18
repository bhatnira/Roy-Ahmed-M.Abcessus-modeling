#!/usr/bin/env python3
"""
Analyze RosettaCM output models for dimeric homology modeling.

This script:
1. Parses the score file
2. Ranks models by total score
3. Identifies the best models
4. Generates analysis report

Usage:
    python analyze_models.py --scores output/scores.sc --models output/models/
"""

import os
import argparse
from pathlib import Path
from typing import List, Dict, Tuple
import statistics

def parse_score_file(score_file: str) -> List[Dict]:
    """
    Parse Rosetta score file.
    
    Returns:
        List of dictionaries with score data for each model
    """
    scores = []
    header = None
    
    with open(score_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('SEQUENCE:'):
                continue
            
            parts = line.split()
            
            if parts[0] == 'SCORE:':
                if header is None:
                    # This is the header line
                    header = parts[1:]  # Skip 'SCORE:'
                else:
                    # This is a data line
                    values = parts[1:]
                    if len(values) == len(header):
                        model_data = {}
                        for i, col in enumerate(header):
                            try:
                                model_data[col] = float(values[i])
                            except ValueError:
                                model_data[col] = values[i]
                        scores.append(model_data)
    
    return scores

def rank_models(scores: List[Dict], sort_by: str = 'total_score') -> List[Dict]:
    """Rank models by specified score term."""
    return sorted(scores, key=lambda x: x.get(sort_by, float('inf')))

def calculate_statistics(scores: List[Dict], term: str) -> Dict:
    """Calculate statistics for a score term."""
    values = [s[term] for s in scores if term in s and isinstance(s[term], (int, float))]
    
    if not values:
        return {}
    
    return {
        'mean': statistics.mean(values),
        'stdev': statistics.stdev(values) if len(values) > 1 else 0,
        'min': min(values),
        'max': max(values),
        'median': statistics.median(values)
    }

def identify_interface_residues(pdb_file: str, distance_cutoff: float = 8.0) -> List[int]:
    """
    Identify interface residues between chains A and B.
    
    Returns residue numbers that are within distance_cutoff of the other chain.
    """
    # Parse coordinates
    chain_a_atoms = []
    chain_b_atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name != 'CA':  # Use only CA atoms for speed
                    continue
                
                chain = line[21]
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                if chain == 'A':
                    chain_a_atoms.append((res_num, x, y, z))
                elif chain == 'B':
                    chain_b_atoms.append((res_num, x, y, z))
    
    # Find interface residues
    interface_a = set()
    interface_b = set()
    
    for res_a, xa, ya, za in chain_a_atoms:
        for res_b, xb, yb, zb in chain_b_atoms:
            dist = ((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2) ** 0.5
            if dist < distance_cutoff:
                interface_a.add(res_a)
                interface_b.add(res_b)
    
    return sorted(interface_a), sorted(interface_b)

def calculate_interface_score(pdb_file: str, scores: Dict) -> float:
    """
    Calculate approximate interface score.
    Uses dG_separated if available, otherwise estimates from total score.
    """
    # If we have interface scores from Rosetta InterfaceAnalyzer
    if 'dG_separated' in scores:
        return scores['dG_separated']
    
    # Otherwise return total score as proxy
    return scores.get('total_score', float('inf'))

def generate_report(scores: List[Dict], output_dir: str, top_n: int = 10):
    """Generate analysis report."""
    
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("RosettaCM Dimeric Homology Modeling - Analysis Report")
    report_lines.append("=" * 70)
    report_lines.append("")
    
    # Overall statistics
    report_lines.append("OVERALL STATISTICS")
    report_lines.append("-" * 40)
    report_lines.append(f"Total models analyzed: {len(scores)}")
    
    if scores:
        total_stats = calculate_statistics(scores, 'total_score')
        report_lines.append(f"\nTotal Score Statistics:")
        report_lines.append(f"  Mean:   {total_stats.get('mean', 'N/A'):.2f}")
        report_lines.append(f"  StdDev: {total_stats.get('stdev', 'N/A'):.2f}")
        report_lines.append(f"  Min:    {total_stats.get('min', 'N/A'):.2f}")
        report_lines.append(f"  Max:    {total_stats.get('max', 'N/A'):.2f}")
        report_lines.append(f"  Median: {total_stats.get('median', 'N/A'):.2f}")
    
    # Top models
    report_lines.append("")
    report_lines.append(f"TOP {top_n} MODELS (by total_score)")
    report_lines.append("-" * 40)
    
    ranked = rank_models(scores, 'total_score')
    
    for i, model in enumerate(ranked[:top_n], 1):
        model_name = model.get('description', f'model_{i}')
        total_score = model.get('total_score', 'N/A')
        
        report_lines.append(f"\n{i}. {model_name}")
        report_lines.append(f"   Total Score: {total_score:.2f}" if isinstance(total_score, float) else f"   Total Score: {total_score}")
        
        # Add additional score terms if available
        for term in ['rms', 'fa_atr', 'fa_rep', 'fa_elec', 'hbond_bb_sc']:
            if term in model:
                report_lines.append(f"   {term}: {model[term]:.2f}")
    
    # Recommendations
    report_lines.append("")
    report_lines.append("RECOMMENDATIONS")
    report_lines.append("-" * 40)
    report_lines.append("")
    
    if ranked:
        best = ranked[0]
        report_lines.append(f"Best model: {best.get('description', 'model_1')}")
        report_lines.append(f"Best score: {best.get('total_score', 'N/A'):.2f}" if isinstance(best.get('total_score'), float) else f"Best score: {best.get('total_score', 'N/A')}")
        report_lines.append("")
        report_lines.append("Suggested next steps:")
        report_lines.append("1. Visually inspect top models in PyMOL/Chimera")
        report_lines.append("2. Run InterfaceAnalyzer to assess dimer interface")
        report_lines.append("3. Consider running molecular dynamics for validation")
        report_lines.append("4. Compare with experimental data if available")
    
    report_lines.append("")
    report_lines.append("=" * 70)
    
    # Print report
    report_text = '\n'.join(report_lines)
    print(report_text)
    
    # Save report
    report_file = os.path.join(output_dir, 'analysis_report.txt')
    with open(report_file, 'w') as f:
        f.write(report_text)
    
    print(f"\nReport saved to: {report_file}")
    
    # Save ranked model list
    ranked_file = os.path.join(output_dir, 'ranked_models.txt')
    with open(ranked_file, 'w') as f:
        f.write("Rank\tModel\tTotal_Score\n")
        for i, model in enumerate(ranked, 1):
            name = model.get('description', f'model_{i}')
            score = model.get('total_score', 'N/A')
            f.write(f"{i}\t{name}\t{score}\n")
    
    print(f"Ranked list saved to: {ranked_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze RosettaCM models')
    parser.add_argument('--scores', default='output/scores.sc',
                        help='Path to score file')
    parser.add_argument('--models', default='output/models/',
                        help='Directory containing output models')
    parser.add_argument('--output', default='output/',
                        help='Output directory for analysis')
    parser.add_argument('--top', type=int, default=10,
                        help='Number of top models to report')
    args = parser.parse_args()
    
    # Check if score file exists
    if not os.path.exists(args.scores):
        print(f"Score file not found: {args.scores}")
        print("Please run RosettaCM first.")
        return
    
    # Parse scores
    print(f"Parsing score file: {args.scores}")
    scores = parse_score_file(args.scores)
    
    if not scores:
        print("No scores found in file.")
        return
    
    print(f"Found {len(scores)} models")
    
    # Generate report
    generate_report(scores, args.output, args.top)

if __name__ == '__main__':
    main()
