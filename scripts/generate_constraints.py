#!/usr/bin/env python3
"""
Generate constraint files for RosettaCM dimeric modeling.

Constraints help maintain:
1. Dimer interface contacts
2. Known structural features
3. Template-derived distance restraints

Usage:
    python generate_constraints.py --template template.pdb --output constraints.cst
"""

import os
import argparse
from pathlib import Path
from typing import List, Tuple, Dict
import math

def read_pdb_coordinates(pdb_file: str) -> Dict[str, Dict[int, Dict[str, Tuple[float, float, float]]]]:
    """
    Read coordinates from a PDB file.
    
    Returns:
        Nested dict: chain -> residue_number -> atom_name -> (x, y, z)
    """
    coords = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                res_num = int(line[22:26])
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                if chain not in coords:
                    coords[chain] = {}
                if res_num not in coords[chain]:
                    coords[chain][res_num] = {}
                
                coords[chain][res_num][atom_name] = (x, y, z)
    
    return coords

def calculate_distance(coord1: Tuple[float, float, float], 
                      coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def find_interface_contacts(coords: Dict, distance_cutoff: float = 8.0) -> List[Dict]:
    """
    Find contacts at the dimer interface.
    
    Returns list of contact dictionaries with:
    - chain1, res1, atom1
    - chain2, res2, atom2
    - distance
    """
    contacts = []
    chains = list(coords.keys())
    
    if len(chains) < 2:
        return contacts
    
    chain_a, chain_b = chains[0], chains[1]
    
    for res1 in coords[chain_a]:
        for atom1, coord1 in coords[chain_a][res1].items():
            # Only use backbone and CB atoms for constraints
            if atom1 not in ['CA', 'CB', 'N', 'C', 'O']:
                continue
            
            for res2 in coords[chain_b]:
                for atom2, coord2 in coords[chain_b][res2].items():
                    if atom2 not in ['CA', 'CB', 'N', 'C', 'O']:
                        continue
                    
                    dist = calculate_distance(coord1, coord2)
                    
                    if dist < distance_cutoff:
                        contacts.append({
                            'chain1': chain_a,
                            'res1': res1,
                            'atom1': atom1,
                            'chain2': chain_b,
                            'res2': res2,
                            'atom2': atom2,
                            'distance': dist
                        })
    
    return contacts

def generate_distance_constraints(contacts: List[Dict], 
                                  tolerance: float = 2.0) -> List[str]:
    """
    Generate Rosetta distance constraints from contacts.
    
    Format:
    AtomPair ATOM1 RES1 ATOM2 RES2 HARMONIC distance tolerance
    """
    constraints = []
    
    for contact in contacts:
        # For dimers, need to adjust residue numbering if chains are concatenated
        res1 = contact['res1']
        res2 = contact['res2']
        dist = contact['distance']
        
        # Use BOUNDED constraint: (distance - tolerance) to (distance + tolerance)
        constraint = (
            f"AtomPair {contact['atom1']} {res1} {contact['atom2']} {res2} "
            f"BOUNDED {max(0, dist - tolerance):.2f} {dist + tolerance:.2f} 0.5 NOE"
        )
        constraints.append(constraint)
    
    return constraints

def generate_harmonic_constraints(contacts: List[Dict],
                                   stdev: float = 1.5) -> List[str]:
    """
    Generate harmonic distance constraints.
    
    Format:
    AtomPair ATOM1 RES1 ATOM2 RES2 HARMONIC distance stdev
    """
    constraints = []
    
    for contact in contacts:
        constraint = (
            f"AtomPair {contact['atom1']} {contact['res1']} "
            f"{contact['atom2']} {contact['res2']} "
            f"HARMONIC {contact['distance']:.2f} {stdev:.2f}"
        )
        constraints.append(constraint)
    
    return constraints

def generate_interface_constraints(coords: Dict, 
                                   distance_cutoff: float = 8.0) -> List[str]:
    """
    Generate constraints specifically for the dimer interface.
    Uses CA-CA distances between chains.
    """
    constraints = []
    chains = list(coords.keys())
    
    if len(chains) < 2:
        return constraints
    
    chain_a, chain_b = chains[0], chains[1]
    
    for res1 in coords[chain_a]:
        if 'CA' not in coords[chain_a][res1]:
            continue
        
        for res2 in coords[chain_b]:
            if 'CA' not in coords[chain_b][res2]:
                continue
            
            dist = calculate_distance(
                coords[chain_a][res1]['CA'],
                coords[chain_b][res2]['CA']
            )
            
            if dist < distance_cutoff:
                # Strong constraint for close contacts
                stdev = 1.0 if dist < 5.0 else 2.0
                constraint = (
                    f"AtomPair CA {res1} CA {res2} "
                    f"HARMONIC {dist:.2f} {stdev:.2f}"
                )
                constraints.append(constraint)
    
    return constraints

def write_constraint_file(constraints: List[str], output_file: str):
    """Write constraints to file."""
    with open(output_file, 'w') as f:
        f.write("# RosettaCM Distance Constraints\n")
        f.write("# Generated for dimeric homology modeling\n")
        f.write("#\n")
        for constraint in constraints:
            f.write(constraint + '\n')
    
    print(f"Wrote {len(constraints)} constraints to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate constraints for RosettaCM')
    parser.add_argument('--template', required=True,
                        help='Template PDB file')
    parser.add_argument('--output', default='constraints.cst',
                        help='Output constraint file')
    parser.add_argument('--distance', type=float, default=8.0,
                        help='Distance cutoff for interface contacts')
    parser.add_argument('--tolerance', type=float, default=2.0,
                        help='Tolerance for distance constraints')
    parser.add_argument('--type', choices=['harmonic', 'bounded'], default='harmonic',
                        help='Type of constraint function')
    args = parser.parse_args()
    
    # Read template coordinates
    print(f"Reading template: {args.template}")
    coords = read_pdb_coordinates(args.template)
    
    chains = list(coords.keys())
    print(f"Found {len(chains)} chains: {chains}")
    
    # Find interface contacts
    contacts = find_interface_contacts(coords, args.distance)
    print(f"Found {len(contacts)} interface contacts")
    
    # Generate constraints
    if args.type == 'harmonic':
        constraints = generate_harmonic_constraints(contacts, args.tolerance)
    else:
        constraints = generate_distance_constraints(contacts, args.tolerance)
    
    # Add interface-specific constraints
    interface_constraints = generate_interface_constraints(coords, args.distance)
    
    # Combine and deduplicate
    all_constraints = list(set(constraints + interface_constraints))
    
    # Write output
    write_constraint_file(all_constraints, args.output)
    
    print("\nConstraint generation complete!")
    print(f"Use with: -cst_fa_file {args.output}")

if __name__ == '__main__':
    main()
