#!/usr/bin/env python3
"""
ESMFold Structure Prediction
Meta AI's protein structure prediction model

Usage: python run_esmfold.py
"""

import requests
import os
import sys
from pathlib import Path

# ESMFold API endpoint
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# Sequences (extracted from FASTA)
SEQUENCES = {
    "GyrA": """MTDTTLPPGGDDAVDRIEPVDIQQEMQRSYIDYAMSVIVGRALPEVRDGLKPVHRRVLYA
MYDSGFRPDRSHAKSARSVAETMGNYHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGS
PGNDPAAAMRYTEARLTPLAMEMLREIDEETVDFIPNYDGRVMEPTVLPSRFPNLLANGS
GGIAVGMATNMPPHNLRELAEAVYWALDNHEADEETTLKAVCEKITGPDFPTSGLIVGTQ
GIHDAYTTGRGSIRMRGVAEIEEDSKGRTSLVITELPYQVNHDNFITSIAEQVRDGKIAG
ISNIEDQSSDRVGLRIVVVLKRDAVAKVVLNNLYKHTQLQTSFGANMLSIVDGVPRTLRL
DQLIRLYVNHQLDVIIRRTRYRLRKANERAHILRGLVKALDALDEVIALIRASQTVDIAR
TGLIELLDVDEIQAQAILDMQLRRLAALERQKIIDDLAKIEAEIADLEDILAKPERQRAI
VKDELAEITEKYGDDRRTRIISADGDVADEDLIAREDVVVTITETGYAKRTKTDLYRSQK
RGGKGVQGAGLKQDDIVKHFFVCSTHDWILFFTTKGRVYRAKAYDLPEAARTARGQHVAN
LLAFQPEERIAQVIQIKSYEDAPYLVLATKNGLVKKSKLTEFDSNRSGGLVAVNLRDGDE
LVGAVLCSAEDDLLLVSAHGQSIRFSATDEALRPMGRATSGVQGMRFNGEDDLLSLNVVR
EGTYLLVATSGGYSKRTAIEEYPVQGRGGKGVLTVQYDPRRGSLVGALVVDEESELYAIT
SGGGVIRTIAKQVRKAGRQTKGVRLMNLGEGDTLLAIAHNADEGDADPDEDAAGTTAGE""".replace("\n", ""),
    
    "GyrB": """MAAQKKSAKSEYSADSITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNSVDEAMAG
YATTVEVTMLADGGIQVKDDGRGIPVAMHASGIPTVDVVMTQLHAGGKFDSDSYAVSGGL
HGVGISVVNALSTKLEVEILRDGFEWQQVYTRSEPGTLQKGAATKKTGTTVRFWADPEIF
ETTTYDFETVARRLQEQAFLNKGLTIKLTDERVSDSEVTDEVVSDTAEAPKNAEEQAAES
SAPHKVKNRVFHYPDGLVDFVKHINRTKSAIHTTIVDFSGKGEGHEVEIAMQWNAGYSES
VHTFANTINTHEGGTHEEGFRAALTSVVNKYAKEKKLLKEKDSNLTGDDIREGLAAVISV
KVGEPQFEGQTKTKLGNTEVKSFVQKVCNEQLQHWFDSNPADAKTVVNKAVSSAQARIAA
RKARELVRRKSATDIGGLPGKLADCRSTDPSKSELYVVEGDSAGGSAKSGRDSMFQAILP
LRGKIINVEKARIDRVLKNTEVQAIITALGTGIHDEFDIAKLRYHKIVLMADADVDGQHI
STLLLTLLFRFMRPLVEHGHIFLAQPPLYKLKWQRTQPEFAYSDRERDGLMEAGLKAGKK
INKDDGIQRYKGLGEMDAKELWETTMDPSVRVLRQVTLDDAAAADELFSILMGEDVEARR
SFITRNAKDVRFLDV""".replace("\n", "")
}

def run_esmfold(name, sequence, output_dir):
    """Run ESMFold prediction via API"""
    print(f"\n{'='*60}")
    print(f"Running ESMFold for {name}")
    print(f"Sequence length: {len(sequence)} residues")
    print(f"{'='*60}")
    
    # ESMFold has a limit of ~400 residues for free API
    if len(sequence) > 400:
        print(f"WARNING: Sequence length ({len(sequence)}) exceeds recommended limit (400)")
        print("Truncating to first 400 residues for demo...")
        sequence = sequence[:400]
    
    try:
        print("Submitting to ESMFold API...")
        response = requests.post(
            ESMFOLD_API,
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=300  # 5 minute timeout
        )
        
        if response.status_code == 200:
            output_file = os.path.join(output_dir, f"mabs_{name}_esmfold.pdb")
            with open(output_file, 'w') as f:
                f.write(response.text)
            print(f"SUCCESS: Structure saved to {output_file}")
            return output_file
        else:
            print(f"ERROR: API returned status {response.status_code}")
            print(response.text[:500])
            return None
            
    except requests.exceptions.Timeout:
        print("ERROR: Request timed out (server may be busy)")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        return None

def main():
    output_dir = "output/ai_predictions/esmfold"
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("ESMFold Structure Prediction")
    print("Meta AI - Evolutionary Scale Modeling")
    print("="*60)
    
    results = {}
    for name, seq in SEQUENCES.items():
        result = run_esmfold(name, seq, output_dir)
        results[name] = result
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for name, path in results.items():
        status = "SUCCESS" if path else "FAILED"
        print(f"{name}: {status}")
        if path:
            print(f"  -> {path}")
    
    return results

if __name__ == "__main__":
    main()
