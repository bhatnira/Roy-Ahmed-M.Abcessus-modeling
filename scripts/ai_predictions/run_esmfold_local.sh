#!/bin/bash
# =============================================================================
# ESMFold Local Setup and Run Script
# Meta AI - Evolutionary Scale Modeling (esm)
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="$PROJECT_DIR/output/ai_predictions/esmfold"
ENV_NAME="esmfold_env"

echo "============================================================"
echo "ESMFold Structure Prediction (Local GPU)"
echo "Using: NVIDIA RTX PRO 5000 Blackwell (48GB VRAM)"
echo "============================================================"

# Check if environment exists
if ! conda env list | grep -q "$ENV_NAME"; then
    echo "Creating conda environment: $ENV_NAME"
    conda create -n $ENV_NAME python=3.10 -y
fi

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install dependencies if not present
if ! python -c "import esm" 2>/dev/null; then
    echo "Installing ESMFold dependencies..."
    pip install torch --index-url https://download.pytorch.org/whl/cu121
    pip install fair-esm
    pip install "fair-esm[esmfold]"
    pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
    pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
fi

echo ""
echo "Running ESMFold predictions..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Create Python script for predictions
python << 'PYTHON_SCRIPT'
import torch
import esm
import os
import sys

# Set output directory
output_dir = os.environ.get('OUTPUT_DIR', 'output/ai_predictions/esmfold')
os.makedirs(f"{output_dir}/GyrA", exist_ok=True)
os.makedirs(f"{output_dir}/GyrB", exist_ok=True)

# Sequences
sequences = {
    "GyrA": "MTDTTLPPGGDDAVDRIEPVDIQQEMQRSYIDYAMSVIVGRALPEVRDGLKPVHRRVLYAMYDSGFRPDRSHAKSARSVAETMGNYHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGSPGNDPAAAMRYTEARLTPLAMEMLREIDEETVDFIPNYDGRVMEPTVLPSRFPNLLANGSGGIAVGMATNMPPHNLRELAEAVYWALDNHEADEETTLKAVCEKITGPDFPTSGLIVGTQGIHDAYTTGRGSIRMRGVAEIEEDSKGRTSLVITELPYQVNHDNFITSIAEQVRDGKIAGISNIEDQSSDRVGLRIVVVLKRDAVAKVVLNNLYKHTQLQTSFGANMLSIVDGVPRTLRLDQLIRLYVNHQLDVIIRRTRYRLRKANERAHILRGLVKALDALDEVIALIRASQTVDIARTGLIELLDVDEIQAQAILDMQLRRLAALERQKIIDDLAKIEAEIADLEDILAKPERQRAIVKDELAEITEKYGDDRRTRIISADGDVADEDLIAREDVVVTITETGYAKRTKTDLYRSQKRGGKGVQGAGLKQDDIVKHFFVCSTHDWILFFTTKGRVYRAKAYDLPEAARTARGQHVANLLAFQPEERIAQVIQIKSYEDAPYLVLATKNGLVKKSKLTEFDSNRSGGLVAVNLRDGDELVGAVLCSAEDDLLLVSAHGQSIRFSATDEALRPMGRATSGVQGMRFNGEDDLLSLNVVREGTYLLVATSGGYSKRTAIEEYPVQGRGGKGVLTVQYDPRRGSLVGALVVDEESELYAITSGGGVIRTIAKQVRKAGRQTKGVRLMNLGEGDTLLAIAHNADEGDADPDEDAAGTTAGE",
    "GyrB": "MAAQKKSAKSEYSADSITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNSVDEAMAGYATTVEVTMLADGGIQVKDDGRGIPVAMHASGIPTVDVVMTQLHAGGKFDSDSYAVSGGLHGVGISVVNALSTKLEVEILRDGFEWQQVYTRSEPGTLQKGAATKKTGTTVRFWADPEIFETTTYDFETVARRLQEQAFLNKGLTIKLTDERVSDSEVTDEVVSDTAEAPKNAEEQAAESSAPHKVKNRVFHYPDGLVDFVKHINRTKSAIHTTIVDFSGKGEGHEVEIAMQWNAGYSESВHTFANTINTHEGGTHEEGFRAALTSVVNKYAKEKKLLKEKDSNLTGDDIREGLAAVISVKVGEPQFEGQTKTKLGNTEVKSFVQKVCNEQLQHWFDSNPADAKTVVNKAVSSAQARIAARKARELVRRKSATDIGGLPGKLADCRSTDPSKSELYVVEGDSAGGSAKSGRDSMFQAILPLRGKIINVEKARIDRVLKNTEVQAIITALGTGIHDEFDIAKLRYHKIVLMADADVDGQHISTLLLTLLFRFMRPLVEHGHIFLAQPPLYKLKWQRTQPEFAYSDRERDGLMEAGLKAGKKINKDDGIQRYKGLGEMDAKELWETTMDPSVRVLRQVTLDDAAAADELFSILMGEDVEARRSFITRNAKDVRFLDV"
}

print("Loading ESMFold model...")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# Load model
model = esm.pretrained.esmfold_v1()
model = model.eval()
if torch.cuda.is_available():
    model = model.cuda()

# Set chunk size for long sequences
model.set_chunk_size(128)

for name, seq in sequences.items():
    print(f"\n{'='*60}")
    print(f"Predicting structure for {name}")
    print(f"Sequence length: {len(seq)} residues")
    print(f"{'='*60}")
    
    try:
        with torch.no_grad():
            output = model.infer_pdb(seq)
        
        output_file = f"{output_dir}/{name}/mabs_{name}_esmfold.pdb"
        with open(output_file, 'w') as f:
            f.write(output)
        print(f"SUCCESS: Saved to {output_file}")
        
        # Also extract pLDDT scores
        print(f"Structure predicted successfully!")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*60)
print("ESMFold predictions complete!")
print("="*60)
PYTHON_SCRIPT

echo ""
echo "ESMFold finished!"
echo "Results in: $OUTPUT_DIR"
