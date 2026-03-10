#!/bin/bash
# =============================================================================
# AlphaFold2 Local Installation and Run Script
# DeepMind's AlphaFold2 - State-of-the-art accuracy
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="$PROJECT_DIR/output/ai_predictions/alphafold2"
ALPHAFOLD_DIR="$HOME/alphafold"
DATA_DIR="$HOME/alphafold_data"

echo "============================================================"
echo "AlphaFold2 Structure Prediction (Local GPU)"
echo "DeepMind - State-of-the-art protein structure prediction"
echo "============================================================"

# Check NVIDIA GPU
nvidia-smi

# Create input FASTA
mkdir -p "$OUTPUT_DIR/input"

cat > "$OUTPUT_DIR/input/GyrA.fasta" << 'EOF'
>B1ME58_GyrA_Mabscessus
MTDTTLPPGGDDAVDRIEPVDIQQEMQRSYIDYAMSVIVGRALPEVRDGLKPVHRRVLYAMYDSGFRPDRSHAKSARSVAETMGNYHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGSPGNDPAAAMRYTEARLTPLAMEMLREIDEETVDFIPNYDGRVMEPTVLPSRFPNLLANGSGGIAVGMATNMPPHNLRELAEAVYWALDNHEADEETTLKAVCEKITGPDFPTSGLIVGTQGIHDAYTTGRGSIRMRGVAEIEEDSKGRTSLVITELPYQVNHDNFITSIAEQVRDGKIAGISNIEDQSSDRVGLRIVVVLKRDAVAKVVLNNLYKHTQLQTSFGANMLSIVDGVPRTLRLDQLIRLYVNHQLDVIIRRTRYRLRKANERAHILRGLVKALDALDEVIALIRASQTVDIARTGLIELLDVDEIQAQAILDMQLRRLAALERQKIIDDLAKIEAEIADLEDILAKPERQRAIVKDELAEITEKYGDDRRTRIISADGDVADEDLIAREDVVVTITETGYAKRTKTDLYRSQKRGGKGVQGAGLKQDDIVKHFFVCSTHDWILFFTTKGRVYRAKAYDLPEAARTARGQHVANLLAFQPEERIAQVIQIKSYEDAPYLVLATKNGLVKKSKLTEFDSNRSGGLVAVNLRDGDELVGAVLCSAEDDLLLVSAHGQSIRFSATDEALRPMGRATSGVQGMRFNGEDDLLSLNVVREGTYLLVATSGGYSKRTAIEEYPVQGRGGKGVLTVQYDPRRGSLVGALVVDEESELYAITSGGGVIRTIAKQVRKAGRQTKGVRLMNLGEGDTLLAIAHNADEGDADPDEDAAGTTAGE
EOF

cat > "$OUTPUT_DIR/input/GyrB.fasta" << 'EOF'
>B1ME45_GyrB_Mabscessus
MAAQKKSAKSEYSADSITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNSVDEAMAGYATTVEVTMLADGGIQVKDDGRGIPVAMHASGIPTVDVVMTQLHAGGKFDSDSYAVSGGLHGVGISVVNALSTKLEVEILRDGFEWQQVYTRSEPGTLQKGAATKKTGTTVRFWADPEIFETTTYDFETVARRLQEQAFLNKGLTIKLTDERVSDSEVTDEVVSDTAEAPKNAEEQAAESSAPHKVKNRVFHYPDGLVDFVKHINRTKSAIHTTIVDFSGKGEGHEVEIAMQWNAGYSESВHTFANTINTHEGGTHEEGFRAALTSVVNKYAKEKKLLKEKDSNLTGDDIREGLAAVISVKVGEPQFEGQTKTKLGNTEVKSFVQKVCNEQLQHWFDSNPADAKTVVNKAVSSAQARIAARKARELVRRKSATDIGGLPGKLADCRSTDPSKSELYVVEGDSAGGSAKSGRDSMFQAILPLRGKIINVEKARIDRVLKNTEVQAIITALGTGIHDEFDIAKLRYHKIVLMADADVDGQHISTLLLTLLFRFMRPLVEHGHIFLAQPPLYKLKWQRTQPEFAYSDRERDGLMEAGLKAGKKINKDDGIQRYKGLGEMDAKELWETTMDPSVRVLRQVTLDDAAAADELFSILMGEDVEARRSFITRNAKDVRFLDV
EOF

# Install AlphaFold if not present
if [ ! -d "$ALPHAFOLD_DIR" ]; then
    echo "Cloning AlphaFold repository..."
    cd $HOME
    git clone https://github.com/deepmind/alphafold.git
    cd alphafold
    
    echo "Creating conda environment..."
    conda env create -f docker/conda/environment.yml -n alphafold_env
fi

# Check if databases exist
if [ ! -d "$DATA_DIR" ]; then
    echo ""
    echo "WARNING: AlphaFold databases not found at $DATA_DIR"
    echo "You need to download the databases (requires ~2.5TB disk space):"
    echo ""
    echo "  bash $ALPHAFOLD_DIR/scripts/download_all_data.sh $DATA_DIR"
    echo ""
    echo "Alternative: Use ColabFold which uses MMseqs2 server (no local DB needed)"
    exit 1
fi

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate alphafold_env

echo ""
echo "Running AlphaFold2 predictions..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run predictions using AlphaFold
for protein in GyrA GyrB; do
    echo "=== Predicting ${protein} structure ==="
    python $ALPHAFOLD_DIR/run_alphafold.py \
        --fasta_paths="$OUTPUT_DIR/input/${protein}.fasta" \
        --output_dir="$OUTPUT_DIR/${protein}" \
        --data_dir="$DATA_DIR" \
        --model_preset=monomer \
        --db_preset=full_dbs \
        --max_template_date=2026-02-01 \
        --use_gpu_relax \
        --benchmark=false
done

echo ""
echo "============================================================"
echo "AlphaFold2 predictions complete!"
echo "Results in: $OUTPUT_DIR"
echo "============================================================"
