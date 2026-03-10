#!/bin/bash
# =============================================================================
# ColabFold Local Installation and Run Script
# LocalColabFold - Fast AlphaFold2 predictions with MMseqs2
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="$PROJECT_DIR/output/ai_predictions/colabfold"
ENV_NAME="colabfold_env"
COLABFOLD_DIR="$HOME/localcolabfold"

echo "============================================================"
echo "ColabFold Structure Prediction (Local GPU)"
echo "AlphaFold2 + MMseqs2 for fast MSA generation"
echo "============================================================"

# Create input FASTA files
mkdir -p "$PROJECT_DIR/output/ai_predictions/colabfold/input"

cat > "$PROJECT_DIR/output/ai_predictions/colabfold/input/GyrA.fasta" << 'EOF'
>B1ME58_GyrA_Mabscessus
MTDTTLPPGGDDAVDRIEPVDIQQEMQRSYIDYAMSVIVGRALPEVRDGLKPVHRRVLYAMYDSGFRPDRSHAKSARSVAETMGNYHPHGDASIYDTLVRMAQPWSLRYPLVDGQGNFGSPGNDPAAAMRYTEARLTPLAMEMLREIDEETVDFIPNYDGRVMEPTVLPSRFPNLLANGSGGIAVGMATNMPPHNLRELAEAVYWALDNHEADEETTLKAVCEKITGPDFPTSGLIVGTQGIHDAYTTGRGSIRMRGVAEIEEDSKGRTSLVITELPYQVNHDNFITSIAEQVRDGKIAGISNIEDQSSDRVGLRIVVVLKRDAVAKVVLNNLYKHTQLQTSFGANMLSIVDGVPRTLRLDQLIRLYVNHQLDVIIRRTRYRLRKANERAHILRGLVKALDALDEVIALIRASQTVDIARTGLIELLDVDEIQAQAILDMQLRRLAALERQKIIDDLAKIEAEIADLEDILAKPERQRAIVKDELAEITEKYGDDRRTRIISADGDVADEDLIAREDVVVTITETGYAKRTKTDLYRSQKRGGKGVQGAGLKQDDIVKHFFVCSTHDWILFFTTKGRVYRAKAYDLPEAARTARGQHVANLLAFQPEERIAQVIQIKSYEDAPYLVLATKNGLVKKSKLTEFDSNRSGGLVAVNLRDGDELVGAVLCSAEDDLLLVSAHGQSIRFSATDEALRPMGRATSGVQGMRFNGEDDLLSLNVVREGTYLLVATSGGYSKRTAIEEYPVQGRGGKGVLTVQYDPRRGSLVGALVVDEESELYAITSGGGVIRTIAKQVRKAGRQTKGVRLMNLGEGDTLLAIAHNADEGDADPDEDAAGTTAGE
EOF

cat > "$PROJECT_DIR/output/ai_predictions/colabfold/input/GyrB.fasta" << 'EOF'
>B1ME45_GyrB_Mabscessus
MAAQKKSAKSEYSADSITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNSVDEAMAGYATTVEVTMLADGGIQVKDDGRGIPVAMHASGIPTVDVVMTQLHAGGKFDSDSYAVSGGLHGVGISVVNALSTKLEVEILRDGFEWQQVYTRSEPGTLQKGAATKKTGTTVRFWADPEIFETTTYDFETVARRLQEQAFLNKGLTIKLTDERVSDSEVTDEVVSDTAEAPKNAEEQAAESSAPHKVKNRVFHYPDGLVDFVKHINRTKSAIHTTIVDFSGKGEGHEVEIAMQWNAGYSESВHTFANTINTHEGGTHEEGFRAALTSVVNKYAKEKKLLKEKDSNLTGDDIREGLAAVISVKVGEPQFEGQTKTKLGNTEVKSFVQKVCNEQLQHWFDSNPADAKTVVNKAVSSAQARIAARKARELVRRKSATDIGGLPGKLADCRSTDPSKSELYVVEGDSAGGSAKSGRDSMFQAILPLRGKIINVEKARIDRVLKNTEVQAIITALGTGIHDEFDIAKLRYHKIVLMADADVDGQHISTLLLTLLFRFMRPLVEHGHIFLAQPPLYKLKWQRTQPEFAYSDRERDGLMEAGLKAGKKINKDDGIQRYKGLGEMDAKELWETTMDPSVRVLRQVTLDDAAAADELFSILMGEDVEARRSFITRNAKDVRFLDV
EOF

# Install LocalColabFold if not present
if [ ! -d "$COLABFOLD_DIR" ]; then
    echo "Installing LocalColabFold..."
    cd $HOME
    wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
    bash install_colabbatch_linux.sh
fi

# Add to PATH
export PATH="$COLABFOLD_DIR/colabfold-conda/bin:$PATH"

echo ""
echo "Running ColabFold predictions..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run GyrA prediction
echo "=== Predicting GyrA structure ==="
colabfold_batch \
    "$PROJECT_DIR/output/ai_predictions/colabfold/input/GyrA.fasta" \
    "$OUTPUT_DIR/GyrA" \
    --num-recycle 3 \
    --num-models 5 \
    --amber \
    --use-gpu-relax

# Run GyrB prediction  
echo "=== Predicting GyrB structure ==="
colabfold_batch \
    "$PROJECT_DIR/output/ai_predictions/colabfold/input/GyrB.fasta" \
    "$OUTPUT_DIR/GyrB" \
    --num-recycle 3 \
    --num-models 5 \
    --amber \
    --use-gpu-relax

echo ""
echo "============================================================"
echo "ColabFold predictions complete!"
echo "Results in: $OUTPUT_DIR"
echo "============================================================"
