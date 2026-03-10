#!/bin/bash
# =============================================================================
# RoseTTAFold Installation and Run Script
# Rosetta + Deep Learning Hybrid - Baker Lab
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="$PROJECT_DIR/output/ai_predictions/rosettafold"
ROSETTAFOLD_DIR="$HOME/RoseTTAFold"
ENV_NAME="rosettafold_env"

echo "============================================================"
echo "RoseTTAFold Structure Prediction (Local GPU)"
echo "Baker Lab - Rosetta + Deep Learning Hybrid"
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

# Install RoseTTAFold if not present
if [ ! -d "$ROSETTAFOLD_DIR" ]; then
    echo "Cloning RoseTTAFold repository..."
    cd $HOME
    git clone https://github.com/RosettaCommons/RoseTTAFold.git
    cd RoseTTAFold
    
    echo "Creating conda environment..."
    # Create environment
    conda env create -f RoseTTAFold-linux.yml -n $ENV_NAME
    
    echo "Downloading model weights..."
    # Download weights (requires license agreement)
    wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
    tar xfz weights.tar.gz
    
    # Download databases (UniRef30, BFD)
    echo "Note: Database download requires ~500GB space"
    # ./install_dependencies.sh
fi

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

cd $ROSETTAFOLD_DIR

echo ""
echo "Running RoseTTAFold predictions..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run predictions
for protein in GyrA GyrB; do
    echo "=== Predicting ${protein} structure ==="
    
    # RoseTTAFold run command
    ./run_e2e_ver.sh \
        "$OUTPUT_DIR/input/${protein}.fasta" \
        "$OUTPUT_DIR/${protein}"
done

echo ""
echo "============================================================"
echo "RoseTTAFold predictions complete!"
echo "Results in: $OUTPUT_DIR"
echo "============================================================"
