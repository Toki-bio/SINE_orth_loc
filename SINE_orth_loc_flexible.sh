#!/bin/bash

################################################################################
# SINE_orth_loc - Flexible Input Wrapper
# 
# Modern interface for finding orthologous SINE-containing loci between genomes
# This wrapper provides flexible input handling while calling the original
# SINE_orth_loc.bash script internally.
################################################################################

VERSION="2.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORIGINAL_SCRIPT="${SCRIPT_DIR}/SINE_orth_loc.bash"

################################################################################
# Usage and Help
################################################################################

usage() {
    cat << EOF
SINE_orth_loc v${VERSION} - Find orthologous SINE-containing loci between genomes

USAGE:
    $0 [OPTIONS]

REQUIRED ARGUMENTS:
    -g1, --genome1 FILE         First genome assembly (FASTA format)
    -g2, --genome2 FILE         Second genome assembly (FASTA format)
    -s,  --sine FILE            SINE consensus sequence (FASTA format)
    -b1, --bed1 FILE            BED file with SINE coordinates in genome1
    -b2, --bed2 FILE            BED file with SINE coordinates in genome2

OPTIONAL ARGUMENTS:
    -o,  --output DIR           Output directory (default: sine_results_TIMESTAMP)
    -n1, --name1 NAME           Short name for genome1 (default: derived from filename)
    -n2, --name2 NAME           Short name for genome2 (default: derived from filename)
    -k,  --keep-intermediates   Keep all intermediate files (default: cleanup)
    -t,  --threads NUM          Number of threads (default: auto-detect)
    -h,  --help                 Show this help message
    -v,  --version              Show version

EXAMPLES:
    # Basic usage - let script derive names:
    $0 -g1 mouse_genome.fasta -g2 human_genome.fasta \\
       -s B1_consensus.fa -b1 mouse_B1.bed -b2 human_B1.bed

    # Specify custom names and output directory:
    $0 -g1 GRCm39.fa -g2 GRCh38.fa -s B1.fa \\
       -b1 mouse_elements.bed -b2 human_elements.bed \\
       -n1 mmu -n2 hsa -o results_mouse_vs_human

    # Keep intermediate files for debugging:
    $0 -g1 sp1.fa -g2 sp2.fa -s SINE.fa \\
       -b1 sp1_SINEs.bed -b2 sp2_SINEs.bed --keep-intermediates

OUTPUT FILES:
    The script creates organized output with final results in:
    - {output}/results/statbed_{name1}-{name2}.bed     - Coordinates of validated loci
    - {output}/results/statbed_{name2}-{name1}.bed     - Reciprocal coordinates
    - {output}/results/MP_PM_SINE_{name1}-{name2}.txt  - Summary statistics
    - {output}/results/alignments/                     - Individual alignment files

NOTES:
    - Input FASTA files can have any extension (.fa, .fasta, .fna, etc.)
    - BED files must be 6-column format with strand information
    - SINE consensus header in FASTA must match the filename (without extension)
    - Requires: mafft, esl-alipid, seqkit, bedtools, samtools, bwa, sam2bed, ComPair.sh

EOF
}

################################################################################
# Argument Parsing
################################################################################

GENOME1=""
GENOME2=""
SINE=""
BED1=""
BED2=""
OUTPUT_DIR=""
NAME1=""
NAME2=""
KEEP_INTERMEDIATES=false
THREADS=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g1|--genome1)
            GENOME1="$2"
            shift 2
            ;;
        -g2|--genome2)
            GENOME2="$2"
            shift 2
            ;;
        -s|--sine)
            SINE="$2"
            shift 2
            ;;
        -b1|--bed1)
            BED1="$2"
            shift 2
            ;;
        -b2|--bed2)
            BED2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -n1|--name1)
            NAME1="$2"
            shift 2
            ;;
        -n2|--name2)
            NAME2="$2"
            shift 2
            ;;
        -k|--keep-intermediates)
            KEEP_INTERMEDIATES=true
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -v|--version)
            echo "SINE_orth_loc version ${VERSION}"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

################################################################################
# Input Validation
################################################################################

error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

warning() {
    echo "WARNING: $1" >&2
}

# Check required arguments
[[ -z "$GENOME1" ]] && error_exit "Missing required argument: --genome1"
[[ -z "$GENOME2" ]] && error_exit "Missing required argument: --genome2"
[[ -z "$SINE" ]] && error_exit "Missing required argument: --sine"
[[ -z "$BED1" ]] && error_exit "Missing required argument: --bed1"
[[ -z "$BED2" ]] && error_exit "Missing required argument: --bed2"

# Check if files exist and are readable
[[ ! -f "$GENOME1" ]] && error_exit "Genome1 file not found: $GENOME1"
[[ ! -r "$GENOME1" ]] && error_exit "Genome1 file not readable: $GENOME1"
[[ ! -f "$GENOME2" ]] && error_exit "Genome2 file not found: $GENOME2"
[[ ! -r "$GENOME2" ]] && error_exit "Genome2 file not readable: $GENOME2"
[[ ! -f "$SINE" ]] && error_exit "SINE consensus file not found: $SINE"
[[ ! -r "$SINE" ]] && error_exit "SINE consensus file not readable: $SINE"
[[ ! -f "$BED1" ]] && error_exit "BED1 file not found: $BED1"
[[ ! -r "$BED1" ]] && error_exit "BED1 file not readable: $BED1"
[[ ! -f "$BED2" ]] && error_exit "BED2 file not found: $BED2"
[[ ! -r "$BED2" ]] && error_exit "BED2 file not readable: $BED2"

# Check if original script exists
[[ ! -f "$ORIGINAL_SCRIPT" ]] && error_exit "Original script not found: $ORIGINAL_SCRIPT"
[[ ! -x "$ORIGINAL_SCRIPT" ]] && error_exit "Original script not executable: $ORIGINAL_SCRIPT"

# Validate file formats
validate_fasta() {
    local file=$1
    local name=$2
    if ! head -n1 "$file" | grep -q "^>"; then
        error_exit "$name does not appear to be a valid FASTA file (no '>' header)"
    fi
}

validate_bed() {
    local file=$1
    local name=$2
    local columns=$(head -n1 "$file" | awk '{print NF}')
    if [[ $columns -lt 6 ]]; then
        error_exit "$name must have at least 6 columns (BED6 format with strand info)"
    fi
}

echo "Validating input files..."
validate_fasta "$GENOME1" "Genome1"
validate_fasta "$GENOME2" "Genome2"
validate_fasta "$SINE" "SINE consensus"
validate_bed "$BED1" "BED1"
validate_bed "$BED2" "BED2"

# Check required tools
echo "Checking required tools..."
REQUIRED_TOOLS=("mafft" "esl-alipid" "seqkit" "bedtools" "samtools" "bwa" "sam2bed")
MISSING_TOOLS=()

for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

# Check for ComPair.sh
if ! command -v ComPair.sh &> /dev/null && [[ ! -f "${SCRIPT_DIR}/ComPair.sh" ]]; then
    MISSING_TOOLS+=("ComPair.sh")
fi

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    error_exit "Missing required tools: ${MISSING_TOOLS[*]}"
fi

################################################################################
# Derive Names and Setup
################################################################################

# Auto-generate names if not provided
if [[ -z "$NAME1" ]]; then
    NAME1=$(basename "$GENOME1" | sed 's/\.[^.]*$//' | head -c 3)
    echo "Auto-derived name1: $NAME1"
fi

if [[ -z "$NAME2" ]]; then
    NAME2=$(basename "$GENOME2" | sed 's/\.[^.]*$//' | head -c 3)
    echo "Auto-derived name2: $NAME2"
fi

# Validate names (3 characters, alphanumeric)
if [[ ! "$NAME1" =~ ^[a-zA-Z0-9]{3}$ ]]; then
    error_exit "Name1 must be exactly 3 alphanumeric characters. Current: '$NAME1'"
fi

if [[ ! "$NAME2" =~ ^[a-zA-Z0-9]{3}$ ]]; then
    error_exit "Name2 must be exactly 3 alphanumeric characters. Current: '$NAME2'"
fi

# Extract SINE name from filename
SINE_NAME=$(basename "$SINE" | sed 's/\.[^.]*$//')

# Validate SINE header matches filename
SINE_HEADER=$(head -n1 "$SINE" | sed 's/^>//' | awk '{print $1}')
if [[ "$SINE_HEADER" != "$SINE_NAME" ]]; then
    error_exit "SINE FASTA header (>$SINE_HEADER) must match filename ($SINE_NAME)"
fi

# Setup output directory
if [[ -z "$OUTPUT_DIR" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    OUTPUT_DIR="sine_results_${TIMESTAMP}"
fi

if [[ -e "$OUTPUT_DIR" ]]; then
    warning "Output directory already exists: $OUTPUT_DIR"
    read -p "Continue and potentially overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted by user"
        exit 0
    fi
fi

mkdir -p "$OUTPUT_DIR"
WORK_DIR="${OUTPUT_DIR}/work"
RESULTS_DIR="${OUTPUT_DIR}/results"
ALIGNMENTS_DIR="${RESULTS_DIR}/alignments"
mkdir -p "$WORK_DIR" "$RESULTS_DIR" "$ALIGNMENTS_DIR"

# Get absolute paths
GENOME1=$(readlink -f "$GENOME1")
GENOME2=$(readlink -f "$GENOME2")
SINE=$(readlink -f "$SINE")
BED1=$(readlink -f "$BED1")
BED2=$(readlink -f "$BED2")
OUTPUT_DIR=$(readlink -f "$OUTPUT_DIR")
WORK_DIR=$(readlink -f "$WORK_DIR")
RESULTS_DIR=$(readlink -f "$RESULTS_DIR")
ALIGNMENTS_DIR=$(readlink -f "$ALIGNMENTS_DIR")

################################################################################
# Setup Working Directory
################################################################################

echo ""
echo "=========================================="
echo "SINE Orthologous Loci Finder"
echo "=========================================="
echo "Genome 1:     $GENOME1 (${NAME1})"
echo "Genome 2:     $GENOME2 (${NAME2})"
echo "SINE:         $SINE (${SINE_NAME})"
echo "BED 1:        $BED1"
echo "BED 2:        $BED2"
echo "Output:       $OUTPUT_DIR"
echo "Work dir:     $WORK_DIR"
echo "=========================================="
echo ""

cd "$WORK_DIR" || error_exit "Cannot change to work directory"

# Create symlinks with required naming convention
echo "Setting up working directory..."

ln -sf "$GENOME1" "${NAME1}.bnk"
ln -sf "$GENOME2" "${NAME2}.bnk"
ln -sf "$SINE" "${SINE_NAME}.fa"
ln -sf "$BED1" "${NAME1}-${SINE_NAME}.bed"
ln -sf "$BED2" "${NAME2}-${SINE_NAME}.bed"

# Copy ComPair.sh to work directory if it exists
if [[ -f "${SCRIPT_DIR}/ComPair.sh" ]]; then
    cp "${SCRIPT_DIR}/ComPair.sh" "$WORK_DIR/"
    chmod +x "$WORK_DIR/ComPair.sh"
fi

# Add current directory to PATH for ComPair.sh
export PATH="${WORK_DIR}:${PATH}"

################################################################################
# Run Original Script
################################################################################

echo ""
echo "Starting analysis..."
echo "This may take a while depending on genome sizes..."
echo ""

# Save original directory
ORIGINAL_DIR=$(pwd)

# Set threads if specified
if [[ -n "$THREADS" ]]; then
    export THREADS
fi

# Run the original script
set -o pipefail
if bash "$ORIGINAL_SCRIPT" "${NAME1}.bnk" "${NAME2}.bnk" "${SINE_NAME}.fa" 2>&1 | tee "${OUTPUT_DIR}/pipeline.log"; then
    echo ""
    echo "Analysis completed successfully!"
else
    EXIT_CODE=$?
    error_exit "Pipeline failed with exit code $EXIT_CODE. Check ${OUTPUT_DIR}/pipeline.log"
fi

################################################################################
# Organize Results
################################################################################

echo ""
echo "Organizing results..."

# Move final result files to results directory
[[ -f "statbed_${NAME1}-${NAME2}.bed" ]] && mv "statbed_${NAME1}-${NAME2}.bed" "$RESULTS_DIR/"
[[ -f "statbed_${NAME2}-${NAME1}.bed" ]] && mv "statbed_${NAME2}-${NAME1}.bed" "$RESULTS_DIR/"
[[ -f "MP_PM_SINE_${NAME1}-${NAME2}.txt" ]] && mv "MP_PM_SINE_${NAME1}-${NAME2}.txt" "$RESULTS_DIR/"
[[ -f "stat_doubles_${NAME1}-${NAME2}" ]] && mv "stat_doubles_${NAME1}-${NAME2}" "$RESULTS_DIR/"
[[ -f "stat_multi_${NAME1}-${NAME2}" ]] && mv "stat_multi_${NAME1}-${NAME2}" "$RESULTS_DIR/"

# Move alignment files
for ext in PM MP SINE; do
    if ls *.${ext} &> /dev/null; then
        mv *.${ext} "$ALIGNMENTS_DIR/" 2>/dev/null || true
    fi
done

# Create summary report
SUMMARY_FILE="${RESULTS_DIR}/summary.txt"
cat > "$SUMMARY_FILE" << EOF
SINE Orthologous Loci Analysis Summary
======================================

Run Date: $(date)
Pipeline Version: ${VERSION}

Input Files:
  Genome 1: $GENOME1
  Genome 2: $GENOME2
  SINE:     $SINE
  BED 1:    $BED1
  BED 2:    $BED2

Analysis Names:
  Species 1: ${NAME1}
  Species 2: ${NAME2}
  SINE:      ${SINE_NAME}

Results:
EOF

# Add statistics if available
if [[ -f "${RESULTS_DIR}/MP_PM_SINE_${NAME1}-${NAME2}.txt" ]]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Classification Summary:" >> "$SUMMARY_FILE"
    cat "${RESULTS_DIR}/MP_PM_SINE_${NAME1}-${NAME2}.txt" >> "$SUMMARY_FILE"
fi

# Count alignment files
PM_COUNT=$(ls "$ALIGNMENTS_DIR"/*.PM 2>/dev/null | wc -l)
MP_COUNT=$(ls "$ALIGNMENTS_DIR"/*.MP 2>/dev/null | wc -l)
SINE_COUNT=$(ls "$ALIGNMENTS_DIR"/*.SINE 2>/dev/null | wc -l)

cat >> "$SUMMARY_FILE" << EOF

Alignment Files:
  Plus-Minus (PM):  ${PM_COUNT} files
  Minus-Plus (MP):  ${MP_COUNT} files
  Both have SINE:   ${SINE_COUNT} files

Output Directory Structure:
  ${OUTPUT_DIR}/
    ├── results/
    │   ├── statbed_${NAME1}-${NAME2}.bed    (Genomic coordinates)
    │   ├── statbed_${NAME2}-${NAME1}.bed    (Reciprocal coordinates)
    │   ├── MP_PM_SINE_${NAME1}-${NAME2}.txt (Summary statistics)
    │   ├── stat_doubles_${NAME1}-${NAME2}   (Detailed stats - doubles)
    │   ├── stat_multi_${NAME1}-${NAME2}     (Detailed stats - multis)
    │   ├── alignments/                      (Individual alignments)
    │   └── summary.txt                      (This file)
    ├── work/                                (Intermediate files)
    └── pipeline.log                         (Complete log)

EOF

# Cleanup intermediate files if requested
if [[ "$KEEP_INTERMEDIATES" == false ]]; then
    echo "Cleaning up intermediate files..."
    cat >> "$SUMMARY_FILE" << EOF
Intermediate files were removed (run with --keep-intermediates to preserve)
EOF
    
    # Keep only essential intermediate files
    cd "$WORK_DIR"
    find . -type f ! -name "*.bed" ! -name "*.bnk" ! -name "*.fa" ! -name "*.fai" ! -name "*.log" -delete 2>/dev/null || true
else
    cat >> "$SUMMARY_FILE" << EOF
Intermediate files preserved in: ${WORK_DIR}/
EOF
fi

################################################################################
# Final Report
################################################################################

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
cat "$SUMMARY_FILE"
echo ""
echo "All results saved to: $OUTPUT_DIR"
echo "View summary: cat ${SUMMARY_FILE}"
echo "=========================================="
