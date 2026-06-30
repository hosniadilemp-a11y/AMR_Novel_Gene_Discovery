#!/bin/bash
# =============================================================================
# 04_mge_detection.sh — Mobile Genetic Element (MGE) Detection Pipeline
# =============================================================================
#
# USAGE:
#   bash scripts/04_mge_detection.sh [OPTIONS]
#
# OPTIONS:
#   --contigs PATH     Assembled contigs FASTA file (required)
#   --threads INT      Number of CPU threads to use (default: 8)
#   --output DIR       Output directory (default: results/step4_mge)
#   --help             Show this help message
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
CONTIGS=""
THREADS=8
OUTPUT="results/step4_mge"

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --contigs)   CONTIGS="$2";   shift 2 ;;
        --threads)   THREADS="$2";   shift 2 ;;
        --output)    OUTPUT="$2";    shift 2 ;;
        --help)      head -20 "$0";  exit 0 ;;
        *)           echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---- Check dependencies ----
TOOLS=("isescan.py" "integron_finder" "python3")
for tool in "${TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

# ---- Validate required inputs ----
if [[ -z "$CONTIGS" ]]; then
    echo "ERROR: --contigs is required."
    exit 1
fi

if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Contigs file '$CONTIGS' does not exist."
    exit 1
fi

# ---- Create output directories ----
mkdir -p "$OUTPUT/isescan_out"
mkdir -p "$OUTPUT/integron_finder_out"

# ---- 1. Run ISEScan to detect Insertion Sequences ----
echo "Running ISEScan..."
isescan.py --seqfile "$CONTIGS" \
           --output "$OUTPUT/isescan_out" \
           --nthread "$THREADS" \
           --removeShortIS

# ---- 2. Run Isolated IntegronFinder to prevent contamination ----
echo "Running isolated IntegronFinder to prevent multi-fasta database leakage..."
python3 scripts/run_integron_finder_clean.py \
    --contigs "$CONTIGS" \
    --output "$OUTPUT/integron_finder_out" \
    --threads "$THREADS"

echo "04_mge_detection.sh completed successfully. Outputs saved in: $OUTPUT"
