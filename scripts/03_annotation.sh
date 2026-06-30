#!/bin/bash
# =============================================================================
# 03_annotation.sh — Genome Annotation, AMR & Plasmid Profiling
# =============================================================================
#
# USAGE:
#   bash scripts/03_annotation.sh [OPTIONS]
#
# OPTIONS:
#   --contigs PATH     Assembled contigs FASTA file (required)
#   --threads INT      Number of CPU threads (default: 4)
#   --output DIR       Output directory (default: results/step3_annotation)
#   --help             Show this help message
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
CONTIGS=""
THREADS=4
OUTPUT="results/step3_annotation"

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
TOOLS=("prokka" "abricate" "amrfinder" "mob_recon" "mlst")
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
mkdir -p "$OUTPUT/prokka_out"
mkdir -p "$OUTPUT/amr_virulence"
mkdir -p "$OUTPUT/plasmids"

# ---- 1. Sequence Typing (MLST) ----
echo "Running MLST Sequence Typing..."
mlst "$CONTIGS" > "$OUTPUT/mlst_prediction.txt"

# ---- 2. Prokka Core Genome Annotation ----
echo "Running Prokka for genome annotation..."
prokka --outdir "$OUTPUT/prokka_out" \
       --prefix "QA5221" \
       --kingdom "Bacteria" \
       --genus "Escherichia" \
       --species "coli" \
       --strain "QA5221" \
       --force \
       --cpus "$THREADS" \
       "$CONTIGS"

# ---- 3. ABRicate AMR/Virulence Screening ----
echo "Running ABRicate screens..."
abricate --db resfinder "$CONTIGS" > "$OUTPUT/amr_virulence/resfinder_hits.tab"
abricate --db vfdb "$CONTIGS" > "$OUTPUT/amr_virulence/vfdb_hits.tab"
abricate --db card "$CONTIGS" > "$OUTPUT/amr_virulence/card_hits.tab"
abricate --db plasmidfinder "$CONTIGS" > "$OUTPUT/amr_virulence/plasmidfinder_hits.tab"

echo "Generating ABRicate summary report..."
abricate --summary "$OUTPUT/amr_virulence"/*.tab > "$OUTPUT/amr_virulence/summary.tab"

# ---- 4. NCBI AMRFinderPlus Validation ----
echo "Running AMRFinderPlus validation..."
amrfinder -p "$OUTPUT/prokka_out/QA5221.faa" \
          -O "Escherichia" \
          --plus > "$OUTPUT/amr_virulence/amrfinder_hits.tsv"

# ---- 5. MOB-suite Plasmid Reconstruction ----
echo "Running MOB-recon for plasmid reconstruction..."
mob_recon --infile "$CONTIGS" --outdir "$OUTPUT/plasmids"

echo "03_annotation.sh completed successfully. Output saved in: $OUTPUT"
