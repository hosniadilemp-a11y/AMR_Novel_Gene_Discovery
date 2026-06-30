#!/bin/bash
# =============================================================================
# 03b_pangenome.sh — Pangenome Clustering and Phylogeny Pipeline
# =============================================================================
#
# USAGE:
#   bash scripts/03b_pangenome.sh [OPTIONS]
#
# OPTIONS:
#   --focal-gff PATH   GFF file of the focus strain (required)
#   --accessions PATH  List of NCBI GenBank accessions to download (optional)
#   --threads INT      Number of CPU threads to use (default: 8)
#   --output DIR       Output directory (default: results/step3_pangenome)
#   --help             Show this help message
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
FOCAL_GFF=""
ACCESSIONS=""
THREADS=8
OUTPUT="results/step3_pangenome"

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --focal-gff)  FOCAL_GFF="$2";  shift 2 ;;
        --accessions) ACCESSIONS="$2"; shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
        --output)     OUTPUT="$2";     shift 2 ;;
        --help)       head -20 "$0";   exit 0 ;;
        *)            echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---- Check dependencies ----
TOOLS=("panaroo" "iqtree" "mafft" "python3" "curl")
for tool in "${TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

# ---- Validate required inputs ----
if [[ -z "$FOCAL_GFF" ]]; then
    echo "ERROR: --focal-gff is required."
    exit 1
fi

if [[ ! -f "$FOCAL_GFF" ]]; then
    echo "ERROR: Focus GFF file '$FOCAL_GFF' does not exist."
    exit 1
fi

# ---- Create output directories ----
mkdir -p "$OUTPUT/reference_genomes"
mkdir -p "$OUTPUT/reference_prokka"
mkdir -p "$OUTPUT/panaroo_out"
mkdir -p "$OUTPUT/iqtree_out"

# ---- 1. Download reference genomes ----
echo "Downloading reference assemblies from NCBI..."
DOWNLOAD_OPTS=""
if [[ -n "$ACCESSIONS" ]]; then
    DOWNLOAD_OPTS="--accessions $ACCESSIONS"
fi

python3 scripts/download_genomes.py \
    $DOWNLOAD_OPTS \
    --output "$OUTPUT/reference_genomes"

# ---- 2. Annotate reference genomes uniformly ----
echo "Running parallel Prokka annotation on reference genomes..."
# We run parallel Prokka workers based on available threads
WORKERS=$(( THREADS / 4 ))
if [[ $WORKERS -lt 1 ]]; then WORKERS=1; fi

python3 scripts/annotate_references.py \
    --indir "$OUTPUT/reference_genomes" \
    --outdir "$OUTPUT/reference_prokka" \
    --workers "$WORKERS" \
    --cpus-per-worker 4

# ---- 3. Compile GFF list for Panaroo ----
echo "Compiling GFF files list for Panaroo..."
GFF_LIST="$OUTPUT/gff_list.txt"
> "$GFF_LIST"

# Add focal strain GFF
echo "$FOCAL_GFF" >> "$GFF_LIST"

# Add all reference strain GFF files
for gff in "$OUTPUT/reference_prokka"/*/*.gff; do
    if [[ -f "$gff" ]]; then
        echo "$gff" >> "$GFF_LIST"
    fi
done

GFF_COUNT=$(wc -l < "$GFF_LIST")
echo "Total GFF files to cluster: $GFF_COUNT"

# ---- 4. Run Panaroo ----
echo "Running Panaroo core/accessory gene clustering..."
panaroo -i "$GFF_LIST" \
        -o "$OUTPUT/panaroo_out" \
        --clean-mode strict \
        -a core \
        --aligner mafft \
        -t "$THREADS"

# ---- 5. Run IQ-TREE core-genome phylogeny ----
echo "Running IQ-TREE core-genome phylogeny reconstruction..."
CORE_ALIGNMENT="$OUTPUT/panaroo_out/core_gene_alignment.aln"

if [[ ! -f "$CORE_ALIGNMENT" ]]; then
    echo "ERROR: Panaroo core alignment not found at $CORE_ALIGNMENT"
    exit 1
fi

iqtree -s "$CORE_ALIGNMENT" \
       -m "GTR+F+I+G4" \
       -bb 1000 \
       -nt "$THREADS" \
       -pre "$OUTPUT/iqtree_out/core_phylogeny" \
       -quiet

echo "03b_pangenome.sh completed successfully. Outputs saved in: $OUTPUT"
