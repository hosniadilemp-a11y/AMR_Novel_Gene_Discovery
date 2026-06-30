#!/bin/bash
# =============================================================================
# 02_assembly.sh — De Novo Assembly & Read-Mapping Validation
# =============================================================================
#
# USAGE:
#   bash scripts/02_assembly.sh [OPTIONS]
#
# OPTIONS:
#   --r1 PATH          Forward trimmed reads (required)
#   --r2 PATH          Reverse trimmed reads (required)
#   --threads INT      Number of CPU threads (default: 8)
#   --output DIR       Output directory (default: results/step2_assembly)
#   --reference PATH   Reference genome FASTA for validation (optional)
#   --help             Show this help message
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
R1=""
R2=""
THREADS=8
OUTPUT="results/step2_assembly"
REF_GENOME=""

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)        R1="$2";        shift 2 ;;
        --r2)        R2="$2";        shift 2 ;;
        --threads)   THREADS="$2";   shift 2 ;;
        --output)    OUTPUT="$2";    shift 2 ;;
        --reference) REF_GENOME="$2"; shift 2 ;;
        --help)      head -20 "$0";  exit 0 ;;
        *)           echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---- Check dependencies ----
TOOLS=("spades.py" "quast.py" "minimap2" "samtools" "mosdepth" "curl")
for tool in "${TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

# ---- Validate required inputs ----
if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: --r1 and --r2 are required."
    exit 1
fi

if [[ ! -f "$R1" ]]; then
    echo "ERROR: Trimmed forward read file '$R1' does not exist."
    exit 1
fi

if [[ ! -f "$R2" ]]; then
    echo "ERROR: Trimmed reverse read file '$R2' does not exist."
    exit 1
fi

# ---- Create output directories ----
mkdir -p "$OUTPUT/reference"
mkdir -p "$OUTPUT/spades_output"
mkdir -p "$OUTPUT/quast_report"
mkdir -p "$OUTPUT/mapping"
mkdir -p "$OUTPUT/coverage_validation"

# ---- 1. Download or use Reference Genome ----
if [[ -z "$REF_GENOME" ]]; then
    echo "Downloading E. coli K-12 Reference Genome (NC_000913.3)..."
    REF_PATH="$OUTPUT/reference/ecoli_k12_ref.fasta"
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" > "$REF_PATH"
else
    REF_PATH="$REF_GENOME"
fi

# ---- 2. Run SPAdes Assembly ----
echo "Running SPAdes De Novo Assembly with careful flag..."
spades.py -1 "$R1" -2 "$R2" -o "$OUTPUT/spades_output" --careful -t "$THREADS"

# Verify assembly output
if [[ ! -f "$OUTPUT/spades_output/contigs.fasta" ]]; then
    echo "ERROR: SPAdes assembly failed! contigs.fasta not generated."
    exit 1
fi

# ---- 3. Run QUAST for Evaluation ----
echo "Running QUAST assembly evaluation..."
quast.py "$OUTPUT/spades_output/contigs.fasta" \
    -r "$REF_PATH" \
    -o "$OUTPUT/quast_report"

# ---- 4. Map Contigs to Reference to find NOVEL elements ----
echo "Mapping contigs against reference using minimap2..."
minimap2 -ax asm5 "$REF_PATH" \
    "$OUTPUT/spades_output/contigs.fasta" > "$OUTPUT/mapping/aln.sam"

samtools view -S -b "$OUTPUT/mapping/aln.sam" > "$OUTPUT/mapping/aln_unsorted.bam"
samtools sort -o "$OUTPUT/mapping/aln_sorted.bam" "$OUTPUT/mapping/aln_unsorted.bam"
samtools index "$OUTPUT/mapping/aln_sorted.bam"

# Extract unmapped contigs
echo "Extracting unmapped (potential novel) contigs..."
samtools view -b -f 4 "$OUTPUT/mapping/aln_sorted.bam" > "$OUTPUT/mapping/unmapped.bam"

# ---- 5. Read-Mapping Validation for Coverage ----
echo "Mapping raw reads back to assembly for coverage validation..."
minimap2 -ax sr "$OUTPUT/spades_output/contigs.fasta" "$R1" "$R2" > "$OUTPUT/coverage_validation/reads_to_contigs.sam"

samtools view -S -b "$OUTPUT/coverage_validation/reads_to_contigs.sam" > "$OUTPUT/coverage_validation/reads_unsorted.bam"
samtools sort -o "$OUTPUT/coverage_validation/reads_to_contigs_sorted.bam" "$OUTPUT/coverage_validation/reads_unsorted.bam"
samtools index "$OUTPUT/coverage_validation/reads_to_contigs_sorted.bam"

echo "Calculating coverage depth using mosdepth..."
mosdepth -n "$OUTPUT/coverage_validation/coverage_report" "$OUTPUT/coverage_validation/reads_to_contigs_sorted.bam"

echo "02_assembly.sh completed successfully. Output saved in: $OUTPUT"
