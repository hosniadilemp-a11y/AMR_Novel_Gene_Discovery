#!/bin/bash
# =============================================================================
# 01_qc_trimming.sh — Quality Control & Read Trimming Pipeline
# =============================================================================
#
# USAGE:
#   bash scripts/01_qc_trimming.sh [OPTIONS]
#
# OPTIONS:
#   --r1 PATH          Forward reads FASTQ (required)
#   --r2 PATH          Reverse reads FASTQ (required)
#   --threads INT      Number of CPU threads (default: 4)
#   --output DIR       Output directory (default: results/step1_qc)
#   --check-only       Check if tools are available and exit
#   --help             Show this help message
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
R1=""
R2=""
THREADS=4
OUTPUT="results/step1_qc"
CHECK_ONLY=false

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)         R1="$2";         shift 2 ;;
        --r2)         R2="$2";         shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
        --output)     OUTPUT="$2";     shift 2 ;;
        --check-only) CHECK_ONLY=true; shift ;;
        --help)       head -20 "$0";   exit 0 ;;
        *)            echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---- Check dependencies ----
TOOLS=("fastqc" "cutadapt" "multiqc")
for tool in "${TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

if [[ "$CHECK_ONLY" == "true" ]]; then
    echo "All dependencies for 01_qc_trimming.sh are satisfied."
    exit 0
fi

# ---- Validate required inputs ----
if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: --r1 and --r2 are required."
    exit 1
fi

if [[ ! -f "$R1" ]]; then
    echo "ERROR: Forward read file '$R1' does not exist."
    exit 1
fi

if [[ ! -f "$R2" ]]; then
    echo "ERROR: Reverse read file '$R2' does not exist."
    exit 1
fi

# ---- Create output directories ----
echo "Creating output directories under: $OUTPUT"
mkdir -p "$OUTPUT/fastqc_raw"
mkdir -p "$OUTPUT/trimmed_reads"
mkdir -p "$OUTPUT/fastqc_trimmed"
mkdir -p "$OUTPUT/multiqc_report"

# ---- 1. FastQC on raw reads ----
echo "Running FastQC on raw reads..."
fastqc -o "$OUTPUT/fastqc_raw" -t "$THREADS" "$R1" "$R2"

# ---- 2. Cutadapt Trimming ----
echo "Running Cutadapt (Trimming adapters and low-quality bases Q30, minimum length 50)..."
cutadapt -q 30 -m 50 \
  -o "$OUTPUT/trimmed_reads/trimmed_R1.fastq" \
  -p "$OUTPUT/trimmed_reads/trimmed_R2.fastq" \
  "$R1" "$R2" > "$OUTPUT/cutadapt.log"

# ---- 3. FastQC on trimmed reads ----
echo "Running FastQC on trimmed reads..."
fastqc -o "$OUTPUT/fastqc_trimmed" -t "$THREADS" \
  "$OUTPUT/trimmed_reads/trimmed_R1.fastq" \
  "$OUTPUT/trimmed_reads/trimmed_R2.fastq"

# ---- 4. MultiQC Report ----
echo "Running MultiQC..."
multiqc "$OUTPUT" -o "$OUTPUT/multiqc_report"

echo "01_qc_trimming.sh completed successfully. Results saved in $OUTPUT"
