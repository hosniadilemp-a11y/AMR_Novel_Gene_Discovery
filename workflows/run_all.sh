#!/usr/bin/env bash
# =============================================================================
# run_all.sh — Master Pipeline Orchestrator
# AMR Novel Gene Discovery: Structure-guided prioritization of divergent
# virulence and resistance candidates in clinical E. coli ST354
# =============================================================================
#
# USAGE:
#   bash workflows/run_all.sh [OPTIONS]
#
# OPTIONS:
#   --r1 PATH          Forward reads FASTQ (required)
#   --r2 PATH          Reverse reads FASTQ (required)
#   --threads INT      CPU threads to use (default: 8)
#   --output DIR       Output root directory (default: results/)
#   --skip-md          Skip molecular dynamics steps (Steps 9-10)
#   --skip-gpu         Skip GPU steps (Steps 6-7, use precomputed results)
#   --from-step INT    Resume from a specific step number
#   --help             Show this help message
#
# EXAMPLES:
#   # Full pipeline (no GPU steps):
#   bash workflows/run_all.sh \
#       --r1 data/raw_reads/QA5221_R1.fastq.gz \
#       --r2 data/raw_reads/QA5221_R2.fastq.gz \
#       --threads 16 --skip-gpu --skip-md
#
#   # Resume from Step 5:
#   bash workflows/run_all.sh \
#       --r1 data/raw_reads/QA5221_R1.fastq.gz \
#       --r2 data/raw_reads/QA5221_R2.fastq.gz \
#       --from-step 5
#
# NOTE:
#   GPU steps (6 and 7) must be run separately on a GPU cluster or Kaggle.
#   Upload scripts/06a_esmfold_prediction.py and scripts/07_plm_embeddings.py
#   to a Kaggle notebook with a T4 GPU. See docs/pipeline_overview.md.
#
# ESTIMATED RUNTIMES (16 CPU cores, no GPU steps):
#   Step 1: ~1 h  |  Step 2: ~3 h  |  Step 3: ~2 h  |  Step 3b: ~12 h
#   Step 4: ~1 h  |  Step 5: ~30 min  |  Step 8: ~30 min
# =============================================================================

set -euo pipefail

# ---- Default parameters ----
R1=""
R2=""
THREADS=8
OUTPUT="results"
SKIP_MD=false
SKIP_GPU=false
FROM_STEP=1

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)       R1="$2";       shift 2 ;;
        --r2)       R2="$2";       shift 2 ;;
        --threads)  THREADS="$2";  shift 2 ;;
        --output)   OUTPUT="$2";   shift 2 ;;
        --skip-md)  SKIP_MD=true;  shift ;;
        --skip-gpu) SKIP_GPU=true; shift ;;
        --from-step) FROM_STEP="$2"; shift 2 ;;
        --help)     head -50 "$0"; exit 0 ;;
        *)          echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---- Validate required inputs ----
if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: --r1 and --r2 are required."
    echo "Run: bash workflows/run_all.sh --help"
    exit 1
fi

# ---- Activate conda environment ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate amr_env

# ---- Logging ----
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
echo "Pipeline started: $(date)" | tee -a "$MASTER_LOG"
echo "Resuming from Step $FROM_STEP" | tee -a "$MASTER_LOG"

run_step() {
    local STEP=$1; shift
    local DESC=$1; shift
    if [[ $STEP -ge $FROM_STEP ]]; then
        echo "============================================================" | tee -a "$MASTER_LOG"
        echo "  STEP $STEP: $DESC" | tee -a "$MASTER_LOG"
        echo "  Started: $(date)" | tee -a "$MASTER_LOG"
        echo "============================================================" | tee -a "$MASTER_LOG"
        "$@" 2>&1 | tee -a "$MASTER_LOG"
        echo "  Completed: $(date)" | tee -a "$MASTER_LOG"
    else
        echo "  Skipping Step $STEP (already completed)" | tee -a "$MASTER_LOG"
    fi
}

# =============================================================================
# STEP 1: Quality Control & Trimming
# =============================================================================
run_step 1 "QC & Trimming (FastQC, Cutadapt, MultiQC)" \
    bash scripts/01_qc_trimming.sh \
        --r1 "$R1" --r2 "$R2" \
        --threads "$THREADS" \
        --output "$OUTPUT/step1_qc"

# =============================================================================
# STEP 2: De Novo Assembly & Validation
# =============================================================================
run_step 2 "Assembly (SPAdes, QUAST, Minimap2, mosdepth)" \
    bash scripts/02_assembly.sh \
        --r1 "$OUTPUT/step1_qc/trimmed_reads/trimmed_R1.fastq" \
        --r2 "$OUTPUT/step1_qc/trimmed_reads/trimmed_R2.fastq" \
        --threads "$THREADS" \
        --output "$OUTPUT/step2_assembly"

# =============================================================================
# STEP 3: Genome Annotation & AMR Profiling
# =============================================================================
run_step 3 "Annotation (Prokka, ABRicate, AMRFinderPlus, MOB-suite)" \
    bash scripts/03_annotation.sh \
        --contigs "$OUTPUT/step2_assembly/spades_output/contigs.fasta" \
        --threads "$THREADS" \
        --output "$OUTPUT/step3_annotation"

# =============================================================================
# STEP 3b: Pangenome & Phylogeny
# =============================================================================
run_step 4 "Pangenome & Phylogeny (Panaroo, IQ-TREE)" \
    bash scripts/03b_pangenome.sh \
        --focal-gff "$OUTPUT/step3_annotation/prokka_out/QA5221.gff" \
        --accessions data/reference_accessions/st354_cohort_accessions.txt \
        --threads "$THREADS" \
        --output "$OUTPUT/step3_pangenome"

# =============================================================================
# STEP 4: Mobile Genetic Elements & Synteny
# =============================================================================
run_step 5 "MGE Detection (ISEScan, IntegronFinder, clinker)" \
    bash scripts/04_mge_detection.sh \
        --contigs "$OUTPUT/step2_assembly/spades_output/contigs.fasta" \
        --threads "$THREADS" \
        --output "$OUTPUT/step4_mge"

# =============================================================================
# STEP 5: Novel Candidate Extraction
# =============================================================================
run_step 6 "Novel Candidate Extraction (BLASTp, Pfam)" \
    python3 scripts/05_novel_candidates.py \
        --prokka-dir "$OUTPUT/step3_annotation/prokka_out" \
        --pangenome "$OUTPUT/step3_pangenome/panaroo_out/gene_presence_absence.csv" \
        --min-length 200 \
        --output "$OUTPUT/step5_candidates"

# =============================================================================
# STEP 6: Protein Structure Prediction (GPU — Kaggle)
# =============================================================================
if [[ "$SKIP_GPU" == "false" ]]; then
    echo "================================================================"
    echo "  STEP 6: Structure Prediction — REQUIRES GPU"
    echo "  Upload scripts/06a_esmfold_prediction.py to Kaggle T4 GPU."
    echo "  Download results to: $OUTPUT/step6_structures/esmfold/"
    echo "  Then run Foldseek and TM-align locally:"
    echo "================================================================"
    run_step 7 "Foldseek Search + TM-align Validation" \
        python3 scripts/06b_foldseek_search.py \
            --pdb-dir "$OUTPUT/step6_structures/esmfold" \
            --output "$OUTPUT/step6_structures/foldseek_results"
    run_step 7 "TM-align Local Validation" \
        python3 scripts/06c_tmalign_validation.py \
            --query-pdb "$OUTPUT/step6_structures/esmfold" \
            --target-pdb "$OUTPUT/step6_structures/target_structures" \
            --output "$OUTPUT/step6_structures/tmalign_results.tsv"
else
    echo "  Skipping GPU steps — using precomputed results from results/step6_structures/" | tee -a "$MASTER_LOG"
fi

# =============================================================================
# STEP 7: PLM Embeddings & Docking (GPU + CPU)
# =============================================================================
if [[ "$SKIP_GPU" == "false" ]]; then
    echo "  STEP 7: PLM Embeddings — REQUIRES GPU (Kaggle)"
    echo "  Upload scripts/07_plm_embeddings.py to Kaggle."
    echo "  Download results to: $OUTPUT/step7_plm/"
fi

run_step 8 "Molecular Docking + MM-GBSA (AutoDock Vina)" \
    python3 scripts/10_docking.py \
        --receptor "$OUTPUT/step6_structures/esmfold/KNGPFPPJ_02769.pdb" \
        --ligands data/ligands/ \
        --output "$OUTPUT/step7_docking"

# =============================================================================
# STEP 8: Novelty Scoring
# =============================================================================
run_step 9 "Novelty Scoring (16-point weighted framework)" \
    python3 scripts/08_novelty_scoring.py \
        --candidates "$OUTPUT/step5_candidates/prioritized_candidates.tsv" \
        --foldseek "$OUTPUT/step6_structures/foldseek_results" \
        --plm "$OUTPUT/step7_plm" \
        --mge "$OUTPUT/step4_mge" \
        --coverage "$OUTPUT/step2_assembly/coverage_validation" \
        --output "$OUTPUT/step8_scoring"

# =============================================================================
# STEP 9-10: Molecular Dynamics (Optional — requires GPU + days of runtime)
# =============================================================================
if [[ "$SKIP_MD" == "false" ]]; then
    echo "================================================================"
    echo "  STEP 9: MD Simulations — LONG RUNTIME (days per protein)"
    echo "  Proteins: GNAT_KA27 (127 ns), Ehly_61 (100 ns), OAgP_161 (100 ns)"
    echo "================================================================"
    run_step 10 "MD Simulation — GNAT_KA27 Apo" \
        python3 scripts/09b_md_production.py \
            --config config/md_config.yaml \
            --protein KNGPFPPJ_02769 \
            --state apo \
            --output "$OUTPUT/step9_md/gnat_apo"
    run_step 11 "MD Analysis — All Proteins" \
        python3 scripts/09c_md_analysis.py \
            --md-dir "$OUTPUT/step9_md" \
            --output "$OUTPUT/step9_md/analysis"
else
    echo "  MD steps skipped (--skip-md). Precomputed logs are in logs/." | tee -a "$MASTER_LOG"
fi

# =============================================================================
# STEP 11: Figure Generation
# =============================================================================
run_step 12 "Figure Generation (all manuscript figures)" \
    python3 scripts/generate_figures.py \
        --all \
        --results "$OUTPUT" \
        --output figures/

echo "============================================================" | tee -a "$MASTER_LOG"
echo "  PIPELINE COMPLETE: $(date)"                                  | tee -a "$MASTER_LOG"
echo "  Master log: $MASTER_LOG"                                     | tee -a "$MASTER_LOG"
echo "  Figures: figures/"                                           | tee -a "$MASTER_LOG"
echo "============================================================" | tee -a "$MASTER_LOG"
