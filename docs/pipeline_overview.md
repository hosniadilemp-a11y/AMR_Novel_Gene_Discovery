# Pipeline Overview and Architecture

The AMR Novel Gene Discovery pipeline is organized into 12 structured steps, coordinated by the master script `workflows/run_all.sh`. 

---

## Pipeline Workflow Steps

### Step 1: Quality Control & Trimming
*   **Script:** `scripts/01_qc_trimming.sh`
*   **Tools:** FastQC, Cutadapt, MultiQC
*   **Purpose:** Trims low-quality bases and adapter sequences from raw paired-end reads.

### Step 2: De Novo Assembly & Validation
*   **Script:** `scripts/02_assembly.sh`
*   **Tools:** SPAdes, QUAST, minimap2, mosdepth
*   **Purpose:** Assembles reads into contigs, evaluates N50/L50 assembly statistics, and computes read coverage depth mapping.

### Step 3: Genome Annotation & AMR Profiling
*   **Script:** `scripts/03_annotation.sh`
*   **Tools:** Prokka, ABRicate (ResFinder & VFDB), AMRFinderPlus, MOB-recon
*   **Purpose:** Annotates the assembly, profiles known acquired resistance genes, virulence factors, and plasmids.

### Step 3b: Reference Cohort Download & Pangenome Analysis
*   **Script:** `scripts/03b_pangenome.sh`
*   **Tools:** Panaroo, IQ-TREE
*   **Purpose:** Downloads and annotates the 32 reference genomes, performs orthology clustering to identify core/accessory gene frequencies, and reconstructs a maximum-likelihood phylogeny.

### Step 4: Mobile Genetic Elements Detection
*   **Script:** `scripts/04_mge_detection.sh`
*   **Tools:** ISEScan, IntegronFinder, clinker
*   **Purpose:** Screen for Insertion Sequences, complete integron cassettes, and map syntenic contexts of candidates.

### Step 5: Novel Candidate Extraction
*   **Script:** `scripts/05_novel_candidates.py`
*   **Purpose:** Filters genes of length $\ge 200$ aa, screens out any proteins matching Pfam domains or known CARD entries, and identifies singletons/rare accessory genes.

### Step 6: Protein Structure Prediction & Alignments
*   **Scripts:** 
    *   `scripts/06a_esmfold_prediction.py` (Run on GPU cluster / Kaggle)
    *   `scripts/06b_foldseek_search.py` (Local structure search)
    *   `scripts/06c_tmalign_validation.py` (Local structural alignment)
*   **Purpose:** Predict 3D structures for candidates using ESMFold. Conduct search against PDB database and align candidates to templates using TM-align.

### Step 7: Protein Language Model (PLM) Embeddings
*   **Script:** `scripts/07_plm_embeddings.py` (Run on GPU cluster / Kaggle)
*   **Purpose:** Generate 1280-dimension protein embeddings using ESM-2 (650M) and project them using UMAP/t-SNE to search for cluster membership.

### Step 8: Molecular Docking & MM-GBSA
*   **Script:** `scripts/10_docking.py`
*   **Tools:** AutoDock Vina
*   **Purpose:** Perform receptor-ligand docking against priority candidates and calculate binding free energies.

### Step 9: Novelty Scoring Framework
*   **Script:** `scripts/08_novelty_scoring.py`
*   **Purpose:** Computes weighted novelty score using the 16-point scale combining structure, genomic, and sequence features.

### Steps 10-11: Molecular Dynamics Simulations
*   **Scripts:**
    *   `scripts/09b_md_production.py` (GPU script to execute long MD)
    *   `scripts/run_md_replicate.py` (Execute seed replicates)
    *   `scripts/09c_md_analysis.py` (Parse RMSD, RMSF, Rg)
*   **Purpose:** Verify structural stability of prioritized candidates.

### Step 12: Figure Generation
*   **Script:** `scripts/generate_figures.py`
*   **Purpose:** Automatically compiles all results and generates plots.
