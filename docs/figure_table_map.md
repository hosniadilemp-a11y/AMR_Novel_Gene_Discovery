# Manuscript Figure and Table Mapping Registry

This document lists the files, generation scripts, and source paths for every figure and table in the manuscript and supplementary materials.

---

## 1. Main Text Figures & Tables

### Figure 1: Assembly QC and Coverage Depth
*   **Public Path:** `figures/figure1_assembly_coverage.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step2_assembly/coverage_validation/`

### Figure 2: Pangenome Partitioning
*   **Public Path:** `figures/figure2_pangenome_partition.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step3_pangenome/panaroo_out/`

### Figure 3: Phylogenomic Core Tree
*   **Public Path:** `figures/figure3_phylogeny_tree.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step3_pangenome/iqtree_out/`

### Figure 4: Insertion Sequences Distribution
*   **Public Path:** `figures/figure4_is_distribution.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step4_mge/isescan_out/`

### Figure 5: Feature Statistics Comparison
*   **Public Path:** `figures/figure5_stats_comparison.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step5_candidates/advanced_stats_results.json`

### Figure 6: AMR Profile Abundance Boxplots
*   **Public Path:** `figures/figure6_amr_comparison.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step3_annotation/amr_abricate/`

### Figure 7: CALIN Synteny Context (clinker)
*   **Public Path:** `figures/figure7_amr_genetic_context.pdf` (and `.png`)
*   **Generating Script:** `clinker` run on `results/step4_mge/clinker_out/`

### Figure 8: GC vs GC3 Codon Bias Scatter Plot
*   **Public Path:** `figures/figure8_amr_gc3_scatter.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step5_candidates/prioritized_candidates.tsv`

### Figure 9: PLM UMAP Embeddings Projection
*   **Public Path:** `figures/figure9_plm_umap_clustering.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `results/step7_plm/esm2_650m_coordinates.tsv`

### Figure 10: Molecular Dynamics Trajectory (RMSD & Rg)
*   **Public Path:** `figures/figure10_md_rmsd_rg.pdf` (and `.png`)
*   **Generating Script:** `scripts/generate_figures.py`
*   **Source Data:** `logs/step7_md_apo_gnat.csv`

---

## 2. Supplementary Tables

### Table S1: Assembly Quality Control Statistics
*   **Public Path:** `supplementary/Table_S1_Assembly_QC.tsv`
*   **Generating Script:** `scratch/create_tsv_tables.py`
*   **Source Data:** `results/step2_assembly/quast_report/`

### Table S2: Reference Cohort GenBank Accessions
*   **Public Path:** `supplementary/Table_S2_ST354_Cohort_Accessions.tsv`
*   **Generating Script:** `scratch/create_tsv_tables.py`
*   **Source Data:** `data/reference_accessions/`

### Table S6: Codon Bias Statistics and Mann-Whitney testing
*   **Public Path:** `supplementary/Table_S6_Codon_Bias_Statistics.tsv`
*   **Generating Script:** `scratch/create_tsv_tables.py`
*   **Source Data:** `results/step5_candidates/advanced_stats_results.json`

### Table S8: Foldseek Structural Alignment Validation
*   **Public Path:** `supplementary/Table_S8_Foldseek_Alignments.tsv`
*   **Generating Script:** `scratch/create_tsv_tables.py`
*   **Source Data:** `results/step6_structures/foldseek_results/`

### Table S11: PanGWAS Significant Orthologous Gene Hits
*   **Public Path:** `supplementary/Table_S11_PanGWAS_Results.tsv`
*   **Generating Script:** `scratch/create_tsv_tables.py`
*   **Source Data:** `supplementary/pan_gwas_association_results.csv`
