# Pipeline Execution Logs

This directory contains representative execution logs for each step of the structure-guided pangenome-to-structure prioritization pipeline:

* **[01_qc_cutadapt.log](01_qc_cutadapt.log)**: Raw read adapter trimming and quality filtering execution log.
* **[02_assembly_spades.log](02_assembly_spades.log)**: *De novo* genome assembly path selection and coverage stats.
* **[02_assembly_quast.log](02_assembly_quast.log)**: Assembly evaluation and quality metric report.
* **[03_annotation_prokka.log](03_annotation_prokka.log)**: Genome annotation and feature prediction summary.
* **[03_annotation_iqtree.log](03_annotation_iqtree.log)**: Core-genome Maximum-Likelihood phylogeny reconstruction run log.
* **[04_mge_clinker.log](04_mge_clinker.log)**: Mobilome genomic synteny mapping (clinker) run log.
* **[09_md_progress.log](09_md_progress.log)**: OpenMM molecular dynamics explicit-solvent production simulation trajectory log.
* **10_docking_*.log**: AutoDock Vina molecular docking log files for positive and negative validation panel:
  * **[10_docking_kanamycin.log](10_docking_kanamycin.log)** (Positive Control)
  * **[10_docking_gentamicin.log](10_docking_gentamicin.log)** (Positive Control)
  * **[10_docking_amikacin.log](10_docking_amikacin.log)** (Positive Control)
  * **[10_docking_penicilling.log](10_docking_penicilling.log)** (Negative Decoy)
  * **[10_docking_tetracycline.log](10_docking_tetracycline.log)** (Negative Decoy)
  * **[10_docking_dglucose.log](10_docking_dglucose.log)** (Negative Decoy)
