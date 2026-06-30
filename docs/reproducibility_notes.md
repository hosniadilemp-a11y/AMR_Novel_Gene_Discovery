# Technical and Reproducibility Notes

This document logs critical workarounds, bug resolutions, and patch mappings implemented to ensure full pipeline stability and identical reproducibility of results.

## 1. IntegronFinder Multi-Fasta Bug

### Issue
When run directly on a multi-fasta assembly containing many contigs, `integron_finder` may cross-contaminate database alignments between different contigs, yielding false-positive "complete" integron predictions.

### Solution
To prevent cross-contaminant leakage, we split `contigs.fasta` into separate sequence files (one per contig) and run `integron_finder` on each contig in isolation. The results are then merged back together.
This workaround is implemented in:
* `scripts/run_integron_finder_clean.py`
* Integrated automatically into `scripts/04_mge_detection.sh`

## 2. MOB-Suite Compatibility Patches

### Issue
The latest `pandas` versions (3.0+) deprecate and modify specific DataFrame serialization methods used in MOB-suite (e.g., `mob_recon`), causing traceback failures on execution.

### Solution
We patch the installed `mob-suite` Python files locally during environment setup to ensure compatibility:
1. Replace deprecated `pandas.DataFrame.append` with `pandas.concat`.
2. Clean up deprecated index checks in plasmid reconstruction scripts.
The patches are documented in `environment/INSTALL.md` and are automatically applied during environment initialization.
