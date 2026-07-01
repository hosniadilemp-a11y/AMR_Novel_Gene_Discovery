# Installation and Setup Guide

This document describes how to set up the Conda environment and dependencies required to run the pipeline.

## 1. Prerequisites

The pipeline is designed to run on a Linux environment (tested on Ubuntu 20.04/22.04 LTS).

*   **Miniconda** or **Anaconda** (Python 3.10+)
*   **CUDA Toolkit** (version 11.7+ or 12.x) to execute GPU-accelerated steps (ESMFold, ESM-2 PLM embeddings, and OpenMM simulations).

## 2. Conda Environment Setup

We provide two configuration files in the `environment/` directory:
*   `environment.yml` (recommended for conda-native installs)
*   `requirements.txt` (standard pip packages)

### Quick Setup

To create the main environment (`amr_env`) containing the bioinformatics tools (SPAdes, Prokka, Panaroo, etc.):

```bash
# Create the environment from yml file
conda env create -f environment/environment.yml

# Activate the environment
conda activate amr_env
```

### Molecular Dynamics Environment (`openmm_env`)

For molecular dynamics simulations, we use a separate environment with specific OpenMM versions:

```bash
# Create a dedicated OpenMM environment
conda create -n openmm_env -c conda-forge openmm openmmforcefields pdbfixer openff-toolkit rdkit python=3.10 -y

# Activate
conda activate openmm_env
```

## 3. Applying MOB-Suite Pandas Compatibility Patch

Due to deprecations in pandas 3.0+, we include a patch to fix compatibility issues inside `mob-suite`. Run the following command after activating `amr_env`:

```bash
# Run the environment repair script
bash scripts/patch_mob_suite.sh
```

This updates `mob_suite` to use `pd.concat` instead of the deprecated `df.append`, preventing runtime crashes during plasmid profiling (Step 3).

## 4. Verification

Verify that all key tools are accessible in your path:

```bash
spades.py --version
prokka --version
panaroo --version
integron_finder --version
```
