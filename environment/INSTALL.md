# Installation Guide

This document provides step-by-step instructions for setting up the computational genomics and molecular dynamics environment required to run the pipeline.

## 1. Conda Environment Setup

We recommend using Miniconda to manage bioinformatics dependencies.

### Step 1.1: Install Miniconda

If you do not have Conda installed, download and run the installer:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### Step 1.2: Import Environment

Create the environment using the provided `environment.yml` configuration:

```bash
conda env create -f environment.yml
conda activate amr_env
```

This installation will compile and link all compiled genomics packages (SPAdes, Prokka, Panaroo, etc.).

## 2. Initialize Bioinformatics Databases

Two primary tools require database setup prior to execution:

### Step 2.1: AMRFinderPlus Database

```bash
amrfinder -u
```

### Step 2.2: MOB-suite Database

```bash
mob_init
```

## 3. GPU Platform Support (OpenMM)

For explicit-solvent molecular dynamics simulations (Step 9):

1. Install OpenMM using Conda (automatically done inside `environment.yml`).
2. Verify GPU/CUDA platform recognition:

```bash
python3 -c "import openmm; print(openmm.Platform.getPlatform(0).getName())"
```

Expected output: `CUDA` or `OpenCL`.
