================================================================================
          MOLECULAR DYNAMICS SIMULATION SETUP FOR WINDOWS (OPENMM)
================================================================================

This folder contains everything you need to run GPU-accelerated molecular dynamics (MD)
simulations for the candidate KNGPFPPJ_02769 on Windows.

--------------------------------------------------------------------------------
1. CONDA vs PYTHON VENV (RECOMMENDATION)
--------------------------------------------------------------------------------
We STRONGLY recommend using CONDA (via Miniconda or Anaconda) instead of Python venv.

* Why Conda?
  OpenMM contains compiled C++ source code and binds directly to graphics driver
  libraries (CUDA and OpenCL). Building and linking these libraries natively on
  Windows inside a standard Python venv requires Visual Studio C++ build tools, 
  manual CUDA SDK linking, and can lead to DLL loading errors.
  Conda solves this by distributing pre-compiled binaries and automatically 
  linking all GPU dependencies (cudatoolkit, ocl-icd, etc.) in a single step.

--------------------------------------------------------------------------------
2. QUICK START INSTRUCTIONS
--------------------------------------------------------------------------------
To run the setup and launch the simulation:

Step 1: Install Miniconda
  If you do not have Conda installed, download and run the Windows 64-bit installer:
  https://docs.anaconda.com/miniconda/

Step 2: Run the automated batch script
  Double-click the file:
      install_and_run.bat
  This script will automatically open the command prompt, create the 'openmm_env'
  environment, activate it, install all dependencies, and run the MD simulation.

--------------------------------------------------------------------------------
3. RUNNING MANUALLY (ALTERNATIVE)
--------------------------------------------------------------------------------
If you prefer to run the commands yourself, open the "Anaconda Prompt" from your
Windows Start Menu, navigate to this folder, and run:

  conda create -n openmm_env -c conda-forge openmm openmmforcefields openff-toolkit rdkit python=3.11 -y
  conda activate openmm_env
  python run_md_windows.py

--------------------------------------------------------------------------------
4. FILES IN THIS FOLDER
--------------------------------------------------------------------------------
* KNGPFPPJ_02769.pdb  : The predicted 3D structure coordinate file.
* run_md_windows.py   : Standalone OpenMM simulation Python script.
* install_and_run.bat : Automated environment creator and execution script.
* README_windows.txt  : This instruction guide.
================================================================================
