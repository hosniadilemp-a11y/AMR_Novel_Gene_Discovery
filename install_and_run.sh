#!/bin/bash
# =======================================================================
#          OpenMM MD Simulation Automated Linux Launcher
# =======================================================================
set -e

echo "Checking for Conda installation..."

# Source conda environment scripts if conda is not in PATH
if ! command -v conda &> /dev/null; then
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        echo "[ERROR] Conda was not found in your PATH or common installation directories."
        echo "Please make sure conda is installed and initialized."
        exit 1
    fi
fi

# Ensure conda function is available
if ! declare -F conda > /dev/null; then
    eval "$(conda shell.bash hook)"
fi

# Check if openmm_env already exists
if conda env list | grep -q "openmm_env"; then
    echo "Environment 'openmm_env' already exists."
else
    echo "Creating environment 'openmm_env' (this may take a few minutes)..."
    conda create -n openmm_env -c conda-forge openmm openmmforcefields openff-toolkit rdkit python=3.11 -y
fi

echo "Activating environment 'openmm_env'..."
conda activate openmm_env

echo "Launching Molecular Dynamics simulation..."
python run_md.py

echo ""
echo "Simulation finished or stopped."
