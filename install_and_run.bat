@echo off
echo =======================================================================
echo          OpenMM MD Simulation Automated Windows Launcher
echo =======================================================================
echo.
echo Checking for Conda installation...
where conda >nul 2>nul
if %errorlevel% neq 0 (
    echo [ERROR] Conda was not found in your PATH. 
    echo Please install Miniconda/Anaconda and make sure to run this script 
    echo from the "Anaconda Prompt".
    pause
    exit /b 1
)

echo Conda found. Creating environment 'openmm_env' (this may take a few minutes)...
call conda create -n openmm_env -c conda-forge openmm openmmforcefields openff-toolkit rdkit python=3.11 -y

echo Activating environment 'openmm_env'...
call conda activate openmm_env

echo Launching Molecular Dynamics simulation...
python run_md_windows.py

echo.
echo Simulation finished or stopped.
pause
