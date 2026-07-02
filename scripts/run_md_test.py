# ==============================================================================
# GPU-ACCELERATED OPENMM MD PIPELINE (TEST RUN)
# ==============================================================================
import os
import sys

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmmforcefields.generators import SystemGenerator
    from openff.toolkit.topology import Molecule
    from rdkit import Chem
    from pdbfixer import PDBFixer
except ImportError:
    print("[ERROR] OpenMM or its dependencies are not installed in the active environment.")
    print("Please make sure you activated the 'openmm_env' environment.")
    sys.exit(1)

# Configuration for short verification
PDB_PATH = "KNGPFPPJ_02769.pdb" 
STATE = "apo" # Options: 'apo', 'binary', 'ternary'
TOTAL_NS = 0.0002 # 0.2 ps (100 steps)
TIMESTEP_FS = 2.0 # 2 fs timestep
TOTAL_STEPS = int((TOTAL_NS * unit.nanoseconds) / (TIMESTEP_FS * unit.femtoseconds))
CHECKPOINT_INTERVAL = 50 # Save checkpoint every 50 steps
TRAJECTORY_INTERVAL = 10 # Save trajectory coordinates every 10 steps

output_prefix = f"md_{STATE}_test"
checkpoint_file = f"{output_prefix}_checkpoint.chk"
trajectory_file = f"{output_prefix}_trajectory.dcd"
log_file = f"{output_prefix}_log.csv"

# 1. Setup System Generator (Amber14 for protein, TIP3P for water, GAFF2 for ligands)
print("Initializing forcefield generators...")
ligands = []
system_generator = SystemGenerator(
    forcefields=['amber14-all.xml', 'amber14/tip3p.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=ligands
)

# 2. Prepare and Solvate Structure using PDBFixer
print("Preparing structure with PDBFixer...")
if not os.path.exists(PDB_PATH):
    print(f"[ERROR] PDB file '{PDB_PATH}' not found in the current directory.")
    sys.exit(1)

fixer = PDBFixer(filename=PDB_PATH)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)  # pH 7.0

modeller = app.Modeller(fixer.topology, fixer.positions)

# Add water and neutralizing ions (0.15 M NaCl, 1.0 nm padding)
print("Solvating system...")
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

# Create OpenMM System
print("Creating OpenMM system...")
system = system_generator.create_system(modeller.topology)

# Add NPT Barostat for constant pressure (1 atm)
system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300.0*unit.kelvin, 25))

# Setup Integrator (Langevin Middle Integrator for NVT/NPT)
integrator = mm.LangevinMiddleIntegrator(300.0*unit.kelvin, 1.0/unit.picosecond, TIMESTEP_FS*unit.femtoseconds)

# Platform selection: Use OpenCL (GPU) - CUDA PTX incompatible with RTX 5070 Ti + CUDA 13.2 driver
# OpenCL is fully supported and GPU-accelerated on NVIDIA hardware
platform_name = 'CPU'
platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print(f"Available platforms: {platforms}")

if 'OpenCL' in platforms:
    platform_name = 'OpenCL'

print(f"Selecting hardware platform: {platform_name}")
platform = mm.Platform.getPlatformByName(platform_name)

properties = {}
if platform_name == 'OpenCL':
    properties = {'OpenCLPrecision': 'mixed'}

# Create Simulation object
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# 3. Energy Minimization & Equilibration (NVT/NPT)
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=100)

print("Running short NVT thermalization (100 steps)...")
simulation.step(100)

print("Running short NPT equilibration (100 steps)...")
simulation.step(100)

# Save initial checkpoint
with open(checkpoint_file, 'wb') as f:
    f.write(simulation.context.createCheckpoint())

# 4. Production Run Setup
simulation.reporters.append(app.DCDReporter(trajectory_file, TRAJECTORY_INTERVAL, append=False))
simulation.reporters.append(app.StateDataReporter(log_file, TRAJECTORY_INTERVAL, step=True,
                                                 potentialEnergy=True, kineticEnergy=True, 
                                                 totalEnergy=True, temperature=True, 
                                                 volume=True, speed=True, append=False))

# Run production steps
print(f"Running production MD simulation for {TOTAL_NS} ns ({TOTAL_STEPS} steps)...")
simulation.step(TOTAL_STEPS)

# Save final checkpoint
with open(checkpoint_file, 'wb') as f:
    f.write(simulation.context.createCheckpoint())
print("Simulation completed successfully.")
