import os
import sys
import random
import time
import argparse

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
    sys.exit(1)

parser = argparse.ArgumentParser(description="GNAT Apo-State Replicate 2 MD Simulation")
parser.add_argument("--ns", type=float, default=50.0, help="Total simulation length in ns")
args = parser.parse_args()

# Configuration
PDB_PATH = "/media/adel/Data/Hosni/openmm_windows_setup/KNGPFPPJ_02769.pdb" 
solvated_pdb_path = "/media/adel/Data/Hosni/openmm_windows_setup/solvated_system.pdb"
TOTAL_NS = args.ns  # Simulation length
TIMESTEP_FS = 2.0  # 2 fs timestep
TOTAL_STEPS = int((TOTAL_NS * unit.nanoseconds) / (TIMESTEP_FS * unit.femtoseconds))
CHECKPOINT_INTERVAL = 500000  # Save checkpoint every 1 ns (500,000 steps)
TRAJECTORY_INTERVAL = 5000  # Save trajectory coordinates every 10 ps (5,000 steps)

checkpoint_file = "/media/adel/Data/Hosni/openmm_windows_setup/md_apo_r2_checkpoint.chk"
trajectory_file = "/media/adel/Data/Hosni/openmm_windows_setup/md_apo_r2_trajectory.dcd"
log_file = "/media/adel/Data/Hosni/openmm_windows_setup/md_apo_r2_log.csv"

print("==========================================================")
print(" Starting GNAT Apo-State Replicate 2 MD Simulation")
print(f" Target Length: {TOTAL_NS} ns ({TOTAL_STEPS} steps)")
print("==========================================================")

# 1. Setup System Generator
print("Initializing forcefield generators...")
system_generator = SystemGenerator(
    forcefields=['amber14-all.xml', 'amber14/tip3p.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=[]
)

# 2. Load Solvated Structure
if os.path.exists(solvated_pdb_path):
    print(f"Loading existing solvated topology and positions from {solvated_pdb_path}...")
    pdb = app.PDBFile(solvated_pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
else:
    print(f"[ERROR] Solvated structure {solvated_pdb_path} not found.")
    sys.exit(1)

# Create OpenMM System
print("Creating OpenMM system...")
system = system_generator.create_system(modeller.topology)

# Add NPT Barostat
system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300.0*unit.kelvin, 25))

# Setup Integrator
integrator = mm.LangevinMiddleIntegrator(300.0*unit.kelvin, 1.0/unit.picosecond, TIMESTEP_FS*unit.femtoseconds)

# Platform selection: OpenCL (GPU)
platform_name = 'CPU'
platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
if 'OpenCL' in platforms:
    platform_name = 'OpenCL'

print(f"Selecting hardware platform: {platform_name}")
platform = mm.Platform.getPlatformByName(platform_name)
properties = {'OpenCLPrecision': 'mixed'} if platform_name == 'OpenCL' else {}

# Create Simulation object
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# 3. Energy Minimization & Randomized Thermalization
if not os.path.exists(checkpoint_file):
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1500)
    
    # Thermalize velocities with a random seed
    seed = random.randint(1, 1000000)
    print(f"Thermalizing velocities to 300 K with random seed {seed}...")
    simulation.context.setVelocitiesToTemperature(300.0*unit.kelvin, seed)
    
    print("Running NVT thermalization (100 ps)...")
    simulation.step(50000)
    
    print("Running NPT equilibration (100 ps)...")
    simulation.step(50000)
    
    # Save initial checkpoint
    with open(checkpoint_file, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
    print("Equilibration complete. Replicate 2 checkpoint saved.")
else:
    # Resume from checkpoint
    print(f"Resuming simulation from checkpoint: {checkpoint_file}")
    with open(checkpoint_file, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

# 4. Production Run Setup
resuming = os.path.exists(trajectory_file)
simulation.reporters.append(app.DCDReporter(trajectory_file, TRAJECTORY_INTERVAL, append=resuming))
simulation.reporters.append(app.StateDataReporter(log_file, TRAJECTORY_INTERVAL, step=True,
                                                 potentialEnergy=True, kineticEnergy=True, 
                                                 totalEnergy=True, temperature=True, 
                                                 volume=True, speed=True, append=resuming))

# Run production steps
print(f"Running production MD simulation for Replicate 2...")
steps_completed = 0
if resuming and os.path.exists(log_file):
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                last_line = lines[-1].split(',')
                steps_completed = int(last_line[0])
                print(f"Resuming: {steps_completed:,} steps already completed.")
    except Exception as e:
        print(f"Error checking steps completed from log: {e}")

while steps_completed < TOTAL_STEPS:
    chunk_steps = min(CHECKPOINT_INTERVAL, TOTAL_STEPS - steps_completed)
    simulation.step(chunk_steps)
    steps_completed += chunk_steps
    
    # Save checkpoint
    with open(checkpoint_file, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
    pct = (steps_completed / TOTAL_STEPS) * 100.0
    print(f"Completed {steps_completed:,} / {TOTAL_STEPS:,} steps ({pct:.1f}%). Checkpoint updated.")

print("GNAT Apo Replicate 2 simulation completed successfully.")
