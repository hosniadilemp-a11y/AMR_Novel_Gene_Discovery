# ==============================================================================
# GPU-ACCELERATED OPENMM MD PIPELINE
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

# Configuration
PDB_PATH = "KNGPFPPJ_02769.pdb" 
STATE = "apo" # Options: 'apo', 'binary' (with Acetyl-CoA), 'ternary' (+ kanamycin)
TOTAL_NS = 100.0 # Total simulation time in nanoseconds
TIMESTEP_FS = 2.0 # 2 fs timestep
TOTAL_STEPS = int((TOTAL_NS * unit.nanoseconds) / (TIMESTEP_FS * unit.femtoseconds))
CHECKPOINT_INTERVAL = 500000 # Save checkpoint every 1 ns (500,000 steps)
TRAJECTORY_INTERVAL = 5000 # Save trajectory coordinates every 10 ps (5,000 steps)

output_prefix = f"md_{STATE}"
checkpoint_file = f"{output_prefix}_checkpoint.chk"
trajectory_file = f"{output_prefix}_trajectory.dcd"
log_file = f"{output_prefix}_log.csv"

# 1. Setup System Generator (Amber14 for protein, TIP3P for water, GAFF2 for ligands)
print("Initializing forcefield generators...")
ligands = []
if STATE in ['binary', 'ternary']:
    # Load Acetyl-CoA (co-factor) structure from smiles
    accoa_smiles = "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"
    rdkit_accoa = Chem.MolFromSmiles(accoa_smiles)
    accoa_mol = Molecule.from_rdkit(rdkit_accoa)
    ligands.append(accoa_mol)
    
if STATE == 'ternary':
    # Load Kanamycin substrate SMILES
    kanamycin_smiles = "C1C(C(C(C(O2)CO)O)O)N)O)OC3C(C(C(C(O3)CN)O)O)O)N"
    rdkit_kan = Chem.MolFromSmiles(kanamycin_smiles)
    kan_mol = Molecule.from_rdkit(rdkit_kan)
    ligands.append(kan_mol)

system_generator = SystemGenerator(
    forcefields=['amber14-all.xml', 'amber14/tip3p.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=ligands
)

# 2. Prepare and Solvate Structure using PDBFixer or Load Existing Solvated System
solvated_pdb_path = "solvated_system.pdb"

if os.path.exists(solvated_pdb_path):
    print(f"Loading existing solvated topology and positions from {solvated_pdb_path}...")
    pdb = app.PDBFile(solvated_pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
else:
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
    
    # Save solvated topology and positions so subsequent runs are deterministic and matching
    print(f"Saving solvated topology to {solvated_pdb_path}...")
    with open(solvated_pdb_path, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

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
if not os.path.exists(checkpoint_file):
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    print("Running NVT thermalization (100 ps)...")
    simulation.step(50000) # 100 ps NVT
    
    print("Running NPT equilibration (100 ps)...")
    simulation.step(50000) # 100 ps NPT
    
    # Save initial checkpoint
    with open(checkpoint_file, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
else:
    # Resume from checkpoint
    print(f"Resuming simulation from checkpoint: {checkpoint_file}")
    with open(checkpoint_file, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

# 4. Production Run Setup
# Use append=True only when resuming from an existing trajectory file
resuming = os.path.exists(trajectory_file)
simulation.reporters.append(app.DCDReporter(trajectory_file, TRAJECTORY_INTERVAL, append=resuming))
simulation.reporters.append(app.StateDataReporter(log_file, TRAJECTORY_INTERVAL, step=True,
                                                 potentialEnergy=True, kineticEnergy=True, 
                                                 totalEnergy=True, temperature=True, 
                                                 volume=True, speed=True, append=resuming))

# Run production steps
print(f"Running production MD simulation for {TOTAL_NS} ns ({TOTAL_STEPS} steps)...")
steps_completed = 0
while steps_completed < TOTAL_STEPS:
    chunk_steps = min(CHECKPOINT_INTERVAL, TOTAL_STEPS - steps_completed)
    simulation.step(chunk_steps)
    steps_completed += chunk_steps
    
    # Save checkpoint
    with open(checkpoint_file, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
    print(f"Completed {steps_completed} steps. Checkpoint saved.")

print("Simulation completed successfully.")
