#!/usr/bin/env python3
# =============================================================================
# 09b_md_production.py — OpenMM Molecular Dynamics Production Run
# =============================================================================
#
# USAGE:
#   python3 scripts/09b_md_production.py [OPTIONS]
#
# =============================================================================
import os
import sys
import argparse
import openmm as mm
import openmm.app as app
import openmm.unit as unit

def main():
    parser = argparse.ArgumentParser(description="Run OpenMM MD production simulation")
    parser.add_argument("--system", required=True, help="Path to solvated system PDB")
    parser.add_argument("--steps", type=int, default=50000000, help="Total simulation steps")
    parser.add_argument("--temp", type=float, default=310.0, help="Simulation temperature in Kelvin")
    parser.add_argument("--output", default="results/step9_md", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    basename = os.path.basename(args.system)
    prefix = basename.replace("_solvated.pdb", "")

    checkpoint_file = os.path.join(args.output, f"{prefix}_checkpoint.chk")
    trajectory_file = os.path.join(args.output, f"{prefix}_trajectory.dcd")
    log_file = os.path.join(args.output, f"{prefix}_log.csv")

    print(f"Loading solvated system PDB: {args.system}...")
    pdb = app.PDBFile(args.system)

    print("Creating OpenMM system...")
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds
    )

    # Add NPT Barostat
    system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, args.temp*unit.kelvin, 25))

    # Integrator
    integrator = mm.LangevinMiddleIntegrator(args.temp*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)

    # Platform Selection
    platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    print(f"Available hardware platforms: {platforms}")
    platform_name = 'CPU'
    if 'CUDA' in platforms:
        platform_name = 'CUDA'
    elif 'OpenCL' in platforms:
        platform_name = 'OpenCL'
        
    print(f"Selected platform: {platform_name}")
    platform = mm.Platform.getPlatformByName(platform_name)
    properties = {'OpenCLPrecision': 'mixed'} if platform_name == 'OpenCL' else {}

    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)

    # Checkpoint restore or Minimize/Equilibrate
    if not os.path.exists(checkpoint_file):
        print("Running energy minimization...")
        simulation.minimizeEnergy(maxIterations=1000)
        
        print("Running NVT equilibration (100 ps)...")
        simulation.step(50000)
        
        print("Running NPT equilibration (100 ps)...")
        simulation.step(50000)
        
        with open(checkpoint_file, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
        print(f"Initial checkpoint saved to {checkpoint_file}")
    else:
        print(f"Resuming from checkpoint {checkpoint_file}...")
        with open(checkpoint_file, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())

    # Reporters
    resuming = os.path.exists(trajectory_file)
    simulation.reporters.append(app.DCDReporter(trajectory_file, 5000, append=resuming))
    simulation.reporters.append(app.StateDataReporter(
        log_file, 5000, step=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, speed=True, append=resuming
    ))

    # Production loop
    print(f"Starting production simulation for {args.steps} steps...")
    steps_completed = 0
    checkpoint_interval = 500000
    while steps_completed < args.steps:
        chunk = min(checkpoint_interval, args.steps - steps_completed)
        simulation.step(chunk)
        steps_completed += chunk
        with open(checkpoint_file, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
        print(f"Completed {steps_completed} steps. Checkpoint saved.")

    print("MD simulation complete.")

if __name__ == "__main__":
    main()
