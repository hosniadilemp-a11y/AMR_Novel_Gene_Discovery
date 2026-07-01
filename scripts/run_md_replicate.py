#!/usr/bin/env python3
import os
import sys
import argparse
import time
import datetime
import subprocess

# Configure OpenMM and dependencies
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmmforcefields.generators import SystemGenerator
    from openff.toolkit.topology import Molecule
    from rdkit import Chem
    from pdbfixer import PDBFixer
except ImportError as e:
    print(f"[ERROR] Required packages not found: {e}", file=sys.stderr)
    print("Please ensure the 'openmm_env' environment is activated.", file=sys.stderr)
    sys.exit(1)

def get_gpu_utilization():
    try:
        cmd = ["nvidia-smi", "--query-gpu=utilization.gpu,temperature.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            parts = [x.strip() for x in res.stdout.strip().split(',')]
            if len(parts) >= 4:
                return f"GPU Util: {parts[0]}%, Temp: {parts[1]}C, VRAM: {parts[2]}/{parts[3]} MB"
    except Exception:
        pass
    return "GPU Util: N/A"

def main():
    parser = argparse.ArgumentParser(description="GPU-Accelerated Replicate OpenMM MD Simulation Script")
    parser.add_argument("--pdb", required=True, help="Input ESMFold PDB file path")
    parser.add_argument("--outdir", required=True, help="Output directory for MD files")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for velocity and integrator initialization")
    parser.add_argument("--progress-file", required=True, help="Path to progress log file")
    parser.add_argument("--ns", type=float, default=20.0, help="Total simulation time in ns")
    args = parser.parse_args()
    
    # 1. Initialize Directories and Progress File
    os.makedirs(args.outdir, exist_ok=True)
    pdb_basename = os.path.basename(args.pdb)
    
    # Configure output files with seed markers
    checkpoint_file = os.path.join(args.outdir, f"md_replicate_seed_{args.seed}_checkpoint.chk")
    trajectory_file = os.path.join(args.outdir, f"md_replicate_seed_{args.seed}_trajectory.dcd")
    log_file = os.path.join(args.outdir, f"md_replicate_seed_{args.seed}_log.csv")
    solvated_pdb_path = os.path.join(args.outdir, "solvated_system.pdb")
    
    is_resume = os.path.exists(checkpoint_file)
    
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    mode = "a" if os.path.exists(args.progress_file) else "w"
    
    with open(args.progress_file, mode) as prog_f:
        if not is_resume:
            prog_f.write(f"==========================================================\n")
            prog_f.write(f"MD SIMULATION REPLICATE RUN FOR {pdb_basename} (SEED: {args.seed})\n")
            prog_f.write(f"==========================================================\n")
            prog_f.write(f"Start Timestamp       : {timestamp}\n")
            prog_f.write(f"Input ESMFold Structure: {args.pdb}\n")
            prog_f.write(f"Output Directory      : {args.outdir}\n")
            prog_f.write(f"Random Seed           : {args.seed}\n")
            prog_f.write(f"Target Simulation Time: {args.ns} ns\n")
            prog_f.write(f"Time Step             : 2.0 fs\n")
            prog_f.write(f"Force Field           : Amber14-all + TIP3P water\n")
            prog_f.write(f"Ion Concentration     : 0.15 M NaCl\n")
            prog_f.write(f"Water Padding Size    : 1.0 nm\n")
            prog_f.write(f"GPU Status Check      : {get_gpu_utilization()}\n")
            prog_f.write(f"----------------------------------------------------------\n")
            prog_f.write(f"System Preparation starting...\n")
            prog_f.flush()
            
    # 2. System Preparation & Solvation
    system_generator = SystemGenerator(
        forcefields=['amber14-all.xml', 'amber14/tip3p.xml'],
        small_molecule_forcefield='gaff-2.11',
        molecules=[]
    )
    
    if os.path.exists(solvated_pdb_path):
        print(f"Loading existing solvated topology and positions from {solvated_pdb_path}...")
        pdb = app.PDBFile(solvated_pdb_path)
        modeller = app.Modeller(pdb.topology, pdb.positions)
    else:
        print("Preparing structure with PDBFixer...")
        if not os.path.exists(args.pdb):
            print(f"[ERROR] PDB file '{args.pdb}' not found.", file=sys.stderr)
            sys.exit(1)
            
        fixer = PDBFixer(filename=args.pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)  # pH 7.0
        
        modeller = app.Modeller(fixer.topology, fixer.positions)
        
        # Add water and neutralizing ions (0.15 M NaCl, 1.0 nm padding)
        print("Solvating system...")
        modeller.addSolvent(
            system_generator.forcefield, 
            model='tip3p', 
            padding=1.0*unit.nanometer, 
            ionicStrength=0.15*unit.molar
        )
        
        print(f"Saving solvated topology to {solvated_pdb_path}...")
        with open(solvated_pdb_path, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Solvation completed. Atoms: {len(modeller.positions)}\n")
            prog_f.flush()

    # Create OpenMM System
    print("Creating OpenMM system...")
    system = system_generator.create_system(modeller.topology)
    system.addForce(mm.MonteCarloBarostat(1.0*unit.atmosphere, 300.0*unit.kelvin, 25))
    
    # Setup Integrator with the replicate seed
    TIMESTEP_FS = 2.0
    TOTAL_STEPS = int((args.ns * unit.nanoseconds) / (TIMESTEP_FS * unit.femtoseconds))
    CHECKPOINT_INTERVAL = 500000 # 1 ns
    TRAJECTORY_INTERVAL = 5000 # 10 ps
    
    integrator = mm.LangevinMiddleIntegrator(300.0*unit.kelvin, 1.0/unit.picosecond, TIMESTEP_FS*unit.femtoseconds)
    integrator.setRandomNumberSeed(args.seed)  # Integrator random seed
    
    # Platform selection
    platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    platform_name = 'CPU'
    if 'OpenCL' in platforms:
        platform_name = 'OpenCL'
    elif 'CUDA' in platforms:
        platform_name = 'CUDA'
        
    print(f"Selecting platform: {platform_name}")
    platform = mm.Platform.getPlatformByName(platform_name)
    properties = {'OpenCLPrecision': 'mixed'} if platform_name == 'OpenCL' else {}
    
    # Create Simulation object
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    
    # 3. Minimization and Equilibration
    if not is_resume:
        print("Minimizing energy...")
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Energy Minimization starting...\n")
            prog_f.flush()
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Initialize velocities with seed
        print(f"Initializing velocities at 300 K using seed {args.seed}...")
        simulation.context.setVelocitiesToTemperature(300.0*unit.kelvin, args.seed)
        
        print("Running NVT thermalization (100 ps)...")
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] NVT Equilibration starting (100 ps)...\n")
            prog_f.flush()
        simulation.step(50000) # 100 ps NVT
        
        print("Running NPT equilibration (100 ps)...")
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] NPT Equilibration starting (100 ps)...\n")
            prog_f.flush()
        simulation.step(50000) # 100 ps NPT
        
        # Save initial checkpoint
        with open(checkpoint_file, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
            
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Equilibration completed. Checkpoint saved.\n")
            prog_f.write(f"----------------------------------------------------------\n")
            prog_f.flush()
    else:
        print(f"Resuming from checkpoint: {checkpoint_file}")
        with open(checkpoint_file, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Resuming simulation from checkpoint.\n")
            prog_f.flush()
            
    # 4. Production Run Setup
    resuming = os.path.exists(trajectory_file)
    simulation.reporters.append(app.DCDReporter(trajectory_file, TRAJECTORY_INTERVAL, append=resuming))
    simulation.reporters.append(app.StateDataReporter(log_file, TRAJECTORY_INTERVAL, step=True,
                                                     potentialEnergy=True, kineticEnergy=True, 
                                                     totalEnergy=True, temperature=True, 
                                                     volume=True, speed=True, append=resuming))
    
    # 5. Production MD Loop
    print(f"Running production MD simulation (replicate seed {args.seed}) for {args.ns} ns...")
    steps_completed = 0
    if resuming and os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    last_line = lines[-1].split(',')
                    steps_completed = int(last_line[0])
                    print(f"Detected {steps_completed} steps already completed.")
        except Exception as e:
            print(f"Error reading log file: {e}. Starting steps count from 0.")
            steps_completed = 0
            
    t_start = time.time()
    
    while steps_completed < TOTAL_STEPS:
        chunk_steps = min(CHECKPOINT_INTERVAL, TOTAL_STEPS - steps_completed)
        t_chunk_start = time.time()
        
        simulation.step(chunk_steps)
        
        steps_completed += chunk_steps
        t_chunk_end = time.time()
        
        with open(checkpoint_file, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
            
        elapsed_chunk = t_chunk_end - t_chunk_start
        ns_completed = (steps_completed * TIMESTEP_FS * 1e-6)
        speed_ns_per_day = (chunk_steps * TIMESTEP_FS * 1e-15) / (elapsed_chunk / 86400.0) / 1e-9
        
        remaining_ns = args.ns - ns_completed
        if speed_ns_per_day > 0:
            eta_seconds = (remaining_ns / speed_ns_per_day) * 86400.0
            eta_str = str(datetime.timedelta(seconds=int(eta_seconds)))
        else:
            eta_str = "Calculating..."
            
        log_msg = (
            f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Progress: {ns_completed:.2f}/{args.ns:.2f} ns ({steps_completed:,}/{TOTAL_STEPS:,} steps) | "
            f"Speed: {speed_ns_per_day:.2f} ns/day | ETA: {eta_str} | "
            f"{get_gpu_utilization()}\n"
        )
        print(log_msg.strip())
        sys.stdout.flush()
        
        with open(args.progress_file, "a") as prog_f:
            prog_f.write(log_msg)
            prog_f.flush()
            
    with open(args.progress_file, "a") as prog_f:
        prog_f.write(f"----------------------------------------------------------\n")
        prog_f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] MD Simulation completed successfully.\n")
        prog_f.write(f"Final trajectory written to: {trajectory_file}\n")
        prog_f.write(f"==========================================================\n")
        prog_f.flush()
        
    print("Simulation completed successfully.")

if __name__ == "__main__":
    main()
