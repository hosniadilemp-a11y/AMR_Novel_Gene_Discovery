#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd

try:
    import pytraj as pt
except ImportError:
    print("[ERROR] pytraj is not installed. Please run this script in the openmm_env environment.")
    sys.exit(1)

# Paths
TRAJ_FILE = "md_apo_trajectory.dcd"
TOP_FILE = "solvated_system.pdb"

if not os.path.exists(TRAJ_FILE) or not os.path.exists(TOP_FILE):
    print(f"[ERROR] Trajectory ({TRAJ_FILE}) or Topology ({TOP_FILE}) file not found.")
    sys.exit(0)

print("Loading trajectory...")
# Load trajectory iteratively to save RAM
traj = pt.iterload(TRAJ_FILE, TOP_FILE)
print(f"Loaded trajectory with {traj.n_frames} frames, {traj.n_atoms} atoms per frame.")

# Select protein mask (residues 1 to 351)
# Note: In PyTraj/Amber, protein is usually mask ':1-351'
protein_mask = ":1-351"
print(f"Selecting protein residues ({protein_mask})...")

# Align protein backbone to the first frame to remove global translation/rotation
print("Aligning trajectory to first frame backbone...")
aligned_traj = pt.align(traj, mask=f"{protein_mask}@CA,C,N,O", ref=0)

# 1. Calculate Backbone RMSD
print("Calculating backbone RMSD...")
# pt.rmsd calculates rmsd against the reference frame (first frame by default)
rmsd = pt.rmsd(aligned_traj, mask=f"{protein_mask}@CA,C,N,O", ref=0)

# 2. Calculate C-alpha RMSF
print("Calculating C-alpha RMSF...")
# pt.rmsf calculates fluctuation per atom. We select C-alpha of residues 1-351
rmsf = pt.rmsf(aligned_traj, mask=f"{protein_mask}@CA")

# 3. Calculate Radius of Gyration (Rg)
print("Calculating Radius of Gyration (Rg)...")
rg = pt.radgyr(aligned_traj, mask=protein_mask)

# Save RMSD and Rg to CSV
print("Saving timeseries results to CSV...")
df_timeseries = pd.DataFrame({
    'Frame': np.arange(len(rmsd)),
    'Time_ns': np.arange(len(rmsd)) * 0.01,  # 10 ps interval = 0.01 ns
    'RMSD_Angstrom': rmsd,
    'Rg_Angstrom': rg
})
df_timeseries.to_csv("md_analysis_timeseries.csv", index=False)

# Save RMSF to CSV
df_rmsf = pd.DataFrame({
    'Residue': np.arange(1, len(rmsf) + 1),
    'Atom_Index': rmsf[:, 0].astype(int),
    'RMSF_Angstrom': rmsf[:, 1]
})
df_rmsf.to_csv("md_analysis_rmsf.csv", index=False)

# Print Summary Statistics
print("\n" + "="*40)
print("       MD SIMULATION SUMMARY STATISTICS")
print("="*40)
print(f"Total Simulation Time : {df_timeseries['Time_ns'].iloc[-1]:.2f} ns ({traj.n_frames} frames)")
print(f"Backbone RMSD (overall): {rmsd.mean():.4f} Å ± {rmsd.std():.4f} Å (Min: {rmsd.min():.2f} Å, Max: {rmsd.max():.2f} Å)")
print(f"Backbone RMSD (last 50ns): {rmsd[len(rmsd)//2:].mean():.4f} Å ± {rmsd[len(rmsd)//2:].std():.4f} Å")
print(f"Radius of Gyration (Rg): {rg.mean():.4f} Å ± {rg.std():.4f} Å")
print(f"C-alpha RMSF (Average) : {rmsf[:, 1].mean():.4f} Å (Max: {rmsf[:, 1].max():.2f} Å at Residue {rmsf[:, 1].argmax()+1})")
print("="*40 + "\n")

print("Analysis completed successfully. Output files: md_analysis_timeseries.csv, md_analysis_rmsf.csv")
