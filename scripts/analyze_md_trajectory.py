#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import pandas as pd

# Force pure OpenMM + NumPy implementation to prevent pytraj C++ crashes
PYTRAJ_AVAILABLE = False

# OpenMM import (required for fallback, and useful anyway)
try:
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    print("[ERROR] OpenMM is not installed. Please run this script in an environment with OpenMM installed.", file=sys.stderr)
    sys.exit(1)

def kabsch_align(P, Q):
    """
    Align coordinates P (N x 3) to reference coordinates Q (N x 3)
    using the Kabsch algorithm. Returns rotation matrix R and translation vector T
    such that P_aligned = P @ R + T.
    """
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    
    # Covariance matrix
    C = np.dot(P_centered.T, Q_centered)
    
    # SVD
    V, S, W_t = np.linalg.svd(C)
    
    # Rotation matrix
    d = np.linalg.det(V) * np.linalg.det(W_t.T)
    if d < 0.0:
        V[:, -1] = -V[:, -1]
    
    R = np.dot(V, W_t)
    T = centroid_Q - np.dot(centroid_P, R)
    return R, T

def analyze_with_pytraj(args):
    print("Analyzing trajectory using PyTraj...")
    traj = pt.iterload(args.traj, args.top)
    print(f"Loaded trajectory with {traj.n_frames} frames, {traj.n_atoms} atoms per frame.")

    # Determine residue count
    if args.residues:
        n_residues = args.residues
        print(f"Using user-specified residue count: {n_residues}")
    else:
        ca_indices = traj.top.select("@CA")
        if len(ca_indices) == 0:
            print("[ERROR] No CA atoms found in topology. Cannot auto-detect residues.", file=sys.stderr)
            sys.exit(1)
        n_residues = len(ca_indices)
        print(f"Auto-detected {n_residues} residues based on CA atom count.")

    protein_mask = f":1-{n_residues}"
    print(f"Protein selection mask: {protein_mask}")

    # Align protein backbone to the first frame
    print("Aligning trajectory to first frame backbone (CA, C, N, O)...")
    aligned_traj = pt.align(traj, mask=f"{protein_mask}@CA,C,N,O", ref=0)

    # Calculate metrics
    print("Calculating backbone RMSD...")
    rmsd = pt.rmsd(aligned_traj, mask=f"{protein_mask}@CA,C,N,O", ref=0)

    print("Calculating C-alpha RMSF...")
    rmsf = pt.rmsf(aligned_traj, mask=f"{protein_mask}@CA")

    print("Calculating Radius of Gyration (Rg)...")
    rg = pt.radgyr(aligned_traj, mask=protein_mask)

    return rmsd, rmsf[:, 0].astype(int), rmsf[:, 1], rg, n_residues, traj.n_frames

def read_dcd_coordinates_generator(traj_path, n_atoms):
    import struct
    with open(traj_path, "rb") as f:
        # Skip header block (4 bytes size, 84 bytes variables, 4 bytes size)
        f.seek(92, 0)
        
        # Skip title block
        title_size = struct.unpack("i", f.read(4))[0]
        f.seek(title_size + 4, 1) # Skip title + 4 bytes end marker
        
        # Skip atoms block
        f.seek(12, 1) # Skip size (4), n_atoms (4), size (4)
        
        # Frame loop
        while True:
            try:
                # Read box size marker
                box_size_data = f.read(4)
                if not box_size_data or len(box_size_data) < 4:
                    break
                box_size = struct.unpack("i", box_size_data)[0]
                if box_size == 48:
                    f.seek(52, 1) # Skip box (48) + size marker (4)
                else:
                    f.seek(-4, 1) # Go back box size marker
                    
                # Read X coords
                x_size_data = f.read(4)
                if len(x_size_data) < 4:
                    break
                x_size = struct.unpack("i", x_size_data)[0]
                x_coords = np.fromfile(f, dtype=np.float32, count=n_atoms)
                if len(x_coords) < n_atoms:
                    break
                f.seek(4, 1) # Skip end marker
                
                # Read Y coords
                y_size_data = f.read(4)
                if len(y_size_data) < 4:
                    break
                y_size = struct.unpack("i", y_size_data)[0]
                y_coords = np.fromfile(f, dtype=np.float32, count=n_atoms)
                if len(y_coords) < n_atoms:
                    break
                f.seek(4, 1) # Skip end marker
                
                # Read Z coords
                z_size_data = f.read(4)
                if len(z_size_data) < 4:
                    break
                z_size = struct.unpack("i", z_size_data)[0]
                z_coords = np.fromfile(f, dtype=np.float32, count=n_atoms)
                if len(z_coords) < n_atoms:
                    break
                f.seek(4, 1) # Skip end marker
                
                coords = np.column_stack((x_coords, y_coords, z_coords))
                yield coords
            except (struct.error, OSError, ValueError):
                break

def analyze_with_openmm(args):
    print("Analyzing trajectory using pure OpenMM + NumPy...")
    
    # Load topology
    pdb = app.PDBFile(args.top)
    topology = pdb.topology
    n_atoms_total = topology.getNumAtoms()
    print(f"Total atoms in topology: {n_atoms_total}")
    
    # Identify protein residues and atom indices
    amino_acids = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    backbone_atoms = []
    ca_atoms = []
    protein_atoms = []
    
    # Mapping for validation
    for atom in topology.atoms():
        res_name = atom.residue.name
        if res_name in amino_acids:
            protein_atoms.append(atom.index)
            if atom.name in {'CA', 'C', 'N', 'O'}:
                backbone_atoms.append(atom.index)
            if atom.name == 'CA':
                ca_atoms.append(atom.index)
                
    n_residues = len(ca_atoms)
    print(f"Detected {n_residues} protein residues in topology.")
    
    if args.residues and args.residues != n_residues:
        print(f"[WARNING] Overriding detected residue count {n_residues} with specified {args.residues}.")
        n_residues = args.residues
        # Slice lists to match user-specified count
        ca_atoms = ca_atoms[:n_residues]
        # Re-filter backbone and protein atoms to only include residues up to n_residues
        backbone_atoms = [
            atom.index for atom in topology.atoms() 
            if atom.residue.name in amino_acids and atom.residue.index < n_residues and atom.name in {'CA', 'C', 'N', 'O'}
        ]
        protein_atoms = [
            atom.index for atom in topology.atoms() 
            if atom.residue.name in amino_acids and atom.residue.index < n_residues
        ]

    # Lists for metrics
    rmsd = []
    rg = []
    aligned_ca_coords = []
    
    # Initialize generator
    gen = read_dcd_coordinates_generator(args.traj, n_atoms_total)
    
    # First frame coordinates for reference
    try:
        ref_coords = next(gen)
    except StopIteration:
        print("[ERROR] Trajectory has no frames.", file=sys.stderr)
        sys.exit(1)
        
    ref_backbone = ref_coords[backbone_atoms]
    
    # Calculate metrics for first frame
    ref_aligned_frame = ref_coords
    ref_rmsd = 0.0
    rmsd.append(ref_rmsd)
    
    # Radius of gyration (protein atoms only)
    prot_coords = ref_aligned_frame[protein_atoms]
    centroid = np.mean(prot_coords, axis=0)
    prot_centered = prot_coords - centroid
    ref_rg = np.sqrt(np.mean(np.sum(prot_centered ** 2, axis=1)))
    rg.append(ref_rg)
    aligned_ca_coords.append(ref_aligned_frame[ca_atoms])
    
    print("Aligning frames and calculating RMSD & Rg iteratively...")
    frame_idx = 1
    for frame_coords in gen:
        frame_backbone = frame_coords[backbone_atoms]
        
        # Calculate Kabsch alignment based on backbone
        R, T = kabsch_align(frame_backbone, ref_backbone)
        
        # Align protein coordinates (only rotate protein atoms to save cpu time)
        aligned_protein = np.dot(frame_coords[protein_atoms], R) + T
        aligned_backbone = np.dot(frame_backbone, R) + T
        
        # Backbone RMSD
        diff = aligned_backbone - ref_backbone
        frame_rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
        rmsd.append(frame_rmsd)
        
        # Radius of gyration (protein atoms only)
        centroid = np.mean(aligned_protein, axis=0)
        prot_centered = aligned_protein - centroid
        frame_rg = np.sqrt(np.mean(np.sum(prot_centered ** 2, axis=1)))
        rg.append(frame_rg)
        
        # Store aligned CA for RMSF
        aligned_ca = np.dot(frame_coords[ca_atoms], R) + T
        aligned_ca_coords.append(aligned_ca)
        
        frame_idx += 1
        if frame_idx % 1000 == 0:
            print(f"Processed {frame_idx} frames...")
            
    n_frames = frame_idx
    print(f"Completed processing of {n_frames} frames.")
    
    rmsd = np.array(rmsd)
    rg = np.array(rg)
    aligned_ca_coords = np.array(aligned_ca_coords) # Shape: (n_frames, n_residues, 3)
    
    # Calculate C-alpha RMSF
    print("Calculating C-alpha RMSF...")
    mean_ca = np.mean(aligned_ca_coords, axis=0) # Shape: (n_residues, 3)
    fluctuations = aligned_ca_coords - mean_ca # Shape: (n_frames, n_residues, 3)
    rmsf_vals = np.sqrt(np.mean(np.sum(fluctuations ** 2, axis=2), axis=0)) # Shape: (n_residues,)
    
    # Pytraj atom indices are 1-based, we map them here for consistency
    atom_indices = np.array(ca_atoms) + 1
    
    return rmsd, atom_indices, rmsf_vals, rg, n_residues, n_frames

def main():
    parser = argparse.ArgumentParser(description="Analyze OpenMM MD Trajectories (RMSD, RMSF, Rg)")
    parser.add_argument("--traj", required=True, help="Path to trajectory DCD file")
    parser.add_argument("--top", required=True, help="Path to topology solvated PDB file")
    parser.add_argument("--outdir", required=True, help="Directory to save output CSVs")
    parser.add_argument("--residues", type=int, default=None, help="Number of protein residues (optional)")
    args = parser.parse_args()

    if not os.path.exists(args.traj):
        print(f"[ERROR] Trajectory file '{args.traj}' not found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.top):
        print(f"[ERROR] Topology file '{args.top}' not found.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # Perform analysis
    if PYTRAJ_AVAILABLE:
        try:
            rmsd, atom_indices, rmsf, rg, n_residues, n_frames = analyze_with_pytraj(args)
        except Exception as e:
            print(f"[WARNING] PyTraj analysis failed: {e}. Falling back to OpenMM.")
            rmsd, atom_indices, rmsf, rg, n_residues, n_frames = analyze_with_openmm(args)
    else:
        rmsd, atom_indices, rmsf, rg, n_residues, n_frames = analyze_with_openmm(args)

    # Output paths
    timeseries_csv = os.path.join(args.outdir, "md_analysis_timeseries.csv")
    rmsf_csv = os.path.join(args.outdir, "md_analysis_rmsf.csv")

    # Save RMSD and Rg to CSV
    print(f"Saving timeseries results to {timeseries_csv}...")
    df_timeseries = pd.DataFrame({
        'Frame': np.arange(len(rmsd)),
        'Time_ns': np.arange(len(rmsd)) * 0.01,  # 10 ps interval = 0.01 ns
        'RMSD_Angstrom': rmsd,
        'Rg_Angstrom': rg
    })
    df_timeseries.to_csv(timeseries_csv, index=False)

    # Save RMSF to CSV
    print(f"Saving RMSF results to {rmsf_csv}...")
    df_rmsf = pd.DataFrame({
        'Residue': np.arange(1, len(rmsf) + 1),
        'Atom_Index': atom_indices,
        'RMSF_Angstrom': rmsf
    })
    df_rmsf.to_csv(rmsf_csv, index=False)

    # Print Summary Statistics
    print("\n" + "="*40)
    print("       MD SIMULATION SUMMARY STATISTICS")
    print("="*40)
    print(f"Total Simulation Time : {df_timeseries['Time_ns'].iloc[-1]:.2f} ns ({n_frames} frames)")
    print(f"Backbone RMSD (overall): {rmsd.mean():.4f} Å ± {rmsd.std():.4f} Å (Min: {rmsd.min():.2f} Å, Max: {rmsd.max():.2f} Å)")
    print(f"Backbone RMSD (last 50ns): {rmsd[len(rmsd)//2:].mean():.4f} Å ± {rmsd[len(rmsd)//2:].std():.4f} Å")
    print(f"Radius of Gyration (Rg): {rg.mean():.4f} Å ± {rg.std():.4f} Å")
    print(f"C-alpha RMSF (Average) : {rmsf.mean():.4f} Å (Max: {rmsf.max():.2f} Å at Residue {rmsf.argmax()+1})")
    print("="*40 + "\n")

    print("Analysis completed successfully.")

if __name__ == "__main__":
    main()
