#!/usr/bin/env python3
# =============================================================================
# 10_docking.py — Molecular Docking & Specificity Validation
# =============================================================================
#
# USAGE:
#   python3 scripts/10_docking.py [OPTIONS]
#
# =============================================================================
import os
import sys
import subprocess
import re
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
import pandas as pd

LIGANDS = {
    "Kanamycin": {
        "smiles": "C1C(C(C(C(C1N)OC2C(C(C(C(O2)CO)O)O)N)O)OC3C(C(C(C(O3)CN)O)O)O)N",
        "type": "Positive Control"
    },
    "Gentamicin": {
        "smiles": "CC(C1CCC(C(O1)OC2C(CC(C(C2O)OC3C(C(CO3)(O)C)NC)N)N)N)NC",
        "type": "Positive Control"
    },
    "Amikacin": {
        "smiles": "C1C(C(C(C(C1N)OC2C(C(C(C(O2)CO)O)O)N)O)OC3C(C(C(C(O3)CN)O)O)O)NC(=O)C(CCN)O",
        "type": "Positive Control"
    },
    "Penicillin_G": {
        "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
        "type": "Negative Decoy"
    },
    "Tetracycline": {
        "smiles": "CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O",
        "type": "Negative Decoy"
    },
    "D_Glucose": {
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "type": "Negative Decoy"
    }
}

def prepare_receptor(pdb_in, pdbqt_out):
    """Convert heavy-atom PDB structure to PDBQT format."""
    aromatic_residues = {"PHE", "TYR", "TRP", "HIS"}
    aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2", "NE2", "ND1"}
    
    lines = []
    with open(pdb_in, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element = line[76:78].strip()
                if not element:
                    element = atom_name[0]
                    
                ad_type = element
                if element == "C":
                    if res_name in aromatic_residues and atom_name in aromatic_atoms:
                        ad_type = "A"
                
                occ = line[54:60] if len(line) >= 60 else "  1.00"
                temp = line[60:66] if len(line) >= 66 else "  0.00"
                
                pdbqt_line = f"{line[:54]}{occ}{temp}      0.000 {ad_type:<2}\n"
                lines.append(pdbqt_line)
            elif line.startswith("TER") or line.startswith("END"):
                lines.append(line)
                
    with open(pdbqt_out, "w") as f:
        f.writelines(lines)

def prepare_ligand(name, smiles, pdbqt_out):
    """Generate 3D conformations and convert SMILES to PDBQT using RDKit and Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    AllChem.EmbedMolecule(mol, ps)
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparer = MoleculePreparation()
    preparer.prepare(mol)
    preparer.write_pdbqt_file(pdbqt_out)

def run_vina(vina_bin, receptor_path, ligand_path, out_path, log_path, cx, cy, cz, sx, sy, sz):
    cmd = [
        vina_bin,
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--out", out_path,
        "--exhaustiveness", "8"
    ]
    with open(log_path, "w") as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.PIPE, text=True, check=True)

def parse_vina_affinity(log_path):
    affinity = 0.0
    with open(log_path, "r") as f:
        for line in f:
            match = re.search(r"^\s*1\s+([\-\d\.]+)", line)
            if match:
                affinity = float(match.group(1))
                break
    return affinity

def main():
    parser = argparse.ArgumentParser(description="AutoDock Vina specificity docking script")
    parser.add_argument("--receptor", required=True, help="Input receptor PDB file")
    parser.add_argument("--vina", default="vina", help="Path to Vina binary")
    parser.add_argument("--grid-center", default="2.964 2.776 2.588", help="Grid center X Y Z")
    parser.add_argument("--grid-size", default="22 22 22", help="Grid box size X Y Z")
    parser.add_argument("--output", default="results/step7_docking", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    receptor_pdbqt = os.path.join(args.output, "receptor.pdbqt")

    print(f"Preparing receptor PDBQT: {receptor_pdbqt}")
    prepare_receptor(args.receptor, receptor_pdbqt)

    cx, cy, cz = map(float, args.grid_center.split())
    sx, sy, sz = map(float, args.grid_size.split())

    results = []
    for name, cfg in LIGANDS.items():
        print(f"Processing ligand: {name}")
        ligand_pdbqt = os.path.join(args.output, f"{name}.pdbqt")
        out_pdbqt = os.path.join(args.output, f"{name}_docked.pdbqt")
        log_file = os.path.join(args.output, f"{name}_vina.log")

        prepare_ligand(name, cfg["smiles"], ligand_pdbqt)
        
        try:
            run_vina(args.vina, receptor_pdbqt, ligand_pdbqt, out_pdbqt, log_file, cx, cy, cz, sx, sy, sz)
            affinity = parse_vina_affinity(log_file)
            print(f"  Affinity: {affinity} kcal/mol")
        except Exception as e:
            print(f"  Vina failed for {name}: {e}. Skipping Vina step for this ligand.")
            affinity = 0.0

        results.append({
            "Ligand": name,
            "Type": cfg["type"],
            "SMILES": cfg["smiles"],
            "Affinity": affinity
        })

    tsv_path = os.path.join(args.output, "docking_specificity_results.tsv")
    df = pd.DataFrame(results)
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"Results written to {tsv_path}")

if __name__ == "__main__":
    main()
