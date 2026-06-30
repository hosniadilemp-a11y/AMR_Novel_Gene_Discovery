#!/usr/bin/env python3
# =============================================================================
# 09a_md_setup.py — Molecular Dynamics System Preparation
# =============================================================================
#
# USAGE:
#   python3 scripts/09a_md_setup.py [OPTIONS]
#
# =============================================================================
import os
import sys
import argparse
from pdbfixer import PDBFixer
import openmm.app as app
import openmm.unit as unit

def main():
    parser = argparse.ArgumentParser(description="Prepare and solvate protein structures for MD")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--output", default="results/step9_md", help="Output directory")
    parser.add_argument("--padding", type=float, default=1.2, help="Water box padding in nm")
    parser.add_argument("--strength", type=float, default=0.15, help="NaCl ionic strength in M")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    basename = os.path.basename(args.pdb)
    prefix = os.path.splitext(basename)[0]
    solvated_pdb = os.path.join(args.output, f"{prefix}_solvated.pdb")

    if os.path.exists(solvated_pdb):
        print(f"Solvated structure already exists: {solvated_pdb}. Skipping.")
        return

    print(f"Preparing structure {args.pdb} using PDBFixer...")
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # pH 7.0

    print("Solvating system with TIP3P water box and neutralizing ions...")
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = app.Modeller(fixer.topology, fixer.positions)
    
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=args.padding * unit.nanometers,
        ionicStrength=args.strength * unit.molar
    )

    print(f"Saving solvated topology to {solvated_pdb}...")
    with open(solvated_pdb, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print("System setup complete.")

if __name__ == "__main__":
    main()
