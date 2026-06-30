#!/usr/bin/env python3
# =============================================================================
# 06c_tmalign_validation.py — Local TM-align Structural Validation
# =============================================================================
#
# USAGE:
#   python3 scripts/06c_tmalign_validation.py [OPTIONS]
#
# =============================================================================
import os
import sys
import re
import csv
import argparse
import subprocess
import requests

def download_reference_pdb(target_id, target_dir):
    os.makedirs(target_dir, exist_ok=True)
    target_path = os.path.join(target_dir, f"{target_id}.pdb")
    
    if os.path.exists(target_path):
        return target_path
        
    print(f"Downloading reference structure {target_id} from AlphaFold DB...")
    url = f"https://alphafold.ebi.ac.uk/files/{target_id}.pdb"
    
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(target_path, 'wb') as f:
                f.write(response.content)
            print(f"  Saved to {target_path}")
            return target_path
        else:
            print(f"  Failed to download {target_id}: Status {response.status_code}")
            return None
    except Exception as e:
        print(f"  Error downloading {target_id}: {str(e)}")
        return None

def run_tmalign(query_path, target_path):
    # Try finding TMalign in PATH first, else fallback to standard location
    cmd = ["TMalign", query_path, target_path]
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        stdout = result.stdout
        
        aln_match = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d\.]+)", stdout)
        tm_match = re.search(r"TM-score=\s*([\d\.]+)\s*\(normalized by length of Structure_1", stdout)
        
        aln_len = int(aln_match.group(1)) if aln_match else 0
        rmsd = float(aln_match.group(2)) if aln_match else 0.0
        tm_score = float(tm_match.group(1)) if tm_match else 0.0
        
        return tm_score, rmsd, aln_len
    except subprocess.CalledProcessError as e:
        print(f"  TMalign execution failed: {e.stderr}")
        return 0.0, 0.0, 0
    except Exception as e:
        print(f"  Error running TMalign: {str(e)}")
        return 0.0, 0.0, 0

def main():
    parser = argparse.ArgumentParser(description="TM-align structural validation tool")
    parser.add_argument("--query-pdb", required=True, help="Directory containing ESMFold PDB structures")
    parser.add_argument("--target-pdb", default="results/step6_structures/target_structures", help="Directory to save reference PDBs")
    parser.add_argument("--foldseek", default="results/step6_structures/foldseek_summary.tsv", help="Path to Foldseek summary TSV")
    parser.add_argument("--output", default="results/step6_structures/tmalign_results.tsv", help="Output validated summary TSV")
    args = parser.parse_args()

    if not os.path.exists(args.foldseek):
        print(f"Error: Foldseek summary '{args.foldseek}' not found. Run Foldseek script first.")
        sys.exit(1)

    os.makedirs(args.target_pdb, exist_ok=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    rows = []
    with open(args.foldseek) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            rows.append(row)

    print(f"Loaded {len(rows)} candidates for TMalign validation.")
    updated_rows = []

    for idx, row in enumerate(rows):
        cand_id = row[0]
        top_hit_id = row[1]
        db = row[2]
        prob = row[3]
        evalue = row[4]
        identity = row[5]
        aln_len = row[6]
        taxonomy = row[7]
        target_desc = row[8]

        print(f"\n[{idx+1}/{len(rows)}] Candidate: {cand_id} vs {top_hit_id}")

        # Check if we have a valid AlphaFold model ID (typically format AF-XXXXXX-F1-model_vX)
        if top_hit_id == "No structural homology found" or not top_hit_id.startswith("AF-"):
            print("  Skipping TMalign (no valid AlphaFold DB hit).")
            updated_rows.append([cand_id, top_hit_id, db, prob, evalue, identity, aln_len, taxonomy, "0.0000", "0.00", target_desc])
            continue

        query_path = os.path.join(args.query_pdb, f"{cand_id}.pdb")
        target_path = download_reference_pdb(top_hit_id, args.target_pdb)

        if not target_path or not os.path.exists(query_path):
            print("  Required query/target PDB files not found. Skipping.")
            updated_rows.append([cand_id, top_hit_id, db, prob, evalue, identity, aln_len, taxonomy, "0.0000", "0.00", target_desc])
            continue

        tm_score, rmsd, aln_length = run_tmalign(query_path, target_path)
        print(f"  TMalign: TM-score = {tm_score:.4f}, RMSD = {rmsd:.2f} (Aligned: {aln_length})")

        updated_rows.append([
            cand_id, top_hit_id, db, prob, evalue, identity, aln_len, taxonomy,
            f"{tm_score:.4f}", f"{rmsd:.2f}", target_desc
        ])

    new_header = [
        "Candidate_ID", "Top_Hit_ID", "Database", "Probability", "E-Value", 
        "Identity_%", "Aln_Length", "Taxonomy", "TM-score", "RMSD", "Target_Description"
    ]
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(new_header)
        writer.writerows(updated_rows)

    print(f"\nTMalign validation complete. Results written to {args.output}")

if __name__ == "__main__":
    main()
