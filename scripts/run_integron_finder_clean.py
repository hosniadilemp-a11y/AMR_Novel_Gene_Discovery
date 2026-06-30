import os
import sys
import glob
import shutil
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor

def split_contigs(fasta_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    split_files = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Replace special characters in header to avoid path issues
        safe_id = "".join([c if c.isalnum() or c in "._-" else "_" for c in record.id])
        out_path = os.path.join(out_dir, f"{safe_id}.fasta")
        SeqIO.write(record, out_path, "fasta")
        split_files.append(out_path)
    return split_files

def run_integron_finder_on_file(fasta_path, out_dir):
    basename = os.path.basename(fasta_path)
    prefix = os.path.splitext(basename)[0]
    run_out_dir = os.path.join(out_dir, prefix)
    os.makedirs(run_out_dir, exist_ok=True)
    
    cmd = [
        "integron_finder",
        "--outdir", run_out_dir,
        "--cpu", "1",
        "--local-max",
        "--gbk",
        fasta_path
    ]
    
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            print(f"Successfully processed {prefix}")
            return prefix, run_out_dir
        else:
            print(f"Error processing {prefix}: {res.stderr}")
            return prefix, None
    except Exception as e:
        print(f"Exception processing {prefix}: {e}")
        return prefix, None

def merge_results(run_dirs, out_dir):
    summary_dfs = []
    integron_dfs = []
    
    for prefix, run_dir in run_dirs:
        if not run_dir:
            continue
        
        # IntegronFinder creates files matching: results/*.summary and results/*.integrons
        # inside the output directory structure
        summary_files = glob.glob(os.path.join(run_dir, "results", "*.summary"))
        integron_files = glob.glob(os.path.join(run_dir, "results", "*.integrons"))
        
        for sf in summary_files:
            try:
                df = pd.read_csv(sf, sep="\t", comment="#", header=None)
                # Assign column headers typical for IntegronFinder summary
                df.columns = ["id_replicon", "element", "start", "end", "strand", "type", "annotation", "attc_id"][:len(df.columns)]
                df["contig"] = prefix
                summary_dfs.append(df)
            except Exception:
                pass
                
        for inf in integron_files:
            try:
                df = pd.read_csv(inf, sep="\t", comment="#")
                df["contig"] = prefix
                integron_dfs.append(df)
            except Exception:
                pass

    merged_summary_path = os.path.join(out_dir, "merged_clean_summary.csv")
    merged_integrons_path = os.path.join(out_dir, "merged_clean_integrons.csv")
    
    if summary_dfs:
        merged_summary = pd.concat(summary_dfs, ignore_index=True)
        merged_summary.to_csv(merged_summary_path, index=False)
        print(f"Saved merged summary to {merged_summary_path}")
    else:
        # Create empty placeholder file
        pd.DataFrame().to_csv(merged_summary_path, index=False)
        
    if integron_dfs:
        merged_integrons = pd.concat(integron_dfs, ignore_index=True)
        merged_integrons.to_csv(merged_integrons_path, index=False)
        print(f"Saved merged integrons to {merged_integrons_path}")
    else:
        # Create empty placeholder file
        pd.DataFrame().to_csv(merged_integrons_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Clean Parallel IntegronFinder Runner")
    parser.add_argument("--contigs", required=True, help="Input multi-fasta assembly")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel files to process")
    args = parser.parse_args()

    temp_split_dir = os.path.join(args.output, "temp_split_contigs")
    temp_run_dir = os.path.join(args.output, "temp_runs")
    
    print("Splitting multi-fasta contigs into individual files...")
    split_files = split_contigs(args.contigs, temp_split_dir)
    print(f"Split assembly into {len(split_files)} individual contig files.")
    
    print(f"Running IntegronFinder with {args.threads} threads...")
    run_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(run_integron_finder_on_file, f, temp_run_dir) 
            for f in split_files
        ]
        for fut in futures:
            run_results.append(fut.result())
            
    print("Merging individual integron prediction files...")
    merge_results(run_results, args.output)
    
    # Cleanup temporary directories
    print("Cleaning up temporary directories...")
    if os.path.exists(temp_split_dir):
        shutil.rmtree(temp_split_dir)
    if os.path.exists(temp_run_dir):
        shutil.rmtree(temp_run_dir)
        
    print("Done!")

if __name__ == "__main__":
    main()
