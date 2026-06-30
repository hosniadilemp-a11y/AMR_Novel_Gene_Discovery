import os
import glob
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor

def annotate_genome(fasta_path, out_dir, cpus):
    filename = os.path.basename(fasta_path)
    prefix = os.path.splitext(filename)[0]
    out_prefix_dir = os.path.join(out_dir, prefix)
    
    # Check if Prokka GFF already exists
    gff_path = os.path.join(out_prefix_dir, f"{prefix}.gff")
    if os.path.exists(gff_path) and os.path.getsize(gff_path) > 1000:
        print(f"Annotation already exists for {prefix}. Skipping.")
        return gff_path
        
    print(f"Annotating {prefix}...")
    cmd = [
        "prokka",
        "--outdir", out_prefix_dir,
        "--prefix", prefix,
        "--cpus", str(cpus),
        "--force",
        fasta_path
    ]
    
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            print(f"Successfully annotated {prefix}.")
            return gff_path
        else:
            print(f"Error annotating {prefix}: {res.stderr}")
            return None
    except Exception as e:
        print(f"Exception for {prefix}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Parallel Prokka Reference Annotator")
    parser.add_argument("--indir", default="results/step3_annotation/reference_genomes", help="Input FASTA directory")
    parser.add_argument("--outdir", default="results/step3_annotation/reference_prokka", help="Output directory")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel genomes to annotate")
    parser.add_argument("--cpus-per-worker", type=int, default=4, help="CPUs per Prokka process")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    fasta_files = glob.glob(os.path.join(args.indir, "*.fna"))
    if not fasta_files:
        fasta_files = glob.glob(os.path.join(args.indir, "*.fasta"))
        
    print(f"Found {len(fasta_files)} reference genomes to annotate.")
    if not fasta_files:
        return

    print(f"Running Prokka annotation with {args.workers} parallel workers, each using {args.cpus_per_worker} CPUs...")
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(
            lambda f: annotate_genome(f, args.outdir, args.cpus_per_worker), 
            fasta_files
        ))
        
    success_count = sum(1 for r in results if r is not None)
    print(f"Completed! Successfully annotated {success_count}/{len(fasta_files)} reference genomes.")

if __name__ == "__main__":
    main()
