#!/usr/bin/env python3
# =============================================================================
# 05_novel_candidates.py — Novel Candidate Extraction and Homology Screening
# =============================================================================
#
# USAGE:
#   python3 scripts/05_novel_candidates.py [OPTIONS]
#
# =============================================================================
import os
import sys
import shutil
import urllib.request
import urllib.parse
import subprocess
import time
import json
import csv
import argparse
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

def download_file(url, dest):
    print(f"Downloading {url} to {dest}...")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with urllib.request.urlopen(url) as response, open(dest, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    print("Download complete.")

def submit_iprscan_job(email, seq_record):
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
    data = urllib.parse.urlencode({
        'email': email,
        'sequence': f">{seq_record.id}\n{str(seq_record.seq)}",
        'appl': 'PfamA'
    }).encode('utf-8')
    
    req = urllib.request.Request(url, data=data)
    try:
        with urllib.request.urlopen(req) as response:
            job_id = response.read().decode('utf-8').strip()
            return job_id
    except Exception as e:
        print(f"Error submitting sequence {seq_record.id}: {e}")
        return None

def check_iprscan_status(job_id):
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    try:
        with urllib.request.urlopen(url) as response:
            return response.read().decode('utf-8').strip()
    except Exception as e:
        return "ERROR"

def get_iprscan_results(job_id):
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
    try:
        with urllib.request.urlopen(url) as response:
            return json.loads(response.read().decode('utf-8'))
    except Exception as e:
        print(f"Error fetching results for job {job_id}: {e}")
        return None

def scan_pfam_online(seq_records, email="genomics_agent@ebi.ac.uk"):
    print(f"Scanning {len(seq_records)} sequences against Pfam via EBI API...")
    jobs = {}
    for rec in seq_records:
        job_id = submit_iprscan_job(email, rec)
        if job_id:
            jobs[rec.id] = job_id
            print(f"  Submitted {rec.id} -> Job ID: {job_id}")
        else:
            print(f"  Failed to submit {rec.id}")
            
    # Poll for results
    results = {}
    pending = list(jobs.keys())
    
    # Wait loop
    print("Waiting for jobs to finish (polling every 10s)...")
    start_time = time.time()
    while pending and (time.time() - start_time) < 600: # 10 min timeout
        time.sleep(10)
        still_pending = []
        for seq_id in pending:
            job_id = jobs[seq_id]
            status = check_iprscan_status(job_id)
            if status in ("FINISHED", "ERROR", "FAILURE"):
                if status == "FINISHED":
                    res_json = get_iprscan_results(job_id)
                    if res_json:
                        # Extract Pfam signatures
                        pfam_hits = []
                        for match in res_json.get("results", []):
                            sig = match.get("signature", {})
                            lib = sig.get("signatureLibraryRelease", {}).get("library", "")
                            if "pfam" in lib.lower():
                                pfam_hits.append({
                                    "ac": sig.get("accession"),
                                    "name": sig.get("name"),
                                    "desc": sig.get("description")
                                })
                        results[seq_id] = pfam_hits
                        print(f"  Job {seq_id} finished. Found {len(pfam_hits)} Pfam hits.")
                    else:
                        results[seq_id] = []
                else:
                    print(f"  Job {seq_id} failed with status {status}")
                    results[seq_id] = []
            else:
                still_pending.append(seq_id)
        pending = still_pending
        if pending:
            print(f"  {len(pending)} jobs still running...")
            
    # Fill in any timed-out jobs
    for seq_id in pending:
        print(f"  Job {seq_id} timed out.")
        results[seq_id] = []
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract and screen novel pangenome candidates")
    parser.add_argument("--prokka-dir", required=True, help="Directory containing Prokka outputs")
    parser.add_argument("--pangenome", required=True, help="Panaroo presence/absence CSV file")
    parser.add_argument("--min-length", type=int, default=200, help="Minimum amino acid length")
    parser.add_argument("--output", default="results/step5_candidates", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Resolve input Prokka files
    faa_path = glob.glob(os.path.join(args.prokka_dir, "*.faa"))
    ffn_path = glob.glob(os.path.join(args.prokka_dir, "*.ffn"))
    gff_path = glob.glob(os.path.join(args.prokka_dir, "*.gff"))

    if not faa_path or not os.path.exists(faa_path[0]):
        print(f"Error: Prokka protein FASTA (.faa) not found in {args.prokka_dir}")
        sys.exit(1)
    faa_path = faa_path[0]
    ffn_path = ffn_path[0] if ffn_path else None
    gff_path = gff_path[0] if gff_path else None

    # Parse Pangenome CSV to extract private singletons
    print("Loading pangenome presence/absence matrix...")
    if not os.path.exists(args.pangenome):
        print(f"Error: Pangenome CSV '{args.pangenome}' not found.")
        sys.exit(1)

    singletons = set()
    with open(args.pangenome) as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # Check column names. Roary and Panaroo GFF columns differ.
        # Typically the first 14 columns are metadata, rest are genomes.
        # Find the column for the focus strain QA5221
        focal_cols = [i for i, h in enumerate(header) if "QA5221" in h or "Ecoli_isolate" in h]
        if not focal_cols:
            print("ERROR: Focal strain column not found in pangenome headers.")
            sys.exit(1)
        focal_idx = focal_cols[0]
        
        # Identify reference genome columns (all genome columns except focal_idx)
        # Panaroo GFF genome columns start at index 14
        start_idx = 14 if len(header) > 14 else 3
        genome_indices = list(range(start_idx, len(header)))
        ref_indices = [idx for idx in genome_indices if idx != focal_idx]
        
        for row in reader:
            if len(row) < len(header):
                continue
            focal_val = row[focal_idx].strip()
            # If the gene is present in the focal strain and absent in all reference strains
            if focal_val and all(not row[idx].strip() for idx in ref_indices):
                # Handle cases where multiple genes are semicolon-separated
                for gene in focal_val.split(";"):
                    if gene.strip():
                        singletons.add(gene.strip())

    print(f"Found {len(singletons)} unique singleton genes in focal genome.")

    # Filter Prokka annotations
    records = list(SeqIO.parse(faa_path, "fasta"))
    hypo_min = [r for r in records if "hypothetical protein" in r.description.lower() and len(r.seq) >= args.min_length]
    print(f"Found {len(hypo_min)} hypothetical proteins >= {args.min_length} aa.")

    # Intersect
    candidates = [r for r in hypo_min if r.id in singletons]
    print(f"Intersection (Hypothetical + Singleton): {len(candidates)}")

    if not candidates:
        print("No candidates found after intersection. Exiting.")
        sys.exit(0)

    candidates_faa = os.path.join(args.output, "candidates.faa")
    SeqIO.write(candidates, candidates_faa, "fasta")

    # Local Swiss-Prot setting
    db_dir = "/tmp/swiss_prot"
    sprot_gz = os.path.join(db_dir, "uniprot_sprot.fasta.gz")
    sprot_fasta = os.path.join(db_dir, "uniprot_sprot.fasta")
    
    if not os.path.exists(sprot_fasta):
        if not os.path.exists(sprot_gz):
            os.makedirs(db_dir, exist_ok=True)
            sprot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
            try:
                download_file(sprot_url, sprot_gz)
            except Exception as e:
                print(f"Error downloading Swiss-Prot: {e}")
                sys.exit(1)
        print(f"Decompressing {sprot_gz}...")
        subprocess.run(["gzip", "-d", sprot_gz], check=True)
        
    # Build database
    if not os.path.exists(sprot_fasta + ".pdb"):
        print("Building local BLAST database...")
        subprocess.run(["makeblastdb", "-in", sprot_fasta, "-dbtype", "prot"], check=True)

    # Run BLASTp
    print("Running local BLASTp homology filter...")
    blast_out = os.path.join(args.output, "blast_hits.tsv")
    cmd = [
        "blastp",
        "-query", candidates_faa,
        "-db", sprot_fasta,
        "-out", blast_out,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-evalue", "1e-5",
        "-num_threads", "8"
    ]
    subprocess.run(cmd, check=True)

    # Parse BLASTp hits
    blast_hits = set()
    if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
        with open(blast_out) as f:
            for line in f:
                qseqid = line.strip().split("\t")[0]
                blast_hits.add(qseqid)

    print(f"BLASTp filtered out {len(blast_hits)} candidates with database matches.")
    blast_filtered = [c for c in candidates if c.id not in blast_hits]
    print(f"Candidates remaining after BLASTp: {len(blast_filtered)}")

    if not blast_filtered:
        print("No candidates left after BLASTp screening.")
        sys.exit(0)

    # Scan Pfam via EBI REST API
    pfam_hits = scan_pfam_online(blast_filtered)
    
    final_candidates = []
    skipped_pfam = 0
    for c in blast_filtered:
        hits = pfam_hits.get(c.id, [])
        if hits:
            print(f"  Skipping {c.id} due to Pfam domain matches: {', '.join([h['ac'] for h in hits])}")
            skipped_pfam += 1
        else:
            final_candidates.append(c)

    print(f"Filtered out {skipped_pfam} candidates with known Pfam domains.")
    print(f"FINAL NOVEL CANDIDATES REMAINING: {len(final_candidates)}")

    # Write final outputs
    final_faa = os.path.join(args.output, "prioritized_candidates.faa")
    final_tsv = os.path.join(args.output, "prioritized_candidates.tsv")

    if final_candidates:
        SeqIO.write(final_candidates, final_faa, "fasta")
        
        # Compute GC content and locate contig coordinates
        nuc_seqs = {}
        if ffn_path and os.path.exists(ffn_path):
            nuc_seqs = SeqIO.to_dict(SeqIO.parse(ffn_path, "fasta"))
            
        coords = {}
        if gff_path and os.path.exists(gff_path):
            with open(gff_path) as f:
                for line in f:
                    if line.startswith("##FASTA"):
                        break
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 9 or parts[2] != "CDS":
                        continue
                    locus_tag = ""
                    for attr in parts[8].split(";"):
                        if attr.startswith("ID="):
                            locus_tag = attr.split("=")[1]
                    if locus_tag in [c.id for c in final_candidates]:
                        coords[locus_tag] = {
                            "contig": parts[0],
                            "start": int(parts[3]),
                            "end": int(parts[4]),
                            "strand": parts[6]
                        }

        # Write metadata TSV
        with open(final_tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Candidate_ID", "Length_AA", "Contig", "Start", "End", "Strand", "GC_Content", "Description"])
            for c in final_candidates:
                info = coords.get(c.id, {"contig": "NA", "start": "NA", "end": "NA", "strand": "NA"})
                nuc_rec = nuc_seqs.get(c.id)
                gc = 0
                if nuc_rec:
                    seq_str = str(nuc_rec.seq).upper()
                    gc = (seq_str.count("G") + seq_str.count("C")) / len(seq_str) if len(seq_str) > 0 else 0
                writer.writerow([
                    c.id,
                    len(c.seq),
                    info["contig"],
                    info["start"],
                    info["end"],
                    info["strand"],
                    f"{gc*100:.2f}%",
                    c.description
                ])
        print(f"Prioritized candidates metadata saved to {final_tsv}")
    else:
        print("No candidates survived pangenomic and homology filtering.")

if __name__ == "__main__":
    main()
