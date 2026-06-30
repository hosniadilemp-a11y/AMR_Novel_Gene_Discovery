#!/usr/bin/env python3
# =============================================================================
# 06b_foldseek_search.py — Foldseek Structural Homology Search
# =============================================================================
#
# USAGE:
#   python3 scripts/06b_foldseek_search.py [OPTIONS]
#
# =============================================================================
import os
import sys
import time
import json
import csv
import argparse
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

def get_session():
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=2,
        status_forcelist=[500, 502, 503, 504],
        raise_on_status=False
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.mount("http://", HTTPAdapter(max_retries=retries))
    return session

def submit_foldseek_job(session, pdb_path, databases):
    url = "https://search.foldseek.com/api/ticket"
    data = [('mode', '3diaa')]
    for db in databases:
        data.append(('database[]', db))
        
    try:
        with open(pdb_path, 'rb') as f:
            files = {'q': f}
            response = session.post(url, data=data, files=files, timeout=30)
            
        if response.status_code == 200:
            return response.json()
        else:
            print(f"  Failed submission for {os.path.basename(pdb_path)}: Status {response.status_code} - {response.text}")
            return None
    except Exception as e:
        print(f"  Error submitting {os.path.basename(pdb_path)}: {str(e)}")
        return None

def check_job_status(session, ticket_id):
    url = f"https://search.foldseek.com/api/ticket/{ticket_id}"
    try:
        response = session.get(url, timeout=15)
        if response.status_code == 200:
            return response.json()
        else:
            return {"status": "ERROR", "msg": f"Status code {response.status_code}"}
    except Exception as e:
        return {"status": "ERROR", "msg": str(e)}

def download_results(session, ticket_id):
    url = f"https://search.foldseek.com/api/result/{ticket_id}/0"
    try:
        response = session.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"  Failed retrieving results for ticket {ticket_id}: Status {response.status_code}")
            return None
    except Exception as e:
        print(f"  Error retrieving results for ticket {ticket_id}: {str(e)}")
        return None

def process_candidate(session, pdb_path, out_dir, databases):
    cand_id = os.path.basename(pdb_path).replace(".pdb", "")
    json_path = os.path.join(out_dir, f"{cand_id}.json")
    
    if os.path.exists(json_path):
        print(f"[{cand_id}] Results already exist locally. Skipping API query.")
        return True
        
    print(f"[{cand_id}] Submitting PDB to Foldseek API...")
    ticket = submit_foldseek_job(session, pdb_path, databases)
    if not ticket:
        return False
        
    ticket_id = ticket.get("id")
    print(f"[{cand_id}] Ticket created: {ticket_id}. Polling status...")
    
    start_time = time.time()
    timeout_limit = 600
    
    while (time.time() - start_time) < timeout_limit:
        time.sleep(5)
        status_info = check_job_status(session, ticket_id)
        status = status_info.get("status")
        
        if status == "COMPLETE":
            print(f"[{cand_id}] Job complete! Fetching results...")
            results = download_results(session, ticket_id)
            if results:
                with open(json_path, 'w') as out_f:
                    json.dump(results, out_f, indent=2)
                print(f"[{cand_id}] Results written to {os.path.basename(json_path)}")
                return True
            else:
                return False
        elif status in ["ERROR", "FAILED"]:
            print(f"[{cand_id}] Job failed: {status_info.get('msg', 'Unknown Error')}")
            return False
            
    print(f"[{cand_id}] Job timed out after {timeout_limit} seconds.")
    return False

def parse_results(out_dir, summary_path):
    print("\nParsing JSON results to compile summary table...")
    summaries = []
    
    json_files = sorted([f for f in os.listdir(out_dir) if f.endswith(".json") and f != "summary.json"])
    
    for j_file in json_files:
        cand_id = j_file.replace(".json", "")
        json_path = os.path.join(out_dir, j_file)
        
        with open(json_path) as f:
            data = json.load(f)
            
        all_hits = []
        results_list = data.get("results", [])
        
        for db_res in results_list:
            db_name = db_res.get("db", "unknown")
            alignments = db_res.get("alignments", [])
            for chain_alns in alignments:
                for hit in chain_alns:
                    target = hit.get("target", "")
                    parts = target.split(" ", 1)
                    target_id = parts[0]
                    target_desc = parts[1] if len(parts) > 1 else "Uncharacterized protein"
                    
                    e_val = hit.get("eval", 10.0)
                    prob = hit.get("prob", 0.0)
                    seq_id = hit.get("seqId", 0.0)
                    aln_len = hit.get("alnLength", 0)
                    tax_name = hit.get("taxName", "Unknown")
                    score = hit.get("score", 0)
                    
                    all_hits.append({
                        "db": db_name,
                        "target_id": target_id,
                        "target_desc": target_desc,
                        "evalue": e_val,
                        "prob": prob,
                        "identity": seq_id,
                        "aln_len": aln_len,
                        "tax_name": tax_name,
                        "score": score
                    })
                    
        all_hits.sort(key=lambda x: (-x["prob"], x["evalue"]))
        
        if all_hits:
            summaries.append({
                "cand_id": cand_id,
                "has_hits": True,
                "top_hit": all_hits[0],
                "hits_count": len(all_hits)
            })
        else:
            summaries.append({
                "cand_id": cand_id,
                "has_hits": False,
                "hits_count": 0
            })
            
    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Candidate_ID", "Top_Hit_ID", "Database", "Probability", "E-Value", "Identity_%", "Aln_Length", "Taxonomy", "Target_Description"])
        
        for s in summaries:
            if s["has_hits"]:
                t = s["top_hit"]
                writer.writerow([
                    s["cand_id"],
                    t["target_id"],
                    t["db"],
                    f"{t['prob']:.4f}",
                    f"{t['evalue']:.2e}",
                    f"{t['identity']*100:.2f}%" if t['identity'] <= 1.0 else f"{t['identity']:.2f}%",
                    t["aln_len"],
                    t["tax_name"],
                    t["target_desc"]
                ])
            else:
                writer.writerow([s["cand_id"], "No structural homology found", "NA", "0.00", "NA", "0.00%", 0, "NA", "NA"])
                
    print(f"Summary table written to {summary_path}")

def main():
    parser = argparse.ArgumentParser(description="Foldseek API Structural Homology Search")
    parser.add_argument("--pdb-dir", required=True, help="Directory containing ESMFold PDB structures")
    parser.add_argument("--databases", default="pdb100,afdb50,afdb-swissprot", help="Comma-separated databases to search")
    parser.add_argument("--output", default="results/step6_structures/foldseek_results", help="Output directory for JSONs")
    parser.add_argument("--summary", default="results/step6_structures/foldseek_summary.tsv", help="Output summary TSV file")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    os.makedirs(os.path.dirname(args.summary), exist_ok=True)

    pdb_files = sorted([f for f in os.listdir(args.pdb_dir) if f.endswith(".pdb")])
    print(f"Found {len(pdb_files)} candidate structure PDBs.")
    if not pdb_files:
        return

    databases = args.databases.split(",")
    session = get_session()
    success_count = 0
    
    for idx, pdb_name in enumerate(pdb_files):
        pdb_path = os.path.join(args.pdb_dir, pdb_name)
        print(f"\nProcessing [{idx+1}/{len(pdb_files)}]: {pdb_name}")
        
        success = process_candidate(session, pdb_path, args.output, databases)
        if success:
            success_count += 1
            
        time.sleep(2)
        
    print(f"\nFinished structural homology searches. Success: {success_count}/{len(pdb_files)}")
    parse_results(args.output, args.summary)

if __name__ == "__main__":
    main()
