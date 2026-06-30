#!/usr/bin/env python3
# =============================================================================
# 06a_esmfold_prediction.py — High-Throughput ESMFold Prediction (GPU)
# =============================================================================
#
# INSTRUCTIONS FOR RUNNING ON KAGGLE GPU:
#   1. Open Kaggle (https://www.kaggle.com/) and create a "New Notebook".
#   2. On the right-hand Settings pane, set Accelerator to "GPU T4 x2".
#   3. Install dependencies by running a cell with:
#      !pip install transformers accelerate biopython
#   4. Write prioritized_candidates.faa directly into the workspace:
#      Upload or write from script.
#   5. Run this python script to predict 3D structures.
#   6. Download the resulting folded_structures.zip file.
# =============================================================================

import os
import time
import zipfile
import argparse
import torch
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForProteinFolding

def main():
    parser = argparse.ArgumentParser(description="ESMFold 3D Structure Prediction")
    parser.add_argument("--fasta", default="prioritized_candidates.faa", help="Input FASTA file")
    parser.add_argument("--output", default="pdb_outputs", help="PDB output directory")
    parser.add_argument("--zip", default="folded_structures.zip", help="Output ZIP package name")
    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: Input FASTA '{args.fasta}' not found. Please verify the path.")
        return
        
    os.makedirs(args.output, exist_ok=True)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    if device == "cpu":
        print("WARNING: Running on CPU will be extremely slow. Enable GPU accelerator.")

    print("Loading ESMFold tokenizer and model...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    
    # Cast only the language model stem (model.esm) to float16 on CPU first to save memory
    if device == "cuda":
        model.esm = model.esm.half()
        print("Language model stem cast to float16 precision.")
        
    model = model.to(device)
    
    # Set chunk size to reduce memory usage from O(L^2) to O(L)
    if hasattr(model, "trunk"):
        model.trunk.set_chunk_size(64)
        print("Set ESMFold trunk chunk size to 64.")
        
    model.eval()
    
    print(f"Parsing sequences from {args.fasta}...")
    records = list(SeqIO.parse(args.fasta, "fasta"))
    print(f"Found {len(records)} sequences to fold.")
    
    start_time = time.time()
    for idx, record in enumerate(records):
        seq_id = record.id
        sequence = str(record.seq)
        pdb_out_path = os.path.join(args.output, f"{seq_id}.pdb")
        
        print(f"[{idx+1}/{len(records)}] Folding {seq_id} ({len(sequence)} aa)...")
        
        if len(sequence) > 1000:
            print(f"  Skipping {seq_id} — sequence length exceeds memory limits.")
            continue
            
        try:
            with torch.no_grad():
                inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
                input_ids = inputs["input_ids"].to(device)
                outputs = model(input_ids)
                pdb_string = model.output_to_pdb(outputs)[0]
                
            with open(pdb_out_path, "w") as f:
                f.write(pdb_string)
            print(f"  Saved PDB structure to {pdb_out_path}")
            
        except Exception as e:
            print(f"  Error folding {seq_id}: {str(e)}")
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                
    end_time = time.time()
    print(f"\nFolding complete in {((end_time - start_time)/60):.2f} minutes.")
    
    print(f"Packaging structures into {args.zip}...")
    with zipfile.ZipFile(args.zip, 'w') as zipf:
        for file in os.listdir(args.output):
            if file.endswith(".pdb"):
                zipf.write(os.path.join(args.output, file), file)
                
    print(f"ZIP package created successfully at {args.zip}")

if __name__ == "__main__":
    main()
