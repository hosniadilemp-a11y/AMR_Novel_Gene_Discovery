#!/usr/bin/env python3
# =============================================================================
# 07_plm_embeddings.py — Protein Language Model (ESM-2) Embedding Extraction
# =============================================================================
#
# USAGE:
#   python3 scripts/07_plm_embeddings.py [OPTIONS]
#
# =============================================================================
import os
import sys
import time
import json
import csv
import argparse
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from transformers import AutoTokenizer, AutoModel
import umap

def read_fasta(fasta_path):
    records = []
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            records.append((current_id, "".join(current_seq)))
    return records

def main():
    parser = argparse.ArgumentParser(description="ESM-2 Protein Language Model Embeddings")
    parser.add_argument("--candidates", required=True, help="Path to novel prioritized candidates FASTA")
    parser.add_argument("--reference", required=True, help="Path to known acquired AMR genes FASTA")
    parser.add_argument("--model", default="facebook/esm2_t33_650M_UR50D", help="ESM-2 model identifier")
    parser.add_argument("--output", default="results/step7_plm", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    amr_locus_tags = {
        "KNGPFPPJ_04581": "Tetracycline [tet(A)]",
        "KNGPFPPJ_04687": "Macrolide [estX]",
        "KNGPFPPJ_04689": "Aminoglycoside [aadA2]",
        "KNGPFPPJ_04690": "Phenicol [cmlA1]",
        "KNGPFPPJ_04691": "Aminoglycoside [aadA1]",
        "KNGPFPPJ_04692": "Biocide [qacL]",
        "KNGPFPPJ_04694": "Sulfonamide [sul3]",
        "KNGPFPPJ_04716": "Trimethoprim [dfrA7]",
        "KNGPFPPJ_04717": "Biocide [qacE\\Delta1]",
        "KNGPFPPJ_04718": "Sulfonamide [sul1]",
        "KNGPFPPJ_04720": "Beta-lactam [blaTEM-1B]",
        "KNGPFPPJ_04727": "Aminoglycoside [aph(3')-Ia]"
    }

    print("Reading FASTA sequences...")
    candidates = read_fasta(args.candidates)
    known_amrs = read_fasta(args.reference)

    all_data = []
    # Add candidates
    for seq_id, seq in candidates:
        all_data.append({
            "id": seq_id,
            "sequence": seq,
            "group": "Novel Candidate",
            "label": seq_id,
            "is_candidate": True
        })
    # Add known AMR genes
    for seq_id, seq in known_amrs:
        label = amr_locus_tags.get(seq_id, seq_id)
        all_data.append({
            "id": seq_id,
            "sequence": seq,
            "group": label,
            "label": label.split(" ")[0],
            "is_candidate": False
        })

    print(f"Total sequences to embed: {len(all_data)}")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    print(f"Loading ESM-2 model: {args.model}...")
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = AutoModel.from_pretrained(args.model)
    model = model.to(device)
    model.eval()

    embeddings = []
    for idx, item in enumerate(all_data):
        seq = item["sequence"]
        seq_id = item["id"]
        
        inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            last_hidden = outputs.last_hidden_state.squeeze(0).cpu().numpy()
            mean_embedding = np.mean(last_hidden[1:-1, :], axis=0)
            embeddings.append(mean_embedding)
            
        if (idx+1) % 5 == 0 or (idx+1) == len(all_data):
            print(f"  [{idx+1}/{len(all_data)}] Embedded {seq_id}")

    embeddings = np.array(embeddings)
    print(f"Embeddings shape: {embeddings.shape}")
    np.save(os.path.join(args.output, "esm2_embeddings.npy"), embeddings)

    # Dimensionality Reduction
    print("Running t-SNE reduction...")
    tsne = TSNE(n_components=2, perplexity=7, random_state=42, n_iter=2000)
    coords_tsne = tsne.fit_transform(embeddings)

    print("Running UMAP reduction...")
    reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, n_components=2, random_state=42)
    coords_umap = reducer.fit_transform(embeddings)

    df_coords = pd.DataFrame({
        "ID": [item["id"] for item in all_data],
        "Label": [item["label"] for item in all_data],
        "Group": [item["group"] for item in all_data],
        "Is_Candidate": [item["is_candidate"] for item in all_data],
        "tSNE_1": coords_tsne[:, 0],
        "tSNE_2": coords_tsne[:, 1],
        "UMAP_1": coords_umap[:, 0],
        "UMAP_2": coords_umap[:, 1]
    })

    coords_path = os.path.join(args.output, "esm2_plm_coordinates.tsv")
    df_coords.to_csv(coords_path, sep="\t", index=False)
    print(f"Saved coordinates to {coords_path}")

    # Compute Euclidean distance matrix
    from scipy.spatial.distance import pdist, squareform
    dist_matrix = squareform(pdist(embeddings))
    dist_df = pd.DataFrame(dist_matrix, index=df_coords["ID"], columns=df_coords["ID"])
    dist_path = os.path.join(args.output, "plm_distance_matrix.tsv")
    dist_df.to_csv(dist_path, sep="\t")
    print(f"Distance matrix saved to {dist_path}")

if __name__ == "__main__":
    main()
