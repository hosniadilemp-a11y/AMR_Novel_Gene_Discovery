import os
import urllib.request
import urllib.parse
import json
import gzip
import shutil
import csv
import argparse

def get_st354_assembly_ids():
    print("Searching NCBI for E. coli ST354 assemblies...")
    term = "Escherichia coli[Organism] AND (ST354 OR \"sequence type 354\") AND latest[filter]"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={urllib.parse.quote(term)}&retmax=100&retmode=json"
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            id_list = data.get("esearchresult", {}).get("idlist", [])
            print(f"Found {len(id_list)} assembly IDs.")
            return id_list
    except Exception as e:
        print(f"Error searching NCBI: {e}")
        return []

def get_assembly_ftp_details(id_list):
    if not id_list:
        return {}
    
    print("Fetching assembly summary metadata...")
    id_str = ",".join(id_list)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={id_str}&retmode=json"
    
    details = {}
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            results = data.get("result", {})
            for uid in results.keys():
                if uid == "uids":
                    continue
                acc = results[uid].get("assemblyaccession")
                ftp_gb = results[uid].get("ftppath_genbank")
                ftp_rs = results[uid].get("ftppath_refseq")
                # Prefer GenBank/RefSeq paths
                ftp = ftp_rs if ftp_rs else ftp_gb
                if acc and ftp:
                    details[acc] = ftp
        return details
    except Exception as e:
        print(f"Error fetching metadata: {e}")
        return {}

def download_and_extract_genome(acc, ftp_path, outdir, prefix="ST354"):
    # Convert ftp:// to https://
    https_path = ftp_path.replace("ftp://", "https://")
    basename = os.path.basename(ftp_path)
    fna_gz_name = f"{basename}_genomic.fna.gz"
    download_url = f"{https_path}/{fna_gz_name}"
    
    safe_acc = acc.replace(".", "_")
    target_gz = os.path.join(outdir, f"{prefix}_{safe_acc}.fna.gz")
    target_fna = os.path.join(outdir, f"{prefix}_{safe_acc}.fna")
    
    # Check if already downloaded
    if os.path.exists(target_fna):
        print(f"Genome {acc} already exists. Skipping.")
        return target_fna
        
    print(f"Downloading {acc} from {download_url}...")
    try:
        # Download
        with urllib.request.urlopen(download_url) as response, open(target_gz, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
            
        # Decompress
        print(f"Extracting {target_gz}...")
        with gzip.open(target_gz, 'rb') as f_in, open(target_fna, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
        # Remove gz
        os.remove(target_gz)
        return target_fna
    except Exception as e:
        print(f"Failed to download/extract {acc}: {e}")
        if os.path.exists(target_gz):
            os.remove(target_gz)
        return None

def download_reference(acc, name, outdir):
    print(f"Searching reference assembly details for {acc} ({name})...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={acc}&retmode=json"
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if not id_list:
                print(f"No assembly found for reference {acc}")
                return None
            
            details = get_assembly_ftp_details(id_list)
            if details:
                key = list(details.keys())[0]
                return download_and_extract_genome(key, details[key], outdir, prefix=name)
    except Exception as e:
        print(f"Failed to fetch reference details: {e}")
    return None

def main():
    parser = argparse.ArgumentParser(description="NCBI Genome Batch Downloader")
    parser.add_argument("--accessions", help="Path to text file containing accessions (optional)")
    parser.add_argument("--output", default="results/step3_annotation/reference_genomes", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    metadata = []

    # If accessions file is provided, download those specifically
    if args.accessions and os.path.exists(args.accessions):
        print(f"Loading custom accessions list from: {args.accessions}")
        with open(args.accessions, 'r') as f:
            accessions = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        
        print(f"Downloading {len(accessions)} assemblies...")
        # Since we have specific accessions, we search NCBI for each
        for acc in accessions:
            # Simple name based on accession prefix
            prefix = "REF"
            fna_path = download_reference(acc, prefix, args.output)
            if fna_path:
                metadata.append([acc, prefix, fna_path])
    else:
        # Default: Download ST354 genomes automatically + default references
        REFERENCES = {
            "GCA_000005845.2": "E_coli_K12_MG1655",
            "GCA_000007445.1": "E_coli_CFT073",
            "GCA_000013445.1": "E_coli_UT189",
            "GCA_000008865.2": "E_coli_O157_H7_Sakai"
        }
        
        # 1. Download ST354 Genomes
        st354_ids = get_st354_assembly_ids()
        st354_details = get_assembly_ftp_details(st354_ids)
        
        print(f"Starting download of {len(st354_details)} E. coli ST354 genomes...")
        for acc, ftp in st354_details.items():
            fna_path = download_and_extract_genome(acc, ftp, args.output, prefix="ST354")
            if fna_path:
                metadata.append([acc, "ST354", fna_path])
                
        # 2. Download Reference Outgroups
        print("Starting download of E. coli reference outgroups...")
        for acc, name in REFERENCES.items():
            fna_path = download_reference(acc, name, args.output)
            if fna_path:
                metadata.append([acc, name, fna_path])

    # Write metadata file
    metadata_path = os.path.join(args.output, "assembly_metadata.csv")
    with open(metadata_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Accession", "Group", "FastaPath"])
        writer.writerows(metadata)
        
    print(f"Done! Metadata saved to {metadata_path}. Total downloaded genomes: {len(metadata)}")

if __name__ == "__main__":
    main()
