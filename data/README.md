# Data Directory Overview

This directory contains instructions and metadata lists for obtaining the sequence data used in this study.

## 1. Raw Sequencing Reads

The raw paired-end Illumina MiSeq reads for clinical isolate **QA5221** are archived at the NCBI Sequence Read Archive (SRA):

* **BioProject Accession:** `PRJNA1481519`
* **SRA Run Accession:** `SRR39314025`

### Downloading Reads

You can retrieve the FASTQ files using NCBI SRA Toolkit's `fasterq-dump` (automatically included in the `amr_env` Conda environment):

```bash
# Download reads to raw_reads/ folder
fasterq-dump SRR39314025 --outdir raw_reads/ --split-files

# Compress fastq files
gzip raw_reads/SRR39314025_1.fastq
gzip raw_reads/SRR39314025_2.fastq

# Rename to match standard pipeline conventions
mv raw_reads/SRR39314025_1.fastq.gz raw_reads/QA5221_R1.fastq.gz
mv raw_reads/SRR39314025_2.fastq.gz raw_reads/QA5221_R2.fastq.gz
```

## 2. Reference Genome Accessions

The NCBI GenBank/RefSeq accessions for the 32 ST354 reference cohort and outgroup comparison genomes are listed in:

* `data/reference_accessions/st354_cohort_accessions.txt`
* `data/reference_accessions/outgroup_accessions.txt`

These files are read automatically by `scripts/03b_pangenome.sh` to download and annotate the required reference genomes.
