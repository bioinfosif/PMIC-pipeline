# PMIC Pipeline

**PMIC Pipeline** is a bioinformatics tool for detecting antimicrobial resistance (AMR) genes and virulence factors (VF) from next-generation sequencing data. It supports both **Illumina paired-end** and **Oxford Nanopore** reads, and performs end-to-end analysis including assembly, species identification, mapping, variant calling, and gene detection on chromosomal and plasmid sequences.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Database Setup](#database-setup)
- [Usage](#usage)
- [Input Format](#input-format)
- [Output Structure](#output-structure)
- [Pipeline Overview](#pipeline-overview)
- [Author](#author)

---

## Features

- Supports **Illumina paired-end** and **Oxford Nanopore** sequencing data
- De novo assembly with **SPAdes** (Illumina) and **Flye** (Nanopore)
- Species identification via **Mash** against the NCBI RefSeq database
- Automatic reference genome download from **NCBI**
- Read mapping with **BWA**, duplicate marking with **Picard**, and variant calling with **bcftools**
- Plasmid identification via **BLAST** against PLSDB
- AMR gene detection using **abricate** with the **CARD** database
- Virulence factor detection using **abricate** with the **VFDB** database
- Automated bar chart generation for AMR/VF gene coverage and identity

---

## Requirements

### System dependencies (via Conda)

| Tool | Purpose |
|---|---|
| Python 3.10 | Runtime |
| SPAdes | Illumina de novo assembly |
| Flye | Nanopore de novo assembly |
| BWA | Read mapping |
| Samtools | BAM processing |
| bcftools | Variant calling |
| Picard | Duplicate marking |
| Mash | Species identification |
| BLAST+ | Plasmid BLAST search |
| abricate | AMR/VF gene detection |

### Python dependencies

```
pandas
numpy
matplotlib
biopython
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/bioinfosif/PMIC-pipeline
cd pmic-pipeline
```

### 2. Create the Conda environment

```bash
conda env create -f environment.yml
conda activate pmic_pipeline
```

### 3. Install the Python package

```bash
pip install -e .
```

This exposes the `pmic-pipeline` command in your environment.

---

## Database Setup

Before running the pipeline for the first time, download the required databases:

```bash
pmic-pipeline --download-databases
```

This will:
1. Download the **Mash RefSeq sketch** (`refseq.genomes.k21s1000.msh`) to the current directory
2. Download the **NCBI plasmid FASTA** and build a **BLAST nucleotide database** under `~/blastdb/plsdb`

> **Note:** You can override the BLAST database directory by setting the `BLASTDB` environment variable before running the pipeline.

```bash
export BLASTDB=/path/to/your/blastdb
pmic-pipeline --download-databases
```

---

## Usage

### Illumina paired-end reads

```bash
pmic-pipeline --illumina -in /path/to/input_dir -out /path/to/output_dir --threads 12
```

### Oxford Nanopore reads

```bash
pmic-pipeline --nanopore -in /path/to/input_dir -out /path/to/output_dir --threads 12
```

### All arguments

| Argument | Description | Default |
|---|---|---|
| `-in`, `--input` | Input directory containing sample files | required |
| `-out`, `--output` | Output directory | required |
| `--illumina` | Process Illumina paired-end reads | — |
| `--nanopore` | Process Oxford Nanopore reads | — |
| `--threads` | Number of CPU threads to use | `12` |
| `--download-databases` | Download Mash and plasmid databases | — |

---

## Input Format

### Illumina

Place all paired-end FASTQ files in the input directory. Files must follow this naming convention:

```
{sample_name}_R1.fastq.gz
{sample_name}_R2.fastq.gz
```

Example:
```
input/
├── sample_A_R1.fastq.gz
├── sample_A_R2.fastq.gz
├── sample_B_R1.fastq.gz
└── sample_B_R2.fastq.gz
```

### Nanopore

Place each sample's FASTQ file in the input directory. The filename (without extension) is used as the sample name:

```
input/
├── sample_A.fastq.gz
└── sample_B.fastq.gz
```

---

## Output Structure

Results are written per sample under the output directory:

```
output/
└── {sample_name}/
    ├── scaffolds.fasta / assembly.fasta       # De novo assembly
    ├── {sample_name}_chr_reference.fasta      # Downloaded reference genome
    ├── {sample_name}_chr_distance.csv         # Mash distance results
    ├── {sample_name}_chr_markdup.bam          # Deduplicated BAM
    ├── {sample_name}_chr_consensus.fasta      # Variant-called consensus
    ├── {sample_name}_chr_AMR.csv              # Chromosomal AMR genes
    ├── {sample_name}_chr_VF.csv               # Chromosomal virulence factors
    ├── {sample_name}_chr_AMR.png              # AMR bar chart
    ├── {sample_name}_chr_VF.png               # VF bar chart
    ├── {sample_name}_plasmid_AMR{N}.csv       # Plasmid AMR genes
    ├── {sample_name}_plasmid_VF{N}.csv        # Plasmid virulence factors
    └── {sample_name}_plasmid_AMR{N}.png       # Plasmid AMR bar chart
```

> Gene plots are only generated when at least 2 genes pass the filters: **coverage ≥ 60%** and **identity ≥ 95%**.

---

## Pipeline Overview

### Illumina workflow

```
Paired FASTQ reads
      │
      ▼
SPAdes plasmid assembly
      │
      ├──► Plasmid contigs ──► BLAST vs PLSDB ──► Download reference
      │                                                   │
      │                                            BWA mapping → bcftools
      │                                                   │
      │                                           abricate (CARD / VFDB)
      │                                                   │
      │                                           AMR / VF plots
      │
      └──► Mash sketch ──► RefSeq distance ──► Download best reference
                                                        │
                                                 BWA mapping → bcftools
                                                        │
                                               abricate (CARD / VFDB)
                                                        │
                                               AMR / VF plots
```

### Nanopore workflow

```
Nanopore FASTQ reads
      │
      ▼
Flye assembly
      │
      ├──► abricate (CARD / VFDB / PlasmidFinder) on full assembly
      │
      ├──► Mash ──► RefSeq species identification
      │
      └──► Plasmid contigs ──► BLAST vs PLSDB
                                     │
                              abricate (CARD / VFDB) per contig
                                     │
                              AMR / VF plots
```

---

## Author

**Mamadou Ndao**  
📧 ndaom403@gmail.com