# Genome-Agnostic Classification of RNA Virus Infections
Machine learning on host transcriptomes for virus detection without viral sequences.

This repository contains code and data processing pipelines for the manuscript:
**"Machine Learning Enables Genome-Agnostic Classification of RNA Virus Infections from Host Transcriptomes."**

We implement workflows to:
- Preprocess RNA-seq data from public repositories (SRA/GEO).
- Generate host-based features (gene expression & k-mer spectra).
- Run unsupervised clustering and Random Forest classifiers.
- Evaluate robustness with shuffled, permuted, depleted, and enriched datasets.

## Repository Structure
- `notebooks/` – interactive exploration, data checks, prototyping.  
- `scripts/` – reproducible workflows (data preprocessing, model training, figure generation).  
- `data/` – input/output data; large raw files should be excluded via `.gitignore`.  
- `figures/` – curated figures for the manuscript.
- `assemblies/` - files used for assembly in Geneious prime to create enriched or depleted datasets  
- **LICENSE** – open-source license.  
- **README.md** – documentation for installation, usage, and repo overview.

## Analysis workflow

All analysis scripts are in `scripts/` and can be visualized using code from various Jupyter notebooks in `notebooks/`.

`pipeline.py` can be used to generate kmer frequency spectra and run hierarchical clustering and random forest classification from .fastq or .fasta files given an associated metadata file that has labels corresponding to each sequencing sample.

`shuffle-permute.py` takes a .fastq file and creates shuffled and permuted parallel datasets which can be used as controls.

`figures.ipynb` takes the output data and a metadata file in .csv format to generate all of the figures in the paper and supplement.

`scratch-notes.ipynb` visualizes more data and has notes on other analyses that were not used in this study.

## Enrichment and Depletion

The viral reads were enriched or depleted using Geneious Prime 2023.1.2 using Geneious Assembler to map all reads in a .fastq file to a set of viral reference genomes.

This repo does **not** host full RNA-seq data due to size limits.  
Raw data can be downloaded from SRA (PRJNA1074963) and GEO (GSE49840).  
Metadata tables are provided in `data/` and the files used in the enrichment/depletion assemblies are in `assemblies/`.

If you use this code, please cite:
Benjamin Kaza et al. (2025). Machine Learning Enables Genome-Agnostic Classification of RNA Virus Infections from Host Transcriptomes.

MIT License. See LICENSE file for details.
