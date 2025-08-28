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
- **notebooks/** – interactive exploration, data checks, prototyping.  
- **scripts/** – reproducible workflows (data preprocessing, model training, figure generation).  
- **data/** – input/output data; large raw files should be excluded via `.gitignore`.  
- **results/** – auto-generated results; not manually edited.  
- **fig/** – curated figures for the manuscript.  
- **env/** – reproducible environment definitions.  
- **LICENSE** – open-source license.  
- **README.md** – documentation for installation, usage, and repo overview.  

This repo does **not** host full RNA-seq data due to size limits.  
Raw data can be downloaded from SRA (PRJNA1074963) and GEO (GSE49840).  
Metadata tables are provided in `data/` to link accession numbers with conditions.

If you use this code, please cite:
Benjamin Kaza et al. (2025). Machine Learning Enables Genome-Agnostic Classification of RNA Virus Infections from Host Transcriptomes.

MIT License. See LICENSE file for details.
