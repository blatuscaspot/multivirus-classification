# Genome-Agnostic Classification of RNA Virus Infections
Machine learning on host transcriptomes for virus detection without viral sequences.

This repository contains code and data processing pipelines for the manuscript:
**"Machine Learning Enables Genome-Agnostic Classification of RNA Virus Infections from Host Transcriptomes."**

We implement workflows to:
- Preprocess RNA-seq data from public repositories (SRA/GEO).
- Generate host-based features (gene expression & k-mer spectra).
- Run unsupervised clustering and Random Forest classifiers.
- Evaluate robustness with shuffled, permuted, depleted, and enriched datasets.

## Installation
Clone this repository and create the environment:

```bash
git clone https://github.com/blatuscaspot/iav-rna.git
cd virus-transcriptome-ML
conda env create -f environment.yml
conda activate virusml

---

Usage
Give simple example commands:
```markdown
## Usage

### Preprocess raw FASTQ files
```bash
python scripts/preprocess_fastq.py --input data/SRRXXXXXX.fastq.gz --output results/sample_expression.csv

python scripts/train_rf.py --features results/sample_expression.csv --labels data/metadata.csv

This repo does **not** host full RNA-seq data due to size limits.  
Raw data can be downloaded from SRA (PRJNA1074963) and GEO (GSE49840).  
Metadata tables are provided in `data/` to link accession numbers with conditions.

If you use this code, please cite:
Benjamin Kaza et al. (2025). Machine Learning Enables Genome-Agnostic Classification of RNA Virus Infections from Host Transcriptomes.

MIT License. See LICENSE file for details.
