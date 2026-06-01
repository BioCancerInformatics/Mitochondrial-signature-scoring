# Mitochondrial Signature Scoring

### *Interrogating sequencing data for transcriptional signatures of mitochondrial abundance*

A reproducible R workflow accompanying the methodological chapter:

**Higor Almeida Cordeiro Nogueira¹², Abdelrahman M. Ghazal¹, Spencer Brackett¹, Flavie Naulin¹, Ai Sato¹, Pietro Mancuso¹, Manuel Beltrán-Visiedo¹, Emma Guilbaud¹, Enrique Medina-Acosta², Lorenzo Galluzzi¹**, and **Lukas Bolini¹**

¹ Cancer Signaling and Microenvironment Program, Fox Chase Cancer Center, Philadelphia, PA, USA
² Laboratório de Biotecnologia, Centro de Biociências e Biotecnologia, Universidade Estadual do Norte Fluminense Darcy Ribeiro, Campos dos Goytacazes, RJ, Brazil

**Correspondence:** Lorenzo Galluzzi ([deadoc80@gmail.com](mailto:deadoc80@gmail.com)) and Lukas Bolini ([lukas.goncalves@fccc.edu](mailto:lukas.goncalves@fccc.edu))


## Overview

This repository provides a reproducible R pipeline to calculate **MitoAll** and **MitoOnly** mitochondrial signature scores from bulk transcriptomic data.

The workflow was developed to support a methodological chapter describing a transcriptome-based approach to estimate mitochondrial abundance indirectly using nuclear-encoded mitochondrial gene signatures.

Mitochondria contain their own genome, but most mitochondrial proteins are encoded by nuclear genes. Therefore, the expression of nuclear-encoded mitochondrial genes can be used as a transcriptomic proxy of the mitochondrial compartment.

This workflow scores two mitochondrial signatures:

* **MitoAll**: genes encoding proteins with recognized mitochondrial localization, including proteins that may also localize to other cellular compartments.
* **MitoOnly**: a stricter subset of genes encoding proteins recognized as mitochondrial-only.


## Computational framework

The pipeline starts from gene-level TPM expression matrices and produces sample-level mitochondrial signature scores.

![Mitochondrial signature scoring workflow](Results/Figures/Figure1_MitoSignature_Workflow.png)

**Computational workflow.** Gene-level TPM matrices are imported, gene symbols are harmonized, MitoAll and MitoOnly genes are matched to the expression matrix, expression values are transformed as `log2(TPM + 1)`, gene-wise z-scores are calculated, and sample-level mitochondrial signature scores are obtained by averaging standardized expression across detected genes in each signature.

## Representative result

The representative analysis applies the mitochondrial signature scoring workflow to selected GTEx normal tissues. The resulting figure shows the distribution of **MitoAll** and **MitoOnly** sample-level scores across tissues with distinct expected mitochondrial abundance.

![Figure 1. Representative mitochondrial signature scoring results](Results/Figures/Figure_1.png)

**Figure 1. Representative mitochondrial signature scoring in GTEx normal tissues.** Sample-level MitoAll and MitoOnly scores were calculated after gene-symbol harmonization, `log2(TPM + 1)` transformation, and gene-wise z-score standardization. The figure illustrates how the workflow can be used to compare mitochondrial transcriptional signatures across normal human tissues.

## Pipeline

Run the complete workflow in R with:

```r
source("Rscript/run_pipeline.R")
```

The pipeline performs the following steps:

1. Import gene-level TPM matrices.
2. Harmonize gene symbols.
3. Match MitoAll and MitoOnly genes to the expression matrix.
4. Transform expression values as `log2(TPM + 1)`.
5. Calculate gene-wise z-scores.
6. Calculate sample-level MitoAll and MitoOnly scores.
7. Assess signature gene recovery.
8. Perform statistical comparisons when applicable.
9. Generate output tables, RDS objects, session information, and publication-ready figures.


## Repository structure

```text
Mitochondrial-signature-scoring/
│
├── README.md
│
├── signatures/
│   ├── MitoAll.xlsx
│   └── MitoOnly.xlsx
│
├── Rscript/
│   └── run_pipeline.R
│
├── Results/
│   ├── Tables/
│   │   ├── sample_summary_by_tissue.csv
│   │   ├── detected_signature_genes_by_tissue.csv
│   │   ├── signature_coverage_summary_import.csv
│   │   ├── signature_coverage_by_tissue.csv
│   │   ├── signature_coverage_overall.csv
│   │   ├── gtex_mitoall_mitoonly_sample_scores.csv
│   │   ├── signature_score_summary_by_tissue.csv
│   │   ├── wilcoxon_planned_comparisons_fdr.csv
│   │   └── Figure_1_score_summary_by_tissue.csv
│   │
│   ├── RDS/
│   │   ├── signature_genes_clean.rds
│   │   ├── gtex_mito_signature_genes_log2TPM_long.rds
│   │   ├── gtex_expr_gene_wise_zscores.rds
│   │   └── gtex_mitoall_mitoonly_sample_scores.rds
│   │
│   ├── Figures/
│   │   ├── Figure1_MitoSignature_Workflow.png
│   │   ├── Figure1_MitoSignature_Workflow.pdf
│   │   ├── Figure_1.png
│   │   ├── Figure_1.pdf
│   │   └── Figure_1.tiff
│   │
│   └── SessionInfo/
│       └── run_pipeline_sessionInfo.txt
│
└── docs/
    └── GTEx_download_instructions.md
```


## Required R packages

```r
data.table
dplyr
tidyr
purrr
stringr
readr
rio
ggplot2
cowplot
grid
scales
```


## Input data

GTEx expression files are not included in this repository because of file size.

Download the gene-level TPM files for the tissues of interest and place them in a local folder named:

```text
GTEx_data/
```

The mitochondrial signatures used for scoring are provided in:

```text
signatures/MitoAll.xlsx
signatures/MitoOnly.xlsx
```

Each signature file must contain a column named:

```text
Gene name
```


## Main outputs

Tables are saved in:

```text
Results/Tables/
```

Figures are saved in:

```text
Results/Figures/
```

The main framework figure is:

```text
Results/Figures/Figure_MitoSignature_Workflow.png
```

The representative result figure is:

```text
Results/Figures/Figure_1.png
Results/Figures/Figure_1.pdf
```

## Interpretation

MitoAll and MitoOnly scores should be interpreted as transcriptome-based proxies of mitochondrial abundance.

They do not directly measure mitochondrial mass, mitochondrial DNA copy number, mitochondrial membrane potential, respiratory capacity, ATP production, mitophagy, mitochondrial morphology, or mitochondrial quality.

Because this workflow uses bulk transcriptomic data, mitochondrial signature scores may reflect mitochondrial abundance, mitochondrial biogenesis, mitochondrial retrograde signaling, tissue-specific transcriptional programs, and cell-type composition.


## Citation

Citation information will be added upon publication of the chapter.



