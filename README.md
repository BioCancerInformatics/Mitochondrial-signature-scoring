# Mitochondrial Signature Scoring

Reproducible R pipeline to calculate **MitoAll** and **MitoOnly** mitochondrial signature scores from bulk transcriptomic data.

This repository accompanies a methodological chapter describing a transcriptome-based approach to estimate mitochondrial abundance indirectly using nuclear-encoded mitochondrial gene signatures.

## Computational framework

The pipeline starts from gene-level TPM expression matrices and produces sample-level mitochondrial signature scores.

![Workflow](https://raw.githubusercontent.com/BioCancerInformatics/Mitochondrial-signature-scoring/main/Results/Figures/Figure1_MitoSignature_Workflow.png)

**Figure. Computational workflow.** GTEx gene-level TPM matrices are imported, gene symbols are harmonized, MitoAll and MitoOnly genes are matched to the expression matrix, expression values are transformed as `log2(TPM + 1)`, gene-wise z-scores are calculated, and sample-level mitochondrial signature scores are obtained by averaging standardized expression across detected genes in each signature.

**Figure 1. Computational framework for mitochondrial signature scoring.** Gene-level TPM matrices are imported, gene symbols are harmonized, MitoAll and MitoOnly genes are matched to the expression matrix, expression values are transformed as `log2(TPM + 1)`, gene-wise z-scores are calculated, and sample-level mitochondrial signature scores are obtained by averaging standardized expression across detected genes in each signature.

## Pipeline

Run the complete workflow in R with:

```r
source("scripts/run_pipeline.R")
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
9. Generate output tables and figures.

## Repository structure

```text
mitochondrial-signature-scoring/
│
├── README.md
├── LICENSE
├── sessionInfo.txt
│
├── data/
│   └── signatures/
│       ├── MitoAll.xlsx
│       └── MitoOnly.xlsx
│
├── scripts/
│   ├── run_pipeline.R
│   ├── 01_import_harmonize_transform.R
│   ├── 02_score_signatures_statistics.R
│   └── 03_generate_framework_figure.R
│
├── results/
│   ├── tables/
│   └── figures/
│       ├── Figure1_MitoSignature_Workflow_README.png
│       ├── Figure1_MitoSignature_Workflow.pdf
│       └── Figure1_MitoSignature_Workflow.tiff
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

Download the gene-level TPM files for the tissues of interest and place them in:

```text
GTEx_data/
```

The mitochondrial signatures used for scoring are provided in:

```text
data/signatures/MitoAll.xlsx
data/signatures/MitoOnly.xlsx
```

Each signature file must contain a column named:

```text
Gene name
```

## Main outputs

Tables are saved in:

```text
results/tables/
```

Figures are saved in:

```text
results/figures/
```

The main framework figure is:

```text
results/figures/Figure1_MitoSignature_Workflow_README.png
```

## Interpretation

MitoAll and MitoOnly scores should be interpreted as transcriptome-based proxies of mitochondrial abundance.

They do not directly measure mitochondrial mass, mitochondrial DNA copy number, mitochondrial membrane potential, respiratory capacity, ATP production, mitophagy, mitochondrial morphology, or mitochondrial quality.

Because this workflow uses bulk transcriptomic data, mitochondrial signature scores may reflect mitochondrial abundance, mitochondrial biogenesis, mitochondrial retrograde signaling, tissue-specific transcriptional programs, and cell-type composition.

## Citation

Citation information will be added upon publication of the chapter.
