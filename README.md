# Mitochondrial Signature Scoring

Reproducible R workflow to calculate **MitoAll** and **MitoOnly** mitochondrial signature scores from bulk transcriptomic data.

This repository accompanies a methodological chapter describing a simple computational approach to estimate mitochondrial abundance indirectly from transcriptomic datasets using nuclear-encoded mitochondrial gene signatures.


## Scientific rationale

Mitochondria contain their own genome, but the vast majority of mitochondrial proteins are encoded by nuclear genes. Therefore, the expression of nuclear-encoded mitochondrial genes can provide a transcriptome-based proxy of the mitochondrial compartment.

This workflow uses two mitochondrial signatures:

- **MitoAll**: genes encoding proteins with recognized mitochondrial localization, including proteins that may also localize to other cellular compartments.
- **MitoOnly**: a stricter subset of genes encoding proteins recognized as mitochondrial-only.

The representative analysis uses **GTEx Analysis V11** normal tissue RNA-seq data to compare mitochondrial signature scores across tissues expected to differ in mitochondrial abundance.


## Workflow

The workflow starts from gene-level TPM expression matrices and produces sample-level mitochondrial signature scores.

![Workflow](results/figures/Figure0_Workflow.png)

**Figure 0. Computational workflow.** GTEx gene-level TPM matrices are imported, gene symbols are harmonized, MitoAll and MitoOnly genes are matched to the expression matrix, expression values are transformed as `log2(TPM + 1)`, gene-wise z-scores are calculated, and sample-level mitochondrial signature scores are obtained by averaging standardized expression across detected genes in each signature.


## Representative analysis

The representative analysis compares four GTEx normal tissues:

- Heart - Atrial Appendage
- Heart - Left Ventricle
- Muscle - Skeletal
- Whole Blood

Cardiac and skeletal muscle tissues were selected as mitochondria-rich reference tissues, whereas Whole Blood was selected as a lower-mitochondrial comparator.

![Representative result](results/figures/Figure1_MitoSignatures.png)

**Figure 1. Representative mitochondrial signature scoring in GTEx normal tissues.** Distribution of sample-level MitoAll and MitoOnly scores across selected GTEx tissues. Box plots depict median and interquartile range; whiskers identify the most extreme non-outlier observations according to standard Tukey fences; individual points represent GTEx samples.


## Repository structure

```text
mitochondrial-signature-scoring/
│
├── README.md
├── LICENSE
├── sessionInfo.txt
│
├── data/
│   ├── signatures/
│   │   ├── MitoAll.xlsx
│   │   └── MitoOnly.xlsx
│   │
│   └── metadata/
│       └── gtex_tissues_used.csv
│
├── scripts/
│   ├── 01_import_harmonize_transform.R
│   ├── 02_score_signatures_statistics.R
│   └── 03_generate_representative_figure.R
│
├── results/
│   ├── tables/
│   │   ├── sample_summary_by_tissue.csv
│   │   ├── signature_coverage_summary.csv
│   │   ├── signature_score_summary_by_tissue.csv
│   │   └── wilcoxon_planned_comparisons_fdr.csv
│   │
│   └── figures/
│       ├── Figure0_Workflow.png
│       ├── Figure0_Workflow.pdf
│       ├── Figure1_MitoSignatures.png
│       ├── Figure1_MitoSignatures.pdf
│       └── Figure1_MitoSignatures.tiff
│
└── docs/
    └── GTEx_download_instructions.md
```


## Required software

The workflow was developed in **R**. The following R packages are required:

```r
data.table
dplyr
tidyr
purrr
stringr
readr
rio
ggplot2
patchwork
scales
```

The computational environment used for the representative analysis is provided in:

```text
sessionInfo.txt
```


## Input data

GTEx expression matrices are not included in this repository because of file size.

Download the following **GTEx Analysis V11 gene-level TPM files** and place them in a local folder named `GTEx_data/`:

```text
gene_tpm_v11_heart_atrial_appendage.gct.gz
gene_tpm_v11_heart_left_ventricle.gct.gz
gene_tpm_v11_muscle_skeletal.gct.gz
gene_tpm_v11_whole_blood.gct.gz
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


## How to run

Run the scripts in numerical order:

```r
source("scripts/01_import_harmonize_transform.R")
source("scripts/02_score_signatures_statistics.R")
source("scripts/03_generate_representative_figure.R")
```

The scripts perform the following steps:

1. Import GTEx gene-level TPM matrices.
2. Harmonize gene symbols.
3. Match MitoAll and MitoOnly genes to the expression matrix.
4. Transform expression values as `log2(TPM + 1)`.
5. Calculate gene-wise z-scores.
6. Calculate sample-level MitoAll and MitoOnly scores.
7. Assess signature gene recovery.
8. Compare tissues using Wilcoxon rank-sum tests.
9. Adjust P values using the Benjamini-Hochberg false discovery rate method.
10. Generate the representative figure.


## Main outputs

The main sample-level output is:

```text
results/tables/signature_score_summary_by_tissue.csv
```

Additional summary tables include:

```text
results/tables/sample_summary_by_tissue.csv
results/tables/signature_coverage_summary.csv
results/tables/wilcoxon_planned_comparisons_fdr.csv
```

The main figure outputs are:

```text
results/figures/Figure0_Workflow.png
results/figures/Figure1_MitoSignatures.png
```

PDF and TIFF versions are also provided for publication use.


## Interpretation

MitoAll and MitoOnly scores should be interpreted as transcriptome-based proxies of mitochondrial abundance.

They do not directly measure:

- mitochondrial mass;
- mitochondrial DNA copy number;
- mitochondrial membrane potential;
- respiratory capacity;
- ATP production;
- mitophagy;
- mitochondrial morphology;
- mitochondrial quality.

Because this workflow uses bulk tissue transcriptomic data, mitochondrial signature scores may reflect mitochondrial abundance, mitochondrial biogenesis, mitochondrial retrograde signaling, tissue-specific transcriptional programs, and cell-type composition.


## Reproducibility notes

- Use gene-level expression matrices, not transcript-level matrices.
- Keep GTEx sample identifiers unchanged.
- Check gene-symbol compatibility before scoring.
- Always report signature coverage.
- Do not compare scores generated from different tissue panels unless the same z-score reference set was used.
- Inspect both effect sizes and score distributions, not only adjusted P values.
- GTEx expression files are large and should be downloaded directly from the GTEx Portal rather than stored in this repository.


## Citation

Citation information will be added upon publication of the chapter.
