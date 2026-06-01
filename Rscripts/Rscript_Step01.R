###############################################################################
# GTEx MITOCHONDRIAL SIGNATURE ANALYSIS
# Step 1 — Import files, harmonize gene symbols, transform expression
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(rio)
})

###############################################################################
# 0) USER SETTINGS
###############################################################################

# Change these paths according to your computer
BASE_DIR <- "C:/Users/CORDEIH/Desktop/Whole_Project/Manuscript/Analyses"

EXPR_DIR <- file.path(BASE_DIR, "GTEx_data")
SIG_DIR  <- file.path(BASE_DIR, "Mitocondrial_signatures")
OUT_DIR  <- file.path(BASE_DIR, "output")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Signature files
# These can be .xlsx, .csv, .tsv, or .rds if rio::import can read them
mitoall_file  <- file.path(SIG_DIR, "MitoAll.xlsx")
mitoonly_file <- file.path(SIG_DIR, "MitoOnly.xlsx")

# GTEx tissue TPM files already downloaded
# Edit the file names if yours are slightly different
gtex_files <- tibble::tribble(
  ~tissue, ~file,
  "Heart - Left Ventricle",
  file.path(EXPR_DIR, "gene_tpm_v11_heart_left_ventricle.gct/gene_tpm_v11_heart_left_ventricle.gct"),
  
  "Heart - Atrial Appendage",
  file.path(EXPR_DIR, "gene_tpm_v11_heart_atrial_appendage.gct/gene_tpm_v11_heart_atrial_appendage.gct"),
  
  "Muscle - Skeletal",
  file.path(EXPR_DIR, "gene_tpm_v11_muscle_skeletal.gct/gene_tpm_v11_muscle_skeletal.gct"),
  
  # "Spleen",
  # file.path(EXPR_DIR, "gene_tpm_v11_spleen.gct/gene_tpm_v11_spleen.gct"),
  
  "Whole Blood",
  file.path(EXPR_DIR, "gene_tpm_v11_whole_blood.gct/gene_tpm_v11_whole_blood.gct")
)

###############################################################################
# 1) HELPER FUNCTIONS
###############################################################################

clean_gene_symbol <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_upper() %>%
    na_if("") %>%
    na_if("NA")
}

check_files_exist <- function(paths) {
  missing_files <- paths[!file.exists(paths)]
  
  if (length(missing_files) > 0) {
    stop(
      "The following files were not found:\n",
      paste(missing_files, collapse = "\n")
    )
  }
}

import_signature <- function(file, signature_name) {
  
  df <- rio::import(file)
  
  if (!"Gene name" %in% colnames(df)) {
    stop(
      "The file ", basename(file),
      " does not contain a column named 'Gene name'."
    )
  }
  
  df %>%
    transmute(
      signature = signature_name,
      gene_symbol = clean_gene_symbol(`Gene name`)
    ) %>%
    filter(!is.na(gene_symbol)) %>%
    distinct(signature, gene_symbol)
}

read_gtex_tpm_gct <- function(file, tissue, target_genes = NULL) {
  
  message("Reading: ", basename(file))
  
  # GTEx GCT files have two metadata lines:
  # line 1: version
  # line 2: matrix dimensions
  # Therefore, skip = 2 is necessary.
  gct <- data.table::fread(file, skip = 2, data.table = FALSE)
  
  if (!all(c("Name", "Description") %in% colnames(gct))) {
    stop(
      "The file ", basename(file),
      " does not look like a standard GTEx GCT file."
    )
  }
  
  sample_cols <- setdiff(colnames(gct), c("Name", "Description"))
  
  gct_clean <- gct %>%
    mutate(
      ensembl_id = str_remove(Name, "\\..*$"),
      gene_symbol = clean_gene_symbol(Description)
    ) %>%
    filter(!is.na(gene_symbol))
  
  # Optional but recommended:
  # keep only genes that belong to MitoAll or MitoOnly.
  # This makes the analysis much lighter.
  if (!is.null(target_genes)) {
    gct_clean <- gct_clean %>%
      filter(gene_symbol %in% target_genes)
  }
  
  # If duplicated gene symbols exist, collapse them.
  # For TPM, using the mean is conservative and avoids artificially inflating
  # expression by summing duplicated symbols.
  gct_long <- gct_clean %>%
    select(gene_symbol, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample_id",
      values_to = "TPM"
    ) %>%
    mutate(
      TPM = as.numeric(TPM),
      tissue = tissue
    ) %>%
    group_by(tissue, sample_id, gene_symbol) %>%
    summarise(
      TPM = mean(TPM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log2_TPM = log2(TPM + 1)
    )
  
  return(gct_long)
}

###############################################################################
# 2) CHECK INPUT FILES
###############################################################################

check_files_exist(c(mitoall_file, mitoonly_file, gtex_files$file))

# ###############################################################################
# # 3) IMPORT MITOALL AND MITOONLY SIGNATURES
# ###############################################################################
# 
# mitoall_genes <- import_signature(
#   file = mitoall_file,
#   signature_name = "MitoAll"
# )
# 
# mitoonly_genes <- import_signature(
#   file = mitoonly_file,
#   signature_name = "MitoOnly"
# )
# 
# signature_genes <- bind_rows(
#   mitoall_genes,
#   mitoonly_genes
# )
# 
# target_genes <- signature_genes %>%
#   pull(gene_symbol) %>%
#   unique()
# 
# message("MitoAll genes: ", n_distinct(mitoall_genes$gene_symbol))
# message("MitoOnly genes: ", n_distinct(mitoonly_genes$gene_symbol))
# message("Total unique mitochondrial signature genes: ", length(target_genes))


###############################################################################
# 3) IMPORT MITOALL, MITOONLY, AND mtRNA SIGNATURES
###############################################################################

mitoall_genes <- import_signature(
  file = mitoall_file,
  signature_name = "MitoAll"
)

mitoonly_genes <- import_signature(
  file = mitoonly_file,
  signature_name = "MitoOnly"
)

# 13 mitochondrial-encoded protein-coding genes
# These genes are encoded by the mitochondrial genome and can be used to build
# an mtRNA expression score for comparison with nuclear mitochondrial signatures.

mtrna_genes <- tibble::tibble(
  signature = "mtRNA_13PCG",
  gene_symbol = c(
    "MT-ND1",
    "MT-ND2",
    "MT-ND3",
    "MT-ND4",
    "MT-ND4L",
    "MT-ND5",
    "MT-ND6",
    "MT-CO1",
    "MT-CO2",
    "MT-CO3",
    "MT-CYB",
    "MT-ATP6",
    "MT-ATP8"
  )
) %>%
  mutate(
    gene_symbol = clean_gene_symbol(gene_symbol)
  ) %>%
  distinct(signature, gene_symbol)

signature_genes <- bind_rows(
  mitoall_genes,
  mitoonly_genes,
  mtrna_genes
)

target_genes <- signature_genes %>%
  pull(gene_symbol) %>%
  unique()

message("MitoAll genes: ", n_distinct(mitoall_genes$gene_symbol))
message("MitoOnly genes: ", n_distinct(mitoonly_genes$gene_symbol))
message("mtRNA_13PCG genes: ", n_distinct(mtrna_genes$gene_symbol))
message("Total unique mitochondrial signature genes: ", length(target_genes))

###############################################################################
# 4) IMPORT GTEx TPM FILES AND TRANSFORM EXPRESSION
###############################################################################

gtex_expr_long <- pmap_dfr(
  .l = list(
    file = gtex_files$file,
    tissue = gtex_files$tissue
  ),
  .f = ~ read_gtex_tpm_gct(
    file = ..1,
    tissue = ..2,
    target_genes = target_genes
  )
)

###############################################################################
# 5) BASIC QUALITY CHECKS
###############################################################################

# Number of samples per tissue
sample_summary <- gtex_expr_long %>%
  distinct(tissue, sample_id) %>%
  count(tissue, name = "n_samples") %>%
  arrange(tissue)

# Number of detected signature genes per tissue
gene_summary <- gtex_expr_long %>%
  distinct(tissue, gene_symbol) %>%
  count(tissue, name = "n_signature_genes_detected") %>%
  arrange(tissue)

# Signature coverage
signature_coverage <- signature_genes %>%
  left_join(
    gtex_expr_long %>%
      distinct(gene_symbol) %>%
      mutate(detected_in_gtex = TRUE),
    by = "gene_symbol"
  ) %>%
  mutate(detected_in_gtex = if_else(is.na(detected_in_gtex), FALSE, detected_in_gtex)) %>%
  group_by(signature) %>%
  summarise(
    genes_in_signature = n_distinct(gene_symbol),
    genes_detected_in_gtex = n_distinct(gene_symbol[detected_in_gtex]),
    coverage = genes_detected_in_gtex / genes_in_signature,
    .groups = "drop"
  )

###############################################################################
# 6) EXPORT CLEANED OBJECTS
###############################################################################

saveRDS(
  signature_genes,
  file.path(OUT_DIR, "signature_genes_clean.rds")
)

saveRDS(
  gtex_expr_long,
  file.path(OUT_DIR, "gtex_mito_signature_genes_log2TPM_long.rds")
)

write_csv(
  sample_summary,
  file.path(OUT_DIR, "sample_summary_by_tissue.csv")
)

write_csv(
  gene_summary,
  file.path(OUT_DIR, "detected_signature_genes_by_tissue.csv")
)

write_csv(
  signature_coverage,
  file.path(OUT_DIR, "signature_coverage_summary.csv")
)

message("Done.")
message("Main object saved as: gtex_mito_signature_genes_log2TPM_long.rds")

