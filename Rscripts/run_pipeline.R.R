###############################################################################
# MITOCHONDRIAL SIGNATURE SCORING PIPELINE
# MitoAll and MitoOnly scoring from GTEx gene-level TPM matrices
#
# This script performs:
# 1. Import of mitochondrial signatures
# 2. Import of GTEx gene-level TPM matrices
# 3. Gene-symbol harmonization
# 4. log2(TPM + 1) transformation
# 5. Gene-wise z-score calculation
# 6. Sample-level MitoAll and MitoOnly scoring
# 7. Signature coverage assessment
# 8. Wilcoxon rank-sum tests with FDR correction
# 9. Export of reproducible tables and R objects
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
# Run this script from the repository root:
# source("Rscript/run_pipeline.R")

BASE_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

EXPR_DIR <- file.path(BASE_DIR, "GTEx_data")
SIG_DIR  <- file.path(BASE_DIR, "signatures")

RESULTS_DIR <- file.path(BASE_DIR, "Results")
TABLE_DIR   <- file.path(RESULTS_DIR, "Tables")
RDS_DIR     <- file.path(RESULTS_DIR, "RDS")
INFO_DIR    <- file.path(RESULTS_DIR, "SessionInfo")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INFO_DIR, recursive = TRUE, showWarnings = FALSE)

mitoall_file  <- file.path(SIG_DIR, "MitoAll.xlsx")
mitoonly_file <- file.path(SIG_DIR, "MitoOnly.xlsx")

gtex_files <- tibble::tribble(
  ~tissue, ~file_name,
  "Heart - Left Ventricle",   "gene_tpm_v11_heart_left_ventricle.gct.gz",
  "Heart - Atrial Appendage", "gene_tpm_v11_heart_atrial_appendage.gct.gz",
  "Muscle - Skeletal",        "gene_tpm_v11_muscle_skeletal.gct.gz",
  "Whole Blood",              "gene_tpm_v11_whole_blood.gct.gz"
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

resolve_gtex_file <- function(file_name, expr_dir) {
  
  uncompressed_name <- str_remove(file_name, "\\.gz$")
  
  candidate_paths <- c(
    file.path(expr_dir, file_name),
    file.path(expr_dir, uncompressed_name),
    file.path(expr_dir, uncompressed_name, basename(uncompressed_name))
  )
  
  existing_path <- candidate_paths[file.exists(candidate_paths)]
  
  if (length(existing_path) == 0) {
    stop(
      "Could not find GTEx file for: ", file_name, "\n",
      "Checked:\n",
      paste(candidate_paths, collapse = "\n")
    )
  }
  
  existing_path[1]
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
  
  message("Reading GTEx file: ", basename(file), " | Tissue: ", tissue)
  
  gct <- data.table::fread(
    file,
    skip = 2,
    data.table = FALSE,
    showProgress = FALSE
  )
  
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
  
  if (!is.null(target_genes)) {
    gct_clean <- gct_clean %>%
      filter(gene_symbol %in% target_genes)
  }
  
  gct_long <- gct_clean %>%
    select(gene_symbol, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample_id",
      values_to = "TPM"
    ) %>%
    mutate(
      tissue = tissue,
      TPM = as.numeric(TPM)
    ) %>%
    group_by(tissue, sample_id, gene_symbol) %>%
    summarise(
      TPM = mean(TPM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log2_TPM = log2(TPM + 1)
    )
  
  gct_long
}

run_one_wilcox <- function(signature_name, comparison_id, group_1, group_2, score_data) {
  
  df_sig <- score_data %>%
    filter(signature == signature_name)
  
  x <- df_sig %>%
    filter(tissue == group_1) %>%
    pull(signature_score) %>%
    na.omit()
  
  y <- df_sig %>%
    filter(tissue == group_2) %>%
    pull(signature_score) %>%
    na.omit()
  
  if (length(x) < 2 || length(y) < 2) {
    return(
      tibble(
        signature = signature_name,
        comparison_id = comparison_id,
        group_1 = group_1,
        group_2 = group_2,
        n_group_1 = length(x),
        n_group_2 = length(y),
        median_group_1 = NA_real_,
        median_group_2 = NA_real_,
        median_difference = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- wilcox.test(
    x = x,
    y = y,
    alternative = "two.sided",
    exact = FALSE
  )
  
  tibble(
    signature = signature_name,
    comparison_id = comparison_id,
    group_1 = group_1,
    group_2 = group_2,
    n_group_1 = length(x),
    n_group_2 = length(y),
    median_group_1 = median(x, na.rm = TRUE),
    median_group_2 = median(y, na.rm = TRUE),
    median_difference = median_group_1 - median_group_2,
    p_value = test$p.value
  )
}

###############################################################################
# 2) RESOLVE AND CHECK INPUT FILES
###############################################################################

gtex_files <- gtex_files %>%
  mutate(
    file = map_chr(file_name, resolve_gtex_file, expr_dir = EXPR_DIR)
  )

check_files_exist(c(mitoall_file, mitoonly_file, gtex_files$file))

message("\nInput files successfully located.\n")

###############################################################################
# 3) IMPORT MITOALL AND MITOONLY SIGNATURES
###############################################################################

mitoall_genes <- import_signature(
  file = mitoall_file,
  signature_name = "MitoAll"
)

mitoonly_genes <- import_signature(
  file = mitoonly_file,
  signature_name = "MitoOnly"
)

signature_genes <- bind_rows(
  mitoall_genes,
  mitoonly_genes
) %>%
  distinct(signature, gene_symbol)

target_genes <- signature_genes %>%
  pull(gene_symbol) %>%
  unique()

message("MitoAll genes: ", n_distinct(mitoall_genes$gene_symbol))
message("MitoOnly genes: ", n_distinct(mitoonly_genes$gene_symbol))
message("Total unique signature genes: ", length(target_genes))

###############################################################################
# 4) IMPORT GTEX TPM FILES AND TRANSFORM EXPRESSION
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

sample_summary <- gtex_expr_long %>%
  distinct(tissue, sample_id) %>%
  count(tissue, name = "n_samples") %>%
  arrange(tissue)

gene_summary_by_tissue <- gtex_expr_long %>%
  distinct(tissue, gene_symbol) %>%
  count(tissue, name = "n_signature_genes_detected") %>%
  arrange(tissue)

signature_coverage_import <- signature_genes %>%
  left_join(
    gtex_expr_long %>%
      distinct(gene_symbol) %>%
      mutate(detected_in_gtex = TRUE),
    by = "gene_symbol"
  ) %>%
  mutate(
    detected_in_gtex = if_else(is.na(detected_in_gtex), FALSE, detected_in_gtex)
  ) %>%
  group_by(signature) %>%
  summarise(
    genes_in_signature = n_distinct(gene_symbol),
    genes_detected_in_gtex = n_distinct(gene_symbol[detected_in_gtex]),
    coverage = genes_detected_in_gtex / genes_in_signature,
    .groups = "drop"
  )

###############################################################################
# 6) COMPUTE GENE-WISE Z-SCORES
###############################################################################
# Z-scores are calculated gene by gene across all selected GTEx samples.

gtex_expr_z <- gtex_expr_long %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_mean_log2TPM = mean(log2_TPM, na.rm = TRUE),
    gene_sd_log2TPM = sd(log2_TPM, na.rm = TRUE),
    z_log2_TPM = if_else(
      is.na(gene_sd_log2TPM) | gene_sd_log2TPM == 0,
      NA_real_,
      (log2_TPM - gene_mean_log2TPM) / gene_sd_log2TPM
    )
  ) %>%
  ungroup()

###############################################################################
# 7) CALCULATE SAMPLE-LEVEL SIGNATURE SCORES
###############################################################################

signature_scores <- gtex_expr_z %>%
  inner_join(signature_genes, by = "gene_symbol") %>%
  group_by(tissue, sample_id, signature) %>%
  summarise(
    n_genes_used = sum(!is.na(z_log2_TPM)),
    signature_score = mean(z_log2_TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    signature_score = if_else(n_genes_used == 0, NA_real_, signature_score)
  )

###############################################################################
# 8) CHECK SIGNATURE COVERAGE FOR SCORING
###############################################################################

tissues <- sort(unique(gtex_expr_long$tissue))

detected_genes_by_tissue <- gtex_expr_z %>%
  filter(!is.na(z_log2_TPM)) %>%
  distinct(tissue, gene_symbol) %>%
  mutate(detected = TRUE)

coverage_by_tissue <- tidyr::expand_grid(
  tissue = tissues,
  signature_genes
) %>%
  left_join(
    detected_genes_by_tissue,
    by = c("tissue", "gene_symbol")
  ) %>%
  mutate(
    detected = if_else(is.na(detected), FALSE, detected)
  ) %>%
  group_by(tissue, signature) %>%
  summarise(
    genes_in_signature = n_distinct(gene_symbol),
    genes_detected_for_scoring = n_distinct(gene_symbol[detected]),
    coverage = genes_detected_for_scoring / genes_in_signature,
    .groups = "drop"
  ) %>%
  arrange(signature, tissue)

coverage_overall <- signature_genes %>%
  left_join(
    gtex_expr_z %>%
      filter(!is.na(z_log2_TPM)) %>%
      distinct(gene_symbol) %>%
      mutate(detected = TRUE),
    by = "gene_symbol"
  ) %>%
  mutate(
    detected = if_else(is.na(detected), FALSE, detected)
  ) %>%
  group_by(signature) %>%
  summarise(
    genes_in_signature = n_distinct(gene_symbol),
    genes_detected_for_scoring = n_distinct(gene_symbol[detected]),
    coverage = genes_detected_for_scoring / genes_in_signature,
    .groups = "drop"
  )

###############################################################################
# 9) SUMMARIZE SIGNATURE SCORES BY TISSUE
###############################################################################

score_summary_by_tissue <- signature_scores %>%
  group_by(signature, tissue) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    median_score = median(signature_score, na.rm = TRUE),
    mean_score = mean(signature_score, na.rm = TRUE),
    q1 = quantile(signature_score, 0.25, na.rm = TRUE),
    q3 = quantile(signature_score, 0.75, na.rm = TRUE),
    iqr = IQR(signature_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(signature, desc(median_score))

###############################################################################
# 10) WILCOXON RANK-SUM TESTS
###############################################################################

planned_comparisons <- tibble::tribble(
  ~comparison_id, ~group_1, ~group_2,
  "Heart_Left_Ventricle_vs_Whole_Blood",   "Heart - Left Ventricle",   "Whole Blood",
  "Heart_Atrial_Appendage_vs_Whole_Blood", "Heart - Atrial Appendage", "Whole Blood",
  "Muscle_Skeletal_vs_Whole_Blood",        "Muscle - Skeletal",        "Whole Blood"
)

wilcox_results <- tidyr::expand_grid(
  signature = sort(unique(signature_scores$signature)),
  planned_comparisons
) %>%
  pmap_dfr(
    function(signature, comparison_id, group_1, group_2) {
      run_one_wilcox(
        signature_name = signature,
        comparison_id = comparison_id,
        group_1 = group_1,
        group_2 = group_2,
        score_data = signature_scores
      )
    }
  ) %>%
  group_by(signature) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH"),
    direction = case_when(
      median_difference > 0 ~ "Higher in group_1",
      median_difference < 0 ~ "Lower in group_1",
      median_difference == 0 ~ "No median difference",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  arrange(signature, p_adj_fdr)

###############################################################################
# 11) EXPORT RESULTS
###############################################################################

saveRDS(
  signature_genes,
  file.path(RDS_DIR, "signature_genes_clean.rds")
)

saveRDS(
  gtex_expr_long,
  file.path(RDS_DIR, "gtex_mito_signature_genes_log2TPM_long.rds")
)

saveRDS(
  gtex_expr_z,
  file.path(RDS_DIR, "gtex_expr_gene_wise_zscores.rds")
)

saveRDS(
  signature_scores,
  file.path(RDS_DIR, "gtex_mitoall_mitoonly_sample_scores.rds")
)

write_csv(
  sample_summary,
  file.path(TABLE_DIR, "sample_summary_by_tissue.csv")
)

write_csv(
  gene_summary_by_tissue,
  file.path(TABLE_DIR, "detected_signature_genes_by_tissue.csv")
)

write_csv(
  signature_coverage_import,
  file.path(TABLE_DIR, "signature_coverage_summary_import.csv")
)

write_csv(
  coverage_by_tissue,
  file.path(TABLE_DIR, "signature_coverage_by_tissue.csv")
)

write_csv(
  coverage_overall,
  file.path(TABLE_DIR, "signature_coverage_overall.csv")
)

write_csv(
  signature_scores,
  file.path(TABLE_DIR, "gtex_mitoall_mitoonly_sample_scores.csv")
)

write_csv(
  score_summary_by_tissue,
  file.path(TABLE_DIR, "signature_score_summary_by_tissue.csv")
)

write_csv(
  wilcox_results,
  file.path(TABLE_DIR, "wilcoxon_planned_comparisons_fdr.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(INFO_DIR, "run_pipeline_sessionInfo.txt")
)

###############################################################################
# 12) PRINT SUMMARY
###############################################################################

message("\n============================================================")
message("PIPELINE COMPLETED SUCCESSFULLY")
message("============================================================")

message("\nSample summary:")
print(sample_summary)

message("\nOverall signature coverage:")
print(coverage_overall)

message("\nScore summary by tissue:")
print(score_summary_by_tissue)

message("\nWilcoxon planned comparisons:")
print(wilcox_results)

message("\nResults saved in:")
message("Tables: ", TABLE_DIR)
message("RDS:    ", RDS_DIR)
message("Session info: ", INFO_DIR)

###############################################################################
# FIGURE 1 — REPRESENTATIVE MITOCHONDRIAL SIGNATURE SCORING RESULT
#
# Purpose:
#   Generate a publication-ready representative result figure showing
#   MitoAll and MitoOnly sample-level scores across selected GTEx tissues.
#
# Input:
#   Results/Tables/gtex_mitoall_mitoonly_sample_scores.csv
#
# Output:
#   Results/Figures/Figure_1.png
#   Results/Figures/Figure_1.pdf
#   Results/Figures/Figure_1.tiff
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
})

###############################################################################
# 0) SETTINGS
###############################################################################

BASE_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

TABLE_DIR <- file.path(BASE_DIR, "Results", "Tables")
FIG_DIR   <- file.path(BASE_DIR, "Results", "Figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(TABLE_DIR, "gtex_mitoall_mitoonly_sample_scores.csv")

if (!file.exists(score_file)) {
  stop(
    "Score file not found:\n",
    score_file,
    "\nRun the main pipeline before generating Figure 1."
  )
}

###############################################################################
# 1) IMPORT SCORES
###############################################################################

scores <- read_csv(score_file, show_col_types = FALSE)

required_cols <- c("tissue", "sample_id", "signature", "signature_score", "n_genes_used")

if (!all(required_cols %in% colnames(scores))) {
  stop(
    "The score table does not contain the required columns:\n",
    paste(required_cols, collapse = ", ")
  )
}

###############################################################################
# 2) HARMONIZE LABELS AND ORDER
###############################################################################

tissue_order <- c(
  "Whole Blood",
  "Heart - Atrial Appendage",
  "Heart - Left Ventricle",
  "Muscle - Skeletal"
)

tissue_labels <- c(
  "Whole Blood" = "Whole\nBlood",
  "Heart - Atrial Appendage" = "Heart\nAtrial appendage",
  "Heart - Left Ventricle" = "Heart\nLeft ventricle",
  "Muscle - Skeletal" = "Skeletal\nmuscle"
)

signature_order <- c("MitoAll", "MitoOnly")

signature_labels <- c(
  "MitoAll" = "MitoAll",
  "MitoOnly" = "MitoOnly"
)

scores <- scores %>%
  filter(
    tissue %in% tissue_order,
    signature %in% signature_order,
    !is.na(signature_score)
  ) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    signature = factor(signature, levels = signature_order)
  )

if (nrow(scores) == 0) {
  stop("No valid MitoAll/MitoOnly scores were found after filtering.")
}

###############################################################################
# 3) COLORS
###############################################################################

tissue_palette <- c(
  "Whole Blood" = "#B2182B",
  "Heart - Atrial Appendage" = "#2166AC",
  "Heart - Left Ventricle" = "#1B9E77",
  "Muscle - Skeletal" = "#D95F02"
)

###############################################################################
# 4) SUMMARY TABLE FOR FIGURE
###############################################################################

figure_summary <- scores %>%
  group_by(signature, tissue) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    median_score = median(signature_score, na.rm = TRUE),
    mean_score = mean(signature_score, na.rm = TRUE),
    q1 = quantile(signature_score, 0.25, na.rm = TRUE),
    q3 = quantile(signature_score, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  figure_summary,
  file.path(TABLE_DIR, "Figure_1_score_summary_by_tissue.csv")
)

###############################################################################
# 5) THEME
###############################################################################

theme_figure <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(
      size = 10,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.0, "lines"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

###############################################################################
# 6) BUILD FIGURE
###############################################################################

figure_1 <- ggplot(
  scores,
  aes(x = tissue, y = signature_score, fill = tissue)
) +
  geom_violin(
    width = 0.90,
    alpha = 0.32,
    color = NA,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.92,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = tissue),
    width = 0.14,
    size = 0.55,
    alpha = 0.16,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.3,
    fill = "white",
    color = "black",
    stroke = 0.35
  ) +
  facet_wrap(
    ~ signature,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(signature_labels)
  ) +
  scale_fill_manual(values = tissue_palette) +
  scale_color_manual(values = tissue_palette) +
  scale_x_discrete(labels = tissue_labels) +
  labs(
    x = NULL,
    y = "Sample-level mitochondrial signature score\n(mean gene-wise z-score)"
  ) +
  theme_figure

###############################################################################
# 7) EXPORT FIGURE
###############################################################################

ggsave(
  filename = file.path(FIG_DIR, "Figure_1.pdf"),
  plot = figure_1,
  width = 10,
  height = 5.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(FIG_DIR, "Figure_1.png"),
  plot = figure_1,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(FIG_DIR, "Figure_1.tiff"),
  plot = figure_1,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################################
# 8) PRINT
###############################################################################

print(figure_1)

message("\nFigure 1 generated successfully.")
message("Saved in: ", FIG_DIR)
