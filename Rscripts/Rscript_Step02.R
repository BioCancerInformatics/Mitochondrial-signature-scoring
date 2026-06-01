###############################################################################
# GTEx MITOCHONDRIAL SIGNATURE ANALYSIS
# Step 2 — Z-scores, MitoAll/MitoOnly scores, coverage, Wilcoxon tests
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
})

###############################################################################
# 0) USER SETTINGS
###############################################################################

BASE_DIR <- "C:/Users/CORDEIH/Desktop/Whole_Project/Manuscript/Analyses"
OUT_DIR  <- file.path(BASE_DIR, "output")

# Objects generated in Step 1
signature_file <- file.path(OUT_DIR, "signature_genes_clean.rds")
expr_file      <- file.path(OUT_DIR, "gtex_mito_signature_genes_log2TPM_long.rds")

# Output directory for Step 2
STEP2_DIR <- file.path(OUT_DIR, "step2_signature_scores")
dir.create(STEP2_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1) LOAD CLEANED OBJECTS FROM STEP 1
###############################################################################

signature_genes <- readRDS(signature_file)
gtex_expr_long  <- readRDS(expr_file)

# Expected columns:
# signature_genes: signature, gene_symbol
# gtex_expr_long: tissue, sample_id, gene_symbol, TPM, log2_TPM

required_signature_cols <- c("signature", "gene_symbol")
required_expr_cols <- c("tissue", "sample_id", "gene_symbol", "TPM", "log2_TPM")

if (!all(required_signature_cols %in% colnames(signature_genes))) {
  stop("signature_genes does not contain the expected columns: signature, gene_symbol")
}

if (!all(required_expr_cols %in% colnames(gtex_expr_long))) {
  stop("gtex_expr_long does not contain the expected columns: tissue, sample_id, gene_symbol, TPM, log2_TPM")
}

###############################################################################
# 2) ORGANIZE SIGNATURE GENE UNIVERSE
###############################################################################

signature_gene_universe <- signature_genes %>%
  mutate(
    signature = as.character(signature),
    gene_symbol = as.character(gene_symbol)
  ) %>%
  filter(!is.na(signature), !is.na(gene_symbol)) %>%
  distinct(signature, gene_symbol)

message("Signatures loaded:")
print(signature_gene_universe %>% count(signature, name = "n_genes"))

###############################################################################
# 3) COMPUTE GENE-WISE Z-SCORES
###############################################################################
# Important:
# Z-scores are computed gene by gene across all selected GTEx samples.
# This makes each gene comparable before averaging genes into a signature score.

gtex_expr_z <- gtex_expr_long %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_mean_log2TPM = mean(log2_TPM, na.rm = TRUE),
    gene_sd_log2TPM   = sd(log2_TPM, na.rm = TRUE),
    z_log2_TPM = if_else(
      is.na(gene_sd_log2TPM) | gene_sd_log2TPM == 0,
      NA_real_,
      (log2_TPM - gene_mean_log2TPM) / gene_sd_log2TPM
    )
  ) %>%
  ungroup()

###############################################################################
# 4) CALCULATE SAMPLE-LEVEL MITOALL AND MITOONLY SCORES
###############################################################################
# Each score is the mean of the gene-wise z-scores for the genes belonging to
# that signature in that sample.

signature_scores <- gtex_expr_z %>%
  inner_join(signature_gene_universe, by = "gene_symbol") %>%
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
# 5) CHECK SIGNATURE COVERAGE
###############################################################################
# This repeats the coverage check, but now based specifically on genes that
# actually entered the scoring step.

tissues <- sort(unique(gtex_expr_long$tissue))

detected_genes_by_tissue <- gtex_expr_z %>%
  filter(!is.na(z_log2_TPM)) %>%
  distinct(tissue, gene_symbol) %>%
  mutate(detected = TRUE)

coverage_by_tissue <- tidyr::expand_grid(
  tissue = tissues,
  signature_gene_universe
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

coverage_overall <- signature_gene_universe %>%
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
# 6) SUMMARIZE SIGNATURE SCORES BY TISSUE
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
# 7) WILCOXON RANK-SUM TESTS
###############################################################################
# Biological rationale:
# group_1 should be the tissue expected to have higher mitochondrial abundance.
# Positive median_difference means group_1 has higher signature score than group_2.

planned_comparisons <- tibble::tribble(
  ~comparison_id, ~group_1, ~group_2,
  # "Heart_Left_Ventricle_vs_Spleen",       "Heart - Left Ventricle",   "Spleen",
  # "Heart_Atrial_Appendage_vs_Spleen",     "Heart - Atrial Appendage", "Spleen",
  # "Muscle_Skeletal_vs_Spleen",            "Muscle - Skeletal",        "Spleen",
  
  "Heart_Left_Ventricle_vs_Whole_Blood",   "Heart - Left Ventricle",   "Whole Blood",
  "Heart_Atrial_Appendage_vs_Whole_Blood", "Heart - Atrial Appendage", "Whole Blood",
  "Muscle_Skeletal_vs_Whole_Blood",        "Muscle - Skeletal",        "Whole Blood"
)

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
  
  if (length(x) < 2 | length(y) < 2) {
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
# 8) EXPORT RESULTS
###############################################################################

saveRDS(
  gtex_expr_z,
  file.path(STEP2_DIR, "gtex_expr_gene_wise_zscores.rds")
)

saveRDS(
  signature_scores,
  file.path(STEP2_DIR, "gtex_mitoall_mitoonly_mtrna_sample_scores.rds")
)

write_csv(
  signature_scores,
  file.path(STEP2_DIR, "gtex_mitoall_mitoonly_mtrna_sample_scores.csv")
)

write_csv(
  coverage_by_tissue,
  file.path(STEP2_DIR, "signature_coverage_by_tissue.csv")
)

write_csv(
  coverage_overall,
  file.path(STEP2_DIR, "signature_coverage_overall.csv")
)

write_csv(
  score_summary_by_tissue,
  file.path(STEP2_DIR, "signature_score_summary_by_tissue.csv")
)

write_csv(
  wilcox_results,
  file.path(STEP2_DIR, "wilcoxon_planned_comparisons_fdr.csv")
)

###############################################################################
# 9) PRINT SHORT SUMMARY
###############################################################################

message("\nOverall signature coverage:")
print(coverage_overall)

message("\nScore summary by tissue:")
print(score_summary_by_tissue)

message("\nWilcoxon planned comparisons:")
print(wilcox_results)

message("\nDone.")
message("Results saved in: ", STEP2_DIR)


###############################################################################
# 7.1) mtRNA VALIDATION ANALYSIS
###############################################################################
# Aim:
# Compare nuclear mitochondrial signatures with mitochondrial-encoded RNA signal.
#
# Tests:
# 1) MitoAll vs mtRNA_13PCG score
# 2) MitoOnly vs mtRNA_13PCG score
# 3) Tissue-level agreement among MitoAll, MitoOnly, and mtRNA_13PCG

MTRNA_DIR <- file.path(STEP2_DIR, "mtrna_validation")
dir.create(MTRNA_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 7.1.1) Prepare wide score matrix
###############################################################################

signature_scores_wide <- signature_scores %>%
  select(tissue, sample_id, signature, signature_score) %>%
  distinct() %>%
  pivot_wider(
    names_from = signature,
    values_from = signature_score
  )

required_score_cols <- c("MitoAll", "MitoOnly", "mtRNA_13PCG")

missing_score_cols <- setdiff(required_score_cols, colnames(signature_scores_wide))

if (length(missing_score_cols) > 0) {
  stop(
    "The following signature score columns are missing: ",
    paste(missing_score_cols, collapse = ", "),
    "\nCheck whether the 13 mtRNA genes were imported in Step 1."
  )
}

###############################################################################
# 7.1.2) Sample-level Spearman correlations
###############################################################################

run_spearman <- function(df, x_col, y_col, comparison_name) {
  
  df_test <- df %>%
    select(all_of(c(x_col, y_col))) %>%
    filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]])
    )
  
  test <- cor.test(
    df_test[[x_col]],
    df_test[[y_col]],
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    comparison = comparison_name,
    x = x_col,
    y = y_col,
    n_samples = nrow(df_test),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value
  )
}

sample_level_correlations <- bind_rows(
  run_spearman(
    df = signature_scores_wide,
    x_col = "MitoAll",
    y_col = "mtRNA_13PCG",
    comparison_name = "MitoAll_vs_mtRNA_13PCG"
  ),
  run_spearman(
    df = signature_scores_wide,
    x_col = "MitoOnly",
    y_col = "mtRNA_13PCG",
    comparison_name = "MitoOnly_vs_mtRNA_13PCG"
  ),
  run_spearman(
    df = signature_scores_wide,
    x_col = "MitoAll",
    y_col = "MitoOnly",
    comparison_name = "MitoAll_vs_MitoOnly"
  )
) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH")
  )

###############################################################################
# 7.1.3) Tissue-level median scores
###############################################################################

tissue_level_scores <- signature_scores %>%
  group_by(tissue, signature) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    median_score = median(signature_score, na.rm = TRUE),
    mean_score = mean(signature_score, na.rm = TRUE),
    q1 = quantile(signature_score, 0.25, na.rm = TRUE),
    q3 = quantile(signature_score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(signature, desc(median_score))

tissue_level_scores_wide <- tissue_level_scores %>%
  select(tissue, signature, median_score) %>%
  pivot_wider(
    names_from = signature,
    values_from = median_score
  )

###############################################################################
# 7.1.4) Tissue-level agreement
###############################################################################
# Caution:
# With only four tissues, this is descriptive and should not be overinterpreted.
# It is useful as a directional agreement analysis.

tissue_level_correlations <- bind_rows(
  run_spearman(
    df = tissue_level_scores_wide,
    x_col = "MitoAll",
    y_col = "mtRNA_13PCG",
    comparison_name = "Tissue_level_MitoAll_vs_mtRNA_13PCG"
  ),
  run_spearman(
    df = tissue_level_scores_wide,
    x_col = "MitoOnly",
    y_col = "mtRNA_13PCG",
    comparison_name = "Tissue_level_MitoOnly_vs_mtRNA_13PCG"
  ),
  run_spearman(
    df = tissue_level_scores_wide,
    x_col = "MitoAll",
    y_col = "MitoOnly",
    comparison_name = "Tissue_level_MitoAll_vs_MitoOnly"
  )
) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH")
  )

###############################################################################
# 7.1.5) Export mtRNA validation results
###############################################################################

write_csv(
  signature_scores_wide,
  file.path(MTRNA_DIR, "sample_level_signature_scores_wide.csv")
)

write_csv(
  sample_level_correlations,
  file.path(MTRNA_DIR, "sample_level_spearman_correlations_with_mtrna.csv")
)

write_csv(
  tissue_level_scores,
  file.path(MTRNA_DIR, "tissue_level_signature_score_summary.csv")
)

write_csv(
  tissue_level_scores_wide,
  file.path(MTRNA_DIR, "tissue_level_median_scores_wide.csv")
)

write_csv(
  tissue_level_correlations,
  file.path(MTRNA_DIR, "tissue_level_spearman_correlations_with_mtrna.csv")
)

message("\nmtRNA validation analysis completed.")
message("mtRNA validation results saved in: ", MTRNA_DIR)

message("\nSample-level correlations:")
print(sample_level_correlations)

message("\nTissue-level median scores:")
print(tissue_level_scores_wide)

message("\nTissue-level correlations:")
print(tissue_level_correlations)

# ###############################################################################
# # DIAGNOSTIC ANALYSIS — CHECK WHETHER THE HIGH COMPARATOR IS A SCORING ARTIFACT
# ###############################################################################
# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(readr)
# })
# 
# DIAG_DIR <- file.path(OUT_DIR, "diagnostics_signature_scoring")
# dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
# 
# ###############################################################################
# # 1) RAW MEAN LOG2TPM SIGNATURE SCORE
# ###############################################################################
# # This score does not use z-scores.
# # If the same tissue remains highest, the pattern is not caused only by z-scoring.
# 
# raw_signature_scores <- gtex_expr_long %>%
#   inner_join(signature_gene_universe, by = "gene_symbol") %>%
#   group_by(tissue, sample_id, signature) %>%
#   summarise(
#     n_genes_used = sum(!is.na(log2_TPM)),
#     raw_mean_log2TPM_score = mean(log2_TPM, na.rm = TRUE),
#     raw_median_log2TPM_score = median(log2_TPM, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# raw_score_summary <- raw_signature_scores %>%
#   group_by(signature, tissue) %>%
#   summarise(
#     n_samples = n_distinct(sample_id),
#     median_raw_mean_score = median(raw_mean_log2TPM_score, na.rm = TRUE),
#     mean_raw_mean_score = mean(raw_mean_log2TPM_score, na.rm = TRUE),
#     q1 = quantile(raw_mean_log2TPM_score, 0.25, na.rm = TRUE),
#     q3 = quantile(raw_mean_log2TPM_score, 0.75, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(signature, desc(median_raw_mean_score))
# 
# write_csv(
#   raw_signature_scores,
#   file.path(DIAG_DIR, "raw_log2TPM_signature_scores_per_sample.csv")
# )
# 
# write_csv(
#   raw_score_summary,
#   file.path(DIAG_DIR, "raw_log2TPM_signature_score_summary_by_tissue.csv")
# )
# 
# ###############################################################################
# # 2) TISSUE-BALANCED GENE-WISE Z-SCORE
# ###############################################################################
# # This avoids giving more weight to tissues with larger sample sizes.
# # First calculate the mean expression of each gene in each tissue.
# # Then calculate gene-level mean and SD across tissue means, not across all samples.
# 
# gene_tissue_reference <- gtex_expr_long %>%
#   group_by(gene_symbol, tissue) %>%
#   summarise(
#     tissue_mean_log2TPM = mean(log2_TPM, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   group_by(gene_symbol) %>%
#   summarise(
#     ref_mean_tissue_balanced = mean(tissue_mean_log2TPM, na.rm = TRUE),
#     ref_sd_tissue_balanced = sd(tissue_mean_log2TPM, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# gtex_expr_tissue_balanced_z <- gtex_expr_long %>%
#   left_join(gene_tissue_reference, by = "gene_symbol") %>%
#   mutate(
#     z_log2_TPM_tissue_balanced = if_else(
#       is.na(ref_sd_tissue_balanced) | ref_sd_tissue_balanced == 0,
#       NA_real_,
#       (log2_TPM - ref_mean_tissue_balanced) / ref_sd_tissue_balanced
#     )
#   )
# 
# tissue_balanced_signature_scores <- gtex_expr_tissue_balanced_z %>%
#   inner_join(signature_gene_universe, by = "gene_symbol") %>%
#   group_by(tissue, sample_id, signature) %>%
#   summarise(
#     n_genes_used = sum(!is.na(z_log2_TPM_tissue_balanced)),
#     tissue_balanced_signature_score = mean(z_log2_TPM_tissue_balanced, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# tissue_balanced_score_summary <- tissue_balanced_signature_scores %>%
#   group_by(signature, tissue) %>%
#   summarise(
#     n_samples = n_distinct(sample_id),
#     median_tissue_balanced_score = median(tissue_balanced_signature_score, na.rm = TRUE),
#     mean_tissue_balanced_score = mean(tissue_balanced_signature_score, na.rm = TRUE),
#     q1 = quantile(tissue_balanced_signature_score, 0.25, na.rm = TRUE),
#     q3 = quantile(tissue_balanced_signature_score, 0.75, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(signature, desc(median_tissue_balanced_score))
# 
# write_csv(
#   tissue_balanced_signature_scores,
#   file.path(DIAG_DIR, "tissue_balanced_z_signature_scores_per_sample.csv")
# )
# 
# write_csv(
#   tissue_balanced_score_summary,
#   file.path(DIAG_DIR, "tissue_balanced_z_signature_score_summary_by_tissue.csv")
# )
# 
# ###############################################################################
# # 3) EXPRESSION-LEVEL COVERAGE
# ###############################################################################
# # This checks whether genes are merely present in the matrix or actually expressed.
# 
# expression_coverage <- gtex_expr_long %>%
#   inner_join(signature_gene_universe, by = "gene_symbol") %>%
#   group_by(tissue, signature, gene_symbol) %>%
#   summarise(
#     median_TPM = median(TPM, na.rm = TRUE),
#     mean_TPM = mean(TPM, na.rm = TRUE),
#     expressed_TPM_gt_0 = median_TPM > 0,
#     expressed_TPM_gt_1 = median_TPM > 1,
#     .groups = "drop"
#   ) %>%
#   group_by(tissue, signature) %>%
#   summarise(
#     genes_in_signature = n_distinct(gene_symbol),
#     genes_with_median_TPM_gt_0 = sum(expressed_TPM_gt_0),
#     genes_with_median_TPM_gt_1 = sum(expressed_TPM_gt_1),
#     fraction_median_TPM_gt_0 = genes_with_median_TPM_gt_0 / genes_in_signature,
#     fraction_median_TPM_gt_1 = genes_with_median_TPM_gt_1 / genes_in_signature,
#     .groups = "drop"
#   )
# 
# write_csv(
#   expression_coverage,
#   file.path(DIAG_DIR, "expression_level_signature_coverage_by_tissue.csv")
# )
# 
# ###############################################################################
# # 4) GENE-LEVEL DRIVERS OF HIGH SCORE
# ###############################################################################
# # This identifies which mitochondrial genes drive the difference between tissues.
# # Example: compare Spleen or Lung against Muscle.
# 
# driver_analysis <- gtex_expr_long %>%
#   inner_join(signature_gene_universe, by = "gene_symbol") %>%
#   group_by(signature, tissue, gene_symbol) %>%
#   summarise(
#     median_log2TPM = median(log2_TPM, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = tissue,
#     values_from = median_log2TPM
#   )
# 
# write_csv(
#   driver_analysis,
#   file.path(DIAG_DIR, "gene_level_median_log2TPM_by_tissue_for_driver_analysis.csv")
# )
# 
# ###############################################################################
# # 5) PRINT DIAGNOSTIC SUMMARIES
# ###############################################################################
# 
# message("\nRaw log2TPM score summary:")
# print(raw_score_summary)
# 
# message("\nTissue-balanced z-score summary:")
# print(tissue_balanced_score_summary)
# 
# message("\nExpression-level coverage:")
# print(expression_coverage)
# 
# message("\nDiagnostic analysis completed.")
# message("Diagnostic results saved in: ", DIAG_DIR)
