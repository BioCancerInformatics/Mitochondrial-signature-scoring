###############################################################################
# GTEx MITOCHONDRIAL SIGNATURE ANALYSIS
# Figure 1 — Distribution of MitoAll, MitoOnly, and mtRNA_13PCG scores
#            + tissue-level median agreement
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
  library(svglite)
})

###############################################################################
# 0) USER SETTINGS
###############################################################################

BASE_DIR  <- "C:/Users/CORDEIH/Desktop/Whole_Project/Manuscript/Analyses"
STEP2_DIR <- file.path(BASE_DIR, "output", "step2_signature_scores")
FIG_DIR   <- file.path(BASE_DIR, "output", "figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(STEP2_DIR, "gtex_mitoall_mitoonly_mtrna_sample_scores.csv")

###############################################################################
# 1) IMPORT DATA
###############################################################################

scores <- read_csv(score_file, show_col_types = FALSE)

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

signature_order <- c(
  "MitoAll",
  "MitoOnly",
  "mtRNA_13PCG"
)

signature_labels <- c(
  "MitoAll" = "MitoAll",
  "MitoOnly" = "MitoOnly",
  "mtRNA_13PCG" = "mtRNA\n13 PCGs"
)

expected_cols <- c("tissue", "sample_id", "signature", "signature_score")

if (!all(expected_cols %in% colnames(scores))) {
  stop("The score file does not contain the expected columns.")
}

missing_signatures <- setdiff(signature_order, unique(scores$signature))

if (length(missing_signatures) > 0) {
  stop(
    "The following signatures are missing from the score file: ",
    paste(missing_signatures, collapse = ", "),
    "\nPlease rerun Step 1 and Step 2 after adding mtRNA_13PCG."
  )
}

scores <- scores %>%
  filter(signature %in% signature_order) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    signature = factor(signature, levels = signature_order)
  ) %>%
  filter(
    !is.na(tissue),
    !is.na(signature),
    !is.na(signature_score)
  )

###############################################################################
# 3) COLORS
###############################################################################

tissue_palette <- c(
  "Whole Blood" = "#B2182B",
  "Heart - Atrial Appendage" = "#2166AC",
  "Heart - Left Ventricle" = "#1B9E77",
  "Muscle - Skeletal" = "#D95F02"
)

signature_palette <- c(
  "MitoAll" = "#2166AC",
  "MitoOnly" = "#B2182B",
  "mtRNA_13PCG" = "#4D4D4D"
)

###############################################################################
# 4) COMMON THEME
###############################################################################

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

###############################################################################
# 5) FUNCTION TO BUILD DISTRIBUTION PANELS
###############################################################################

make_distribution_plot <- function(signature_name, plot_title) {
  
  scores %>%
    filter(signature == signature_name) %>%
    ggplot(aes(x = tissue, y = signature_score, fill = tissue)) +
    geom_violin(
      width = 0.95,
      alpha = 0.35,
      color = NA,
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.85,
      color = "black",
      linewidth = 0.35
    ) +
    geom_jitter(
      aes(color = tissue),
      width = 0.16,
      size = 0.65,
      alpha = 0.16,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.5,
      fill = "white",
      color = "black"
    ) +
    scale_fill_manual(values = tissue_palette, labels = tissue_labels) +
    scale_color_manual(values = tissue_palette, labels = tissue_labels) +
    scale_x_discrete(labels = tissue_labels) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Sample-level score"
    ) +
    theme_pub +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        size = 9,
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      )
    )
}

###############################################################################
# 6) BUILD TOP DISTRIBUTION PANELS
###############################################################################

plot_A <- make_distribution_plot(
  signature_name = "MitoAll",
  plot_title = "MitoAll"
)

plot_B <- make_distribution_plot(
  signature_name = "MitoOnly",
  plot_title = "MitoOnly"
)

plot_C <- make_distribution_plot(
  signature_name = "mtRNA_13PCG",
  plot_title = "mtRNA 13 PCGs"
)

###############################################################################
# 7) TISSUE-LEVEL MEDIAN AGREEMENT
###############################################################################

tissue_median_long <- scores %>%
  group_by(tissue, signature) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    median_score = median(signature_score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  tissue_median_long,
  file.path(FIG_DIR, "Figure1_tissue_level_median_scores.csv")
)

plot_D <- ggplot(
  tissue_median_long,
  aes(x = tissue, y = median_score, group = signature, color = signature)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2) +
  scale_color_manual(values = signature_palette, labels = signature_labels) +
  scale_x_discrete(labels = tissue_labels) +
  labs(
    title = "Tissue-level median agreement",
    x = NULL,
    y = "Median sample-level score",
    color = "Signature"
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      size = 10,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    )
  )

###############################################################################
# 8) COMBINE FIGURE 1
###############################################################################

figure_1 <-
  (plot_A | plot_B | plot_C) /
  plot_D +
  plot_layout(
    heights = c(1.2, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

###############################################################################
# 9) SAVE FIGURE 1
###############################################################################

ggsave(
  filename = file.path(FIG_DIR, "Figure1_MitoSignatures_Distribution_MedianAgreement.pdf"),
  plot = figure_1,
  width = 16,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(FIG_DIR, "Figure1_MitoSignatures_Distribution_MedianAgreement.tiff"),
  plot = figure_1,
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(FIG_DIR, "Figure1_MitoSignatures_Distribution_MedianAgreement.png"),
  plot = figure_1,
  width = 16,
  height = 10,
  units = "in",
  dpi = 600
)

print(figure_1)


###############################################################################
# GTEx MITOCHONDRIAL SIGNATURE ANALYSIS
# Figure 2 — MitoAll/MitoOnly versus mtRNA_13PCG correlations
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

# ###############################################################################
# # 0) USER SETTINGS
# ###############################################################################
# 
# BASE_DIR  <- "C:/Users/CORDEIH/Desktop/Whole_Project/Manuscript/Analyses"
# STEP2_DIR <- file.path(BASE_DIR, "output", "step2_signature_scores")
# FIG_DIR   <- file.path(BASE_DIR, "output", "figures")
# 
# dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
# 
# score_file <- file.path(STEP2_DIR, "gtex_mitoall_mitoonly_sample_scores.csv")

###############################################################################
# 1) IMPORT DATA
###############################################################################

scores <- read_csv(score_file, show_col_types = FALSE)

###############################################################################
# 2) HELPER FUNCTIONS
###############################################################################

format_p_value <- function(p, prefix = "p") {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ paste0(prefix, " < 0.001"),
    TRUE ~ paste0(prefix, " = ", signif(p, 3))
  )
}

run_spearman <- function(df, x_col, y_col, comparison_name) {
  
  df_test <- df %>%
    select(all_of(c(x_col, y_col))) %>%
    filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]])
    )
  
  if (nrow(df_test) < 3) {
    return(
      tibble(
        comparison = comparison_name,
        x = x_col,
        y = y_col,
        n_samples = nrow(df_test),
        spearman_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
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

make_corr_label <- function(rho, p) {
  paste0(
    "Spearman \u03c1 = ", round(rho, 3),
    "\n", format_p_value(p, prefix = "p")
  )
}

###############################################################################
# 3) HARMONIZE LABELS AND ORDER
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

signature_order <- c(
  "MitoAll",
  "MitoOnly",
  "mtRNA_13PCG"
)

missing_signatures <- setdiff(signature_order, unique(scores$signature))

if (length(missing_signatures) > 0) {
  stop(
    "The following signatures are missing from the score file: ",
    paste(missing_signatures, collapse = ", "),
    "\nPlease rerun Step 1 and Step 2 after adding mtRNA_13PCG."
  )
}

scores <- scores %>%
  filter(signature %in% signature_order) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    signature = factor(signature, levels = signature_order)
  ) %>%
  filter(
    !is.na(tissue),
    !is.na(signature_score)
  )

scores_wide <- scores %>%
  select(tissue, sample_id, signature, signature_score) %>%
  distinct() %>%
  pivot_wider(
    names_from = signature,
    values_from = signature_score
  ) %>%
  filter(
    !is.na(MitoAll),
    !is.na(MitoOnly),
    !is.na(mtRNA_13PCG)
  )

###############################################################################
# 4) COLORS AND THEME
###############################################################################

tissue_palette <- c(
  "Whole Blood" = "#B2182B",
  "Heart - Atrial Appendage" = "#2166AC",
  "Heart - Left Ventricle" = "#1B9E77",
  "Muscle - Skeletal" = "#D95F02"
)

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

###############################################################################
# 5) POOLED CORRELATIONS
###############################################################################

pooled_correlations <- bind_rows(
  run_spearman(
    scores_wide,
    "MitoAll",
    "mtRNA_13PCG",
    "Pooled_MitoAll_vs_mtRNA_13PCG"
  ),
  run_spearman(
    scores_wide,
    "MitoOnly",
    "mtRNA_13PCG",
    "Pooled_MitoOnly_vs_mtRNA_13PCG"
  )
)

write_csv(
  pooled_correlations,
  file.path(FIG_DIR, "Figure2_pooled_spearman_correlations.csv")
)

pooled_label_mitoall <- pooled_correlations %>%
  filter(comparison == "Pooled_MitoAll_vs_mtRNA_13PCG") %>%
  mutate(label = make_corr_label(spearman_rho, p_value)) %>%
  pull(label)

pooled_label_mitoonly <- pooled_correlations %>%
  filter(comparison == "Pooled_MitoOnly_vs_mtRNA_13PCG") %>%
  mutate(label = make_corr_label(spearman_rho, p_value)) %>%
  pull(label)

###############################################################################
# 6) PANEL A — POOLED MITOALL VS mtRNA
###############################################################################

plot_A <- ggplot(
  scores_wide,
  aes(x = MitoAll, y = mtRNA_13PCG, color = tissue)
) +
  geom_point(size = 1.2, alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.65) +
  annotate(
    "text",
    x = quantile(scores_wide$MitoAll, 0.05, na.rm = TRUE),
    y = quantile(scores_wide$mtRNA_13PCG, 0.95, na.rm = TRUE),
    label = pooled_label_mitoall,
    hjust = 0,
    vjust = 1,
    size = 3.6
  ) +
  scale_color_manual(values = tissue_palette, labels = tissue_labels) +
  labs(
    title = "Pooled MitoAll vs mtRNA 13 PCGs",
    x = "MitoAll score",
    y = "mtRNA 13 PCGs score",
    color = "Tissue"
  ) +
  theme_pub +
  theme(
    legend.position = "right"
  )

###############################################################################
# 7) PANEL B — POOLED MITOONLY VS mtRNA
###############################################################################

plot_B <- ggplot(
  scores_wide,
  aes(x = MitoOnly, y = mtRNA_13PCG, color = tissue)
) +
  geom_point(size = 1.2, alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.65) +
  annotate(
    "text",
    x = quantile(scores_wide$MitoOnly, 0.05, na.rm = TRUE),
    y = quantile(scores_wide$mtRNA_13PCG, 0.95, na.rm = TRUE),
    label = pooled_label_mitoonly,
    hjust = 0,
    vjust = 1,
    size = 3.6
  ) +
  scale_color_manual(values = tissue_palette, labels = tissue_labels) +
  labs(
    title = "Pooled MitoOnly vs mtRNA 13 PCGs",
    x = "MitoOnly score",
    y = "mtRNA 13 PCGs score",
    color = "Tissue"
  ) +
  theme_pub +
  theme(
    legend.position = "right"
  )

###############################################################################
# 8) TISSUE-SPECIFIC CORRELATIONS
###############################################################################

tissue_specific_correlations <- scores_wide %>%
  group_by(tissue) %>%
  group_modify(~ bind_rows(
    run_spearman(
      .x,
      "MitoAll",
      "mtRNA_13PCG",
      "MitoAll_vs_mtRNA_13PCG"
    ),
    run_spearman(
      .x,
      "MitoOnly",
      "mtRNA_13PCG",
      "MitoOnly_vs_mtRNA_13PCG"
    )
  )) %>%
  ungroup() %>%
  group_by(comparison) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup()

write_csv(
  tissue_specific_correlations,
  file.path(FIG_DIR, "Figure2_tissue_specific_spearman_correlations.csv")
)

make_tissue_label_data <- function(x_col, y_col, comparison_name) {
  
  label_positions <- scores_wide %>%
    group_by(tissue) %>%
    summarise(
      x_pos = quantile(.data[[x_col]], 0.05, na.rm = TRUE),
      y_pos = quantile(.data[[y_col]], 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  tissue_specific_correlations %>%
    filter(comparison == comparison_name) %>%
    left_join(label_positions, by = "tissue") %>%
    mutate(
      label = make_corr_label(spearman_rho, p_value)
    )
}

label_mitoall_by_tissue <- make_tissue_label_data(
  x_col = "MitoAll",
  y_col = "mtRNA_13PCG",
  comparison_name = "MitoAll_vs_mtRNA_13PCG"
)

label_mitoonly_by_tissue <- make_tissue_label_data(
  x_col = "MitoOnly",
  y_col = "mtRNA_13PCG",
  comparison_name = "MitoOnly_vs_mtRNA_13PCG"
)

###############################################################################
# 9) PANEL C — WITHIN-TISSUE MITOALL VS mtRNA
###############################################################################

plot_C <- ggplot(
  scores_wide,
  aes(x = MitoAll, y = mtRNA_13PCG, color = tissue)
) +
  geom_point(size = 0.9, alpha = 0.42) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
  geom_text(
    data = label_mitoall_by_tissue,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0
  ) +
  facet_wrap(~tissue, scales = "free", nrow = 1, labeller = as_labeller(tissue_labels)) +
  scale_color_manual(values = tissue_palette, labels = tissue_labels) +
  labs(
    title = "Within-tissue MitoAll vs mtRNA 13 PCGs",
    x = "MitoAll score",
    y = "mtRNA 13 PCGs score"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

###############################################################################
# 10) PANEL D — WITHIN-TISSUE MITOONLY VS mtRNA
###############################################################################

plot_D <- ggplot(
  scores_wide,
  aes(x = MitoOnly, y = mtRNA_13PCG, color = tissue)
) +
  geom_point(size = 0.9, alpha = 0.42) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
  geom_text(
    data = label_mitoonly_by_tissue,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0
  ) +
  facet_wrap(~tissue, scales = "free", nrow = 1, labeller = as_labeller(tissue_labels)) +
  scale_color_manual(values = tissue_palette, labels = tissue_labels) +
  labs(
    title = "Within-tissue MitoOnly vs mtRNA 13 PCGs",
    x = "MitoOnly score",
    y = "mtRNA 13 PCGs score"
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

###############################################################################
# 11) COMBINE FIGURE 2
###############################################################################

figure_2 <-
  (plot_A | plot_B) /
  plot_C /
  plot_D +
  plot_layout(
    heights = c(1, 1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

###############################################################################
# 12) SAVE FIGURE 2
###############################################################################

ggsave(
  filename = file.path(FIG_DIR, "Figure2_MitoSignatures_mtRNA_Correlations.pdf"),
  plot = figure_2,
  width = 16,
  height = 15,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(FIG_DIR, "Figure2_MitoSignatures_mtRNA_Correlations.tiff"),
  plot = figure_2,
  width = 16,
  height = 15,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(FIG_DIR, "Figure2_MitoSignatures_mtRNA_Correlations.png"),
  plot = figure_2,
  width = 16,
  height = 15,
  units = "in",
  dpi = 600
)

print(figure_2)










