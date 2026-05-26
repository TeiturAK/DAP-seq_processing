#' ---
#' title: "DAP-seq O'Malley Arabidopsis and Rodriguez Aspen comparison"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     code_folding: hide
#'     fig_width: 7
#'     fig_height: 9
#' ---
#'
#' # Description
#'
#' Comparison of DAP-seq sample statistics between the O'Malley Arabidopsis
#' dataset and the Rodriguez Aspen dataset. The analysis compares input library
#' size, read survival, alignment, duplicate rate, FRiP, peak counts, promoter
#' intersections, and motif detection. CentriMo motif detection is called using
#' a top motif E-value cutoff.

# ----

#' # Libraries and functions
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
})

centrimo_evalue_cutoff <- 0.05

read_sample_stats <- function(stats_dir) {
  list.files(
    path = stats_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  ) %>%
    lapply(read_tsv, show_col_types = FALSE) %>%
    bind_rows()
}

make_wide_stats <- function(stats) {
  stats %>%
    pivot_wider(
      id_cols = sample,
      names_from = metric,
      values_from = value
    )
}

get_numeric_column <- function(df, column_name) {
  if (column_name %in% colnames(df)) {
    as.numeric(df[[column_name]])
  } else {
    rep(NA_real_, nrow(df))
  }
}

get_character_column <- function(df, column_name) {
  if (column_name %in% colnames(df)) {
    as.character(df[[column_name]])
  } else {
    rep(NA_character_, nrow(df))
  }
}

frip_to_percent <- function(x) {
  x <- as.numeric(x)
  
  if (max(x, na.rm = TRUE) <= 1) {
    x * 100
  } else {
    x
  }
}

clean_meme_status <- function(x) {
  case_when(
    tolower(x) == "yes" ~ "yes",
    TRUE ~ "no"
  )
}

clean_centrimo_status <- function(evalue, cutoff = centrimo_evalue_cutoff) {
  case_when(
    !is.na(evalue) & evalue < cutoff ~ "yes",
    TRUE ~ "no"
  )
}

make_numeric_stats <- function(stats, dataset_name) {
  trimmomatic_pe_surviving <- get_numeric_column(stats, "trimmomatic_both_surviving")
  trimmomatic_pe_input <- get_numeric_column(stats, "trimmomatic_input_pairs")
  
  trimmomatic_se_surviving <- get_numeric_column(stats, "trimmomatic_surviving_reads")
  trimmomatic_se_input <- get_numeric_column(stats, "trimmomatic_input_reads")
  
  bowtie2_pe_unique <- get_numeric_column(stats, "bowtie2_concordant_1_pct")
  bowtie2_se_unique <- get_numeric_column(stats, "bowtie2_aligned_1_pct")
  
  stats %>%
    transmute(
      sample,
      dataset = dataset_name,
      
      library_size_fragments = coalesce(
        trimmomatic_pe_input,
        trimmomatic_se_input
      ),
      
      trimmomatic_surviving_pct = coalesce(
        100 * trimmomatic_pe_surviving / trimmomatic_pe_input,
        100 * trimmomatic_se_surviving / trimmomatic_se_input
      ),
      
      bowtie2_overall_alignment_rate_pct =
        get_numeric_column(stats, "bowtie2_overall_alignment_rate"),
      
      bowtie2_unique_alignment_pct =
        coalesce(bowtie2_pe_unique, bowtie2_se_unique),
      
      pcr_duplicate_pct =
        100 * get_numeric_column(stats, "samtools_markdup_duplicate_total") /
        get_numeric_column(stats, "samtools_markdup_examined"),
      
      frip_pct =
        frip_to_percent(get_numeric_column(stats, "FRiP")),
      
      macs3_peak_count =
        get_numeric_column(stats, "macs3_peak_count"),
      
      promoter_intersect_count =
        get_numeric_column(stats, "promoter_intersect_count")
    ) %>%
    pivot_longer(
      cols = -c(sample, dataset),
      names_to = "metric",
      values_to = "value"
    )
}

make_motif_stats <- function(stats, dataset_name) {
  stats %>%
    transmute(
      sample,
      dataset = dataset_name,
      meme_has_motif = clean_meme_status(
        get_character_column(stats, "meme_has_motif")
      ),
      centrimo_has_central_motif = clean_centrimo_status(
        get_numeric_column(stats, "centrimo_top_evalue")
      )
    ) %>%
    pivot_longer(
      cols = c(meme_has_motif, centrimo_has_central_motif),
      names_to = "metric",
      values_to = "status"
    )
}

# make_metric_plot <- function(compare_table, metric_name, plot_title, ylab) {
#   plot_table <- compare_table %>%
#     filter(metric == metric_name) %>%
#     filter(!is.na(value)) %>%
#     mutate(
#       dataset = factor(
#         dataset,
#         levels = c("O'Malley Arabidopsis", "Rodriguez Aspen")
#       )
#     )
#   
#   ggplot(plot_table, aes(x = dataset, y = value, fill = dataset)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
#     labs(
#       title = plot_title,
#       x = NULL,
#       y = ylab,
#       fill = NULL
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     )
# }

# make_metric_plot <- function(compare_table, metric_name, plot_title, ylab) {
#   plot_table <- compare_table %>%
#     filter(metric == metric_name) %>%
#     filter(!is.na(value)) %>%
#     mutate(
#       dataset = factor(
#         dataset,
#         levels = c("O'Malley Arabidopsis", "Rodriguez Aspen")
#       )
#     )
#   
#   same_scale_plot <- ggplot(plot_table, aes(x = dataset, y = value, fill = dataset)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
#     labs(
#       title = paste0(plot_title, "\nSame y-axis scale"),
#       x = NULL,
#       y = ylab,
#       fill = NULL
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     )
#   
#   free_scale_plot <- ggplot(plot_table, aes(x = dataset, y = value, fill = dataset)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
#     facet_wrap(
#       ~ dataset,
#       scales = "free_y",
#       nrow = 1
#     ) +
#     labs(
#       title = paste0(plot_title, "\nSeparate y-axis scales"),
#       x = NULL,
#       y = ylab,
#       fill = NULL
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       legend.position = "none"
#     )
#   
#   plot_grid(
#     same_scale_plot,
#     free_scale_plot,
#     ncol = 1,
#     align = "v"
#   )
# }

make_metric_plot <- function(compare_table, metric_name, plot_title, ylab) {
  plot_table <- compare_table %>%
    filter(metric == metric_name) %>%
    filter(!is.na(value)) %>%
    mutate(
      dataset = factor(
        dataset,
        levels = c("O'Malley Arabidopsis", "Rodriguez Aspen")
      )
    )
  
  same_scale_plot <- ggplot(plot_table, aes(x = dataset, y = value, fill = dataset)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(
      title = paste0(plot_title, "\nSame y-axis scale"),
      x = NULL,
      y = ylab,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  free_scale_plot <- ggplot(plot_table, aes(x = "Samples", y = value, fill = dataset)) +
    geom_boxplot(outlier.shape = NA, width = 0.45) +
    geom_jitter(width = 0.08, alpha = 0.4, size = 1) +
    facet_wrap(
      ~ dataset,
      scales = "free_y",
      nrow = 1
    ) +
    labs(
      title = paste0(plot_title, "\nSeparate y-axis scales"),
      x = NULL,
      y = ylab,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  # plot_grid(
  #   same_scale_plot,
  #   free_scale_plot,
  #   ncol = 1,
  #   align = "v",
  #   rel_heights = c(1.2, 1)
  # )
  plot_grid(
    same_scale_plot,
    free_scale_plot,
    ncol = 1,
    rel_heights = c(1.2, 1)
  )
}

make_motif_plot <- function(motif_table, metric_name, plot_title) {
  plot_table <- motif_table %>%
    filter(metric == metric_name) %>%
    mutate(
      dataset = factor(
        dataset,
        levels = c("O'Malley Arabidopsis", "Rodriguez Aspen")
      ),
      status = factor(status, levels = c("yes", "no"))
    )
  
  motif_summary <- plot_table %>%
    count(dataset, status, .drop = FALSE) %>%
    group_by(dataset) %>%
    mutate(
      pct = 100 * n / sum(n),
      label = paste0(n, " (", round(pct, 1), "%)")
    ) %>%
    ungroup()
  
  ggplot(motif_summary, aes(x = dataset, y = pct, fill = status)) +
    geom_col() +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Percent of samples",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' # Data

omalley_stats_dir <- "/mnt/ada/projects/spruce/nstreet/spruce_chromatin_dynamics/DAPseq_bartlett/nf_out/sample_stats"
rodriguez_stats_dir <- "/mnt/ada/projects/spruce/nstreet/spruce_chromatin_dynamics/DAPseq_aspen/nf_out/sample_stats"

omalley_stats <- read_sample_stats(omalley_stats_dir)
rodriguez_stats <- read_sample_stats(rodriguez_stats_dir)

omalley_wide <- make_wide_stats(omalley_stats)
rodriguez_wide <- make_wide_stats(rodriguez_stats)

#' # Analysis

#' Numeric sample statistics
omalley_values <- make_numeric_stats(
  omalley_wide,
  "O'Malley Arabidopsis"
)

rodriguez_values <- make_numeric_stats(
  rodriguez_wide,
  "Rodriguez Aspen"
)

dataset_compare <- bind_rows(
  omalley_values,
  rodriguez_values
)

#' Motif detection statistics
omalley_motifs <- make_motif_stats(
  omalley_wide,
  "O'Malley Arabidopsis"
)

rodriguez_motifs <- make_motif_stats(
  rodriguez_wide,
  "Rodriguez Aspen"
)

dataset_motifs <- bind_rows(
  omalley_motifs,
  rodriguez_motifs
)

#' # Plots

#' Input library size
library_size_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "library_size_fragments",
  plot_title = "Input library size",
  ylab = "Fragments"
)

#' Trimmomatic surviving reads percent
trimmomatic_surviving_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "trimmomatic_surviving_pct",
  plot_title = "Trimmomatic surviving reads",
  ylab = "Surviving reads (%)"
)

#' Bowtie2 overall alignment percent
bowtie2_overall_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "bowtie2_overall_alignment_rate_pct",
  plot_title = "Bowtie2 overall alignment rate",
  ylab = "Overall alignment (%)"
)

#' Bowtie2 unique alignment percent
bowtie2_unique_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "bowtie2_unique_alignment_pct",
  plot_title = "Bowtie2 unique alignment",
  ylab = "Unique alignment (%)"
)

#' PCR duplicate percent
pcr_duplicate_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "pcr_duplicate_pct",
  plot_title = "PCR duplicate rate",
  ylab = "PCR duplicates (%)"
)

#' FRiP percent
frip_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "frip_pct",
  plot_title = "FRiP",
  ylab = "FRiP (%)"
)

#' MACS3 peak count
peak_count_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "macs3_peak_count",
  plot_title = "MACS3 peak count",
  ylab = "Peak count"
)

#' Promoter intersect count
promoter_intersect_plot <- make_metric_plot(
  compare_table = dataset_compare,
  metric_name = "promoter_intersect_count",
  plot_title = "Promoter intersect count",
  ylab = "Promoter intersect count"
)

#' MEME motif status
meme_has_motif_plot <- make_motif_plot(
  motif_table = dataset_motifs,
  metric_name = "meme_has_motif",
  plot_title = "Significant MEME-ChIP motif found"
)

#' CentriMo motif status
centrimo_has_central_motif_plot <- make_motif_plot(
  motif_table = dataset_motifs,
  metric_name = "centrimo_has_central_motif",
  plot_title = paste0("CentriMo centrally enriched motif detected, E-value < ", centrimo_evalue_cutoff)
)

#' # Tables

dataset_compare
dataset_motifs

#' # Library size and read processing

#' ## Input library size
library_size_plot

#' ## Trimmomatic surviving reads
trimmomatic_surviving_plot

#' # Alignment

#' ## Bowtie2 overall alignment
bowtie2_overall_plot

#' ## Bowtie2 unique alignment
bowtie2_unique_plot

#' # Library complexity

#' ## PCR duplicate rate
pcr_duplicate_plot

#' # Peak count and signal strength 

#' ## MACS3 peak count
peak_count_plot

#' ## FRiP
frip_plot

#' # Feature intersect

#' ## Promoter intersect count
promoter_intersect_plot

#' # Motif results

#' ## Significant MEME-ChIP motif found
meme_has_motif_plot

#' ## CentriMo central motif enrichment
centrimo_has_central_motif_plot