#' ---
#' title: "DAP-seq SE and PE comparison"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     code_folding: hide
#' ---
#'
#' # Description
#'
#' Comparison of DAP-seq sample statistics between single-end and paired-end
#' processing runs. The analysis compares read survival, alignment, duplicate
#' rate, FRiP, peak counts, promoter intersections, and motif detection.

# ----

#' # Libraries and functions
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

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

frip_to_percent <- function(x) {
  x <- as.numeric(x)
  
  if (max(x, na.rm = TRUE) <= 1) {
    x * 100
  } else {
    x
  }
}

clean_motif_status <- function(x) {
  case_when(
    tolower(x) == "yes" ~ "yes",
    TRUE ~ "no"
  )
}

make_numeric_stats <- function(stats, mode) {
  if (mode == "PE") {
    stats <- stats %>%
      transmute(
        sample,
        trimmomatic_surviving_pct =
          100 * as.numeric(trimmomatic_both_surviving) /
          as.numeric(trimmomatic_input_pairs),
        bowtie2_overall_alignment_rate_pct =
          as.numeric(bowtie2_overall_alignment_rate),
        bowtie2_unique_alignment_pct =
          as.numeric(bowtie2_concordant_1_pct),
        pcr_duplicate_pct =
          100 * as.numeric(samtools_markdup_duplicate_total) /
          as.numeric(samtools_markdup_examined),
        frip_pct =
          frip_to_percent(FRiP),
        macs3_peak_count =
          as.numeric(macs3_peak_count),
        promoter_intersect_count =
          as.numeric(promoter_intersect_count)
      )
  }
  
  if (mode == "SE") {
    stats <- stats %>%
      transmute(
        sample,
        trimmomatic_surviving_pct =
          100 * as.numeric(trimmomatic_surviving_reads) /
          as.numeric(trimmomatic_input_reads),
        bowtie2_overall_alignment_rate_pct =
          as.numeric(bowtie2_overall_alignment_rate),
        bowtie2_unique_alignment_pct =
          as.numeric(bowtie2_aligned_1_pct),
        pcr_duplicate_pct =
          100 * as.numeric(samtools_markdup_duplicate_total) /
          as.numeric(samtools_markdup_examined),
        frip_pct =
          frip_to_percent(FRiP),
        macs3_peak_count =
          as.numeric(macs3_peak_count),
        promoter_intersect_count =
          as.numeric(promoter_intersect_count)
      )
  }
  
  stats %>%
    pivot_longer(
      cols = -sample,
      names_to = "metric",
      values_to = mode
    )
}

make_motif_stats <- function(stats, mode) {
  stats %>%
    transmute(
      sample,
      meme_has_motif = clean_motif_status(meme_has_motif),
      centrimo_has_central_motif = clean_motif_status(centrimo_has_central_motif)
    ) %>%
    pivot_longer(
      cols = -sample,
      names_to = "metric",
      values_to = mode
    )
}

make_metric_plot <- function(compare_table,
                             metric_name,
                             plot_title,
                             left_label,
                             diff_column,
                             right_label) {
  metric_table <- compare_table %>%
    filter(metric == metric_name)
  
  left_panel <- metric_table %>%
    select(sample, SE, PE) %>%
    pivot_longer(
      cols = c(SE, PE),
      names_to = "mode",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(
      panel = paste0("SE / PE\n", left_label),
      x_value = factor(mode, levels = c("SE", "PE"))
    )
  
  right_panel <- metric_table %>%
    transmute(
      sample,
      value = .data[[diff_column]],
      panel = paste0("Difference\n", right_label),
      x_value = factor("Difference", levels = "Difference")
    ) %>%
    filter(!is.na(value))
  
  panel_levels <- c(
    paste0("SE / PE\n", left_label),
    paste0("Difference\n", right_label)
  )
  
  plot_table <- bind_rows(left_panel, right_panel) %>%
    mutate(panel = factor(panel, levels = panel_levels))
  
  hline_data <- data.frame(
    panel = factor(paste0("Difference\n", right_label), levels = panel_levels),
    yintercept = 0
  )
  
  ggplot(plot_table, aes(x = x_value, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(
      data = hline_data,
      aes(yintercept = yintercept),
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ panel,
      scales = "free",
      nrow = 1
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10)
    )
}

make_motif_plot <- function(compare_table, metric_name, plot_title) {
  motif_table <- compare_table %>%
    filter(metric == metric_name) %>%
    select(sample, SE, PE) %>%
    pivot_longer(
      cols = c(SE, PE),
      names_to = "mode",
      values_to = "status"
    ) %>%
    mutate(
      mode = factor(mode, levels = c("SE", "PE")),
      status = factor(status, levels = c("yes", "no"))
    )
  
  motif_summary <- motif_table %>%
    count(mode, status, .drop = FALSE) %>%
    group_by(mode) %>%
    mutate(
      pct = 100 * n / sum(n),
      label = paste0(n, " (", round(pct, 1), "%)")
    ) %>%
    ungroup()
  
  ggplot(motif_summary, aes(x = mode, y = pct, fill = status)) +
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
    theme_bw()
}

#' # Data

se_stats_dir <- "/mnt/ada/projects/spruce/nstreet/spruce_chromatin_dynamics/DAPseq_test_SE/nf_out/sample_stats"
pe_stats_dir <- "/mnt/ada/projects/spruce/nstreet/spruce_chromatin_dynamics/DAPseq_test_PE/nf_out/sample_stats"

se_stats <- read_sample_stats(se_stats_dir)
pe_stats <- read_sample_stats(pe_stats_dir)

se_wide <- make_wide_stats(se_stats)
pe_wide <- make_wide_stats(pe_stats)

#' # Analysis

#' Numeric sample statistics
pe_values <- make_numeric_stats(pe_wide, "PE")
se_values <- make_numeric_stats(se_wide, "SE")

pe_se_compare <- pe_values %>%
  left_join(se_values, by = c("sample", "metric")) %>%
  mutate(
    diff = SE - PE,
    pct_change = if_else(
      PE == 0 | is.na(PE),
      NA_real_,
      100 * (SE - PE) / PE
    )
  )

#' Motif detection statistics
pe_motifs <- make_motif_stats(pe_wide, "PE")
se_motifs <- make_motif_stats(se_wide, "SE")

pe_se_motifs <- pe_motifs %>%
  left_join(se_motifs, by = c("sample", "metric")) %>%
  mutate(same_status = PE == SE)

#' # Plots

#' Trimmomatic surviving reads percent
trimmomatic_surviving_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "trimmomatic_surviving_pct",
  plot_title = "Trimmomatic surviving reads",
  left_label = "Surviving reads (%)",
  diff_column = "diff",
  right_label = "SE - PE percentage points"
)

#' Bowtie2 overall alignment percent
bowtie2_overall_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "bowtie2_overall_alignment_rate_pct",
  plot_title = "Bowtie2 overall alignment rate",
  left_label = "Overall alignment (%)",
  diff_column = "diff",
  right_label = "SE - PE percentage points"
)

#' Bowtie2 unique alignment percent
bowtie2_unique_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "bowtie2_unique_alignment_pct",
  plot_title = "Bowtie2 unique alignment",
  left_label = "Unique alignment (%)",
  diff_column = "diff",
  right_label = "SE - PE percentage points"
)

#' PCR duplicate percent
pcr_duplicate_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "pcr_duplicate_pct",
  plot_title = "PCR duplicate rate",
  left_label = "PCR duplicates (%)",
  diff_column = "diff",
  right_label = "SE - PE percentage points"
)

#' FRiP percent
frip_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "frip_pct",
  plot_title = "FRiP",
  left_label = "FRiP (%)",
  diff_column = "diff",
  right_label = "SE - PE percentage points"
)

#' MACS3 peak count
peak_count_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "macs3_peak_count",
  plot_title = "MACS3 peak count",
  left_label = "Peak count",
  diff_column = "pct_change",
  right_label = "Relative difference (%)"
)

#' Promoter intersect count
promoter_intersect_plot <- make_metric_plot(
  compare_table = pe_se_compare,
  metric_name = "promoter_intersect_count",
  plot_title = "Promoter intersect count",
  left_label = "Promoter intersect count",
  diff_column = "pct_change",
  right_label = "Relative difference (%)"
)

#' MEME motif status
meme_has_motif_plot <- make_motif_plot(
  compare_table = pe_se_motifs,
  metric_name = "meme_has_motif",
  plot_title = "MEME motif detected"
)

#' CentriMo motif status
centrimo_has_central_motif_plot <- make_motif_plot(
  compare_table = pe_se_motifs,
  metric_name = "centrimo_has_central_motif",
  plot_title = "CentriMo centrally enriched motif detected"
)

#' # Tables and plots

pe_se_compare
pe_se_motifs

trimmomatic_surviving_plot
bowtie2_overall_plot
bowtie2_unique_plot
pcr_duplicate_plot
frip_plot
peak_count_plot
promoter_intersect_plot
meme_has_motif_plot
centrimo_has_central_motif_plot

