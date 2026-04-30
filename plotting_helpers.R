# ===== PLOTTING HELPERS =====
#
# Functions for generating QQ and volcano plots from FAST_omics_WAS results.
# Automatically detects strata, FU levels, and DNAm filtered probes.
#
# Usage:
#   source("plotting_helpers.R")
#   results <- FAST_omics_WAS(pheno, omics)
#   generate_all_plots(results, analysis = "analysis_change")
#   generate_all_plots(results$analysis_change)  # equivalent

library(qqman)
library(ggplot2)

# Accept either:
# - a stratum-level results object (list(all/male/female)), OR
# - the full FAST_omics_WAS return value (list(analysis_change=<...>, analysis_level=<...>)).
.coerce_to_strata_results <- function(results, analysis = c("analysis_change", "analysis_level")) {
  analysis <- match.arg(analysis)

  if (is.list(results) && all(c("all", "male", "female") %in% names(results))) {
    return(results)
  }

  if (is.list(results) && analysis %in% names(results)) {
    return(results[[analysis]])
  }

  stop("Invalid results object. Pass FAST_omics_WAS output (with $analysis_change/$analysis_level) or a stratum-level results list.")
}

# ===== INDIVIDUAL PLOT FUNCTIONS =====

#' QQ plot of p-values using qqman
#'
#' @param results Results object from FAST_omics_WAS
#' @param stratum "all", "male", or "female"
#' @param fu Which FU level to plot
#' @param probe_set "full" (all probes) or "filtered" (filtered probes only, DNAm)
#' @param analysis If passing the full FAST_omics_WAS output, which analysis to plot ("analysis_change" or "analysis_level")
#' @param title Optional custom title (auto-generated if NULL)
plot_qq <- function(results, stratum = "all", fu = 1, probe_set = "full",
                    analysis = c("analysis_change", "analysis_level"), title = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  te <- results[[stratum]]$treatment_effects
  if (is.null(te)) stop(paste0("No results for stratum '", stratum, "'"))

  # Subset to requested FU
  te <- te[te$FU == fu, ]
  if (nrow(te) == 0) stop(paste0("No results for FU=", fu))

  # Subset to probe set
  if (probe_set == "filtered") {
    te <- te[!is.na(te$BH_P_VALUE_FILTERED), ]
    if (nrow(te) == 0) stop("No filtered probes found")
  }

  pvals <- te$P_VALUE[!is.na(te$P_VALUE)]

  # Auto title
  if (is.null(title)) {
    title <- paste0("QQ Plot - ", stratum, ", FU=", fu, ", ", probe_set)
  }

  qq(pvals, main = title)
}

#' Volcano plot of effect size vs -log10 p-value
#'
#' @param results Results object from FAST_omics_WAS
#' @param stratum "all", "male", or "female"
#' @param fu Which FU level to plot
#' @param probe_set "full" (all probes) or "filtered" (filtered probes only, DNAm)
#' @param p_threshold Significance threshold for coloring (default 0.05, applied to BH)
#' @param analysis If passing the full FAST_omics_WAS output, which analysis to plot ("analysis_change" or "analysis_level")
#' @param title Optional custom title (auto-generated if NULL)
plot_volcano <- function(results, stratum = "all", fu = 1, probe_set = "full",
                         p_threshold = 0.05, analysis = c("analysis_change", "analysis_level"), title = NULL) {
  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  te <- results[[stratum]]$treatment_effects
  if (is.null(te)) stop(paste0("No results for stratum '", stratum, "'"))
  te <- te[te$FU == fu, ]
  if (nrow(te) == 0) stop(paste0("No results for FU=", fu))
  if (probe_set == "filtered") {
    te <- te[!is.na(te$BH_P_VALUE_FILTERED), ]
    if (nrow(te) == 0) stop("No filtered probes found")
  }

  bh_col <- if (probe_set == "filtered") "BH_P_VALUE_FILTERED" else "BH_P_VALUE"
  sig <- !is.na(te[[bh_col]]) & te[[bh_col]] < p_threshold

  te$neg_log_p <- -log10(te$P_VALUE)
  te$direction <- ifelse(!sig, "Not significant",
                         ifelse(te$EFFECT_SIZE > 0, "Up", "Down"))

  all_levels <- c("Down", "Not significant", "Up")
  te$direction <- factor(te$direction, levels = all_levels)

  for (lvl in setdiff(all_levels, as.character(unique(te$direction)))) {
    phantom <- te[NA, ]  # one-row NA frame with all columns of te
    phantom$direction <- factor(lvl, levels = all_levels)
    te <- rbind(te, phantom)
  }

  n_up   <- sum(te$direction == "Up",   na.rm = TRUE)
  n_down <- sum(te$direction == "Down", na.rm = TRUE)

  color_values <- c("Down" = "blue", "Not significant" = "grey50", "Up" = "red")
  color_labels <- c(
    "Down"            = paste0("Down"),
    "Not significant" = "Not significant",
    "Up"              = paste0("Up")
  )

  if (is.null(title)) {
    title <- paste0("Volcano Plot - ", stratum, ", FU=", fu, ", ", probe_set)
  }

  ggplot(te, aes(x = EFFECT_SIZE, y = neg_log_p, color = direction)) +
    geom_point(size = 1, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = color_values, labels = color_labels) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Effect Size", y = expression(-log[10](p)),
         title = title, color = NULL) +
    theme_minimal() +
    theme(legend.position = "right")
}

# ===== BATCH GENERATION =====

#' Generate all QQ and volcano plots for a results object
#'
#' Detects strata, FU levels, and DNAm (filtered probes) automatically.
#' Outputs one PDF per plot type x probe set combination.
#'
#' @param results Results object from FAST_omics_WAS
#' @param figures_dir Directory to save PDFs (default "Figures", created if needed)
generate_all_plots <- function(results, figures_dir = NULL, analysis = c("analysis_change", "analysis_level")) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  # Default directory with user prompt
  if (is.null(figures_dir)) {
    figures_dir <- "Figures"
    if (!dir.exists(figures_dir)) {
      response <- readline(paste0("'", figures_dir, "' does not exist. Create it? (y/n): "))
      if (!tolower(response) %in% c("y", "yes")) {
        stop("No output directory. Provide figures_dir or create 'Figures/'.")
      }
      dir.create(figures_dir)
    }
  } else if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }

  # Detect dimensions
  strata <- c("all", "male", "female")[!sapply(results[c("all", "male", "female")], is.null)]
  fu_levels <- sort(unique(results[[strata[1]]]$treatment_effects$FU))

  has_filtered <- "BH_P_VALUE_FILTERED" %in% colnames(results[[strata[1]]]$treatment_effects) &&
    any(!is.na(results[[strata[1]]]$treatment_effects$BH_P_VALUE_FILTERED))

  probe_sets <- if (has_filtered) c("full", "filtered") else c("full")

  cat("Generating plots:\n")
  cat("  Strata:     ", paste(strata, collapse = ", "), "\n")
  cat("  FU levels:  ", paste(fu_levels, collapse = ", "), "\n")
  cat("  Probe sets: ", paste(probe_sets, collapse = ", "), "\n")

  for (ps in probe_sets) {

    # QQ plots — one PDF per probe set, one page per stratum x FU
    qq_file <- file.path(figures_dir, paste0("qq_", ps, ".pdf"))
    pdf(qq_file, width = 7, height = 6)
    for (s in strata) {
      for (fu in fu_levels) {
        plot_qq(results, stratum = s, fu = fu, probe_set = ps)
      }
    }
    dev.off()
    cat("  Saved:", qq_file, "\n")

    # Volcano plots — one PDF per probe set, one page per stratum x FU
    volcano_file <- file.path(figures_dir, paste0("volcano_", ps, ".pdf"))
    pdf(volcano_file, width = 8, height = 6)
    for (s in strata) {
      for (fu in fu_levels) {
        print(plot_volcano(results, stratum = s, fu = fu, probe_set = ps))
      }
    }
    dev.off()
    cat("  Saved:", volcano_file, "\n")
  }

  cat("Done.\n")
}
