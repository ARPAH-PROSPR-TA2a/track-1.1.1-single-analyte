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

.select_treatment_effects <- function(treatment_effects, treatment = NULL) {
  if (is.null(treatment_effects) || nrow(treatment_effects) == 0L) {
    stop("No treatment effects available to plot")
  }
  if (!"CONTRAST" %in% names(treatment_effects)) {
    if (any(c("TREATMENT", "REFERENCE") %in% names(treatment_effects)) ||
        anyDuplicated(treatment_effects[, c("ANALYTE_NAME", "FU"), drop = FALSE])) {
      stop("Missing CONTRAST: treatment labels or repeated analyte/follow-up rows cannot be interpreted as legacy binary results")
    }
    treatment_effects$CONTRAST <- "1_vs_0"
    treatment_effects$TREATMENT <- "Treatment"
    treatment_effects$REFERENCE <- "Control"
  }
  contrasts <- unique(as.character(treatment_effects$CONTRAST))
  if (anyNA(contrasts) || any(!nzchar(contrasts)) ||
      anyDuplicated(treatment_effects[, c("ANALYTE_NAME", "FU", "CONTRAST"), drop = FALSE])) {
    stop("Invalid or duplicate treatment contrast keys")
  }
  if (!all(c("TREATMENT", "REFERENCE") %in% names(treatment_effects))) {
    stop("Contrast results must contain TREATMENT and REFERENCE labels")
  }
  for (column in c("TREATMENT", "REFERENCE")) {
    identities <- split(as.character(treatment_effects[[column]]), treatment_effects$CONTRAST)
    if (any(vapply(identities, function(x) anyNA(x) || length(unique(x)) != 1L, logical(1)))) {
      stop("Each CONTRAST must identify exactly one ", column)
    }
  }
  treatments <- unique(as.character(treatment_effects$TREATMENT))
  if (length(treatments) != length(contrasts)) {
    stop("Each treatment label must identify exactly one CONTRAST")
  }
  if (is.null(treatment)) {
    if (length(treatments) != 1L) {
      stop("Multiple treatments are available. Supply treatment = one of: ",
           paste(treatments, collapse = ", "))
    }
    treatment <- treatments[1L]
  }
  if (!is.character(treatment) || length(treatment) != 1L || is.na(treatment) ||
      !treatment %in% treatments) {
    stop("Unknown treatment. Available treatments: ", paste(treatments, collapse = ", "))
  }
  treatment_effects[treatment_effects$TREATMENT == treatment, , drop = FALSE]
}


.treatment_contrast_label <- function(treatment_effects) {
  if (all(c("TREATMENT", "REFERENCE") %in% names(treatment_effects))) {
    return(paste0(treatment_effects$TREATMENT[1L], " vs ", treatment_effects$REFERENCE[1L]))
  }
  as.character(treatment_effects$CONTRAST[1L])
}

.treatment_filtered_membership <- function(treatment_effects) {
  if (is.null(treatment_effects)) return(logical())
  if ("IN_FILTERED_SET" %in% names(treatment_effects)) {
    return(treatment_effects$IN_FILTERED_SET %in% TRUE)
  }
  # Legacy results encoded membership only through available filtered inference.
  if ("BH_P_VALUE_FILTERED" %in% names(treatment_effects)) {
    p <- treatment_effects$BH_P_VALUE_FILTERED
    return(is.finite(p) & p >= 0 & p <= 1)
  }
  rep(FALSE, nrow(treatment_effects))
}

.treatment_plot_estimable <- function(treatment_effects) {
  keep <- is.finite(treatment_effects$P_VALUE) & is.finite(treatment_effects$EFFECT_SIZE)
  if ("STATUS" %in% names(treatment_effects)) {
    keep <- keep & treatment_effects$STATUS %in% "ok"
  }
  keep
}

.treatment_empty_plot_label <- function(treatment_effects, probe_set) {
  if (nrow(treatment_effects) == 0L && probe_set == "filtered") {
    "No probes in filtered set"
  } else {
    "No estimable results"
  }
}

#' QQ plot of p-values using qqman
#'
#' @param results Results object from FAST_omics_WAS
#' @param stratum "all", "male", or "female"
#' @param fu Which FU level to plot
#' @param probe_set "full" (all probes) or "filtered" (filtered probes only, DNAm)
#' @param analysis If passing the full FAST_omics_WAS output, which analysis to plot ("analysis_change" or "analysis_level")
#' @param title Optional custom title (auto-generated if NULL)
#' @param treatment Treatment label; required when multiple treatments are present
plot_qq <- function(results, stratum = "all", fu = 1, probe_set = "full",
                    analysis = c("analysis_change", "analysis_level"), title = NULL,
                    treatment = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  te <- results[[stratum]]$treatment_effects
  if (is.null(te)) stop(paste0("No results for stratum '", stratum, "'"))
  te <- .select_treatment_effects(te, treatment)

  # Subset to requested FU
  te <- te[te$FU == fu, ]
  if (nrow(te) == 0) stop(paste0("No results for FU=", fu))

  # Preserve the contrast label even when the selected probe set has no members.
  if (is.null(title)) {
    title <- paste0("QQ Plot - ", .treatment_contrast_label(te), ", ",
                    stratum, ", FU=", fu, ", ", probe_set)
  }

  # Subset to probe set
  if (probe_set == "filtered") {
    te <- te[.treatment_filtered_membership(te), , drop = FALSE]
  }

  pvals <- te$P_VALUE[.treatment_plot_estimable(te)]
  if (length(pvals) == 0L) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, .treatment_empty_plot_label(te, probe_set))
    return(invisible(NULL))
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
#' @param treatment Treatment label; required when multiple treatments are present
plot_volcano <- function(results, stratum = "all", fu = 1, probe_set = "full",
                         p_threshold = 0.05, analysis = c("analysis_change", "analysis_level"), title = NULL,
                         treatment = NULL) {
  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  te <- results[[stratum]]$treatment_effects
  if (is.null(te)) stop(paste0("No results for stratum '", stratum, "'"))
  te <- .select_treatment_effects(te, treatment)
  te <- te[te$FU == fu, ]
  if (nrow(te) == 0) stop(paste0("No results for FU=", fu))
  if (is.null(title)) {
    title <- paste0("Volcano Plot - ", .treatment_contrast_label(te), ", ",
                    stratum, ", FU=", fu, ", ", probe_set)
  }
  if (probe_set == "filtered") {
    te <- te[.treatment_filtered_membership(te), , drop = FALSE]
  }
  keep <- .treatment_plot_estimable(te)
  if (!any(keep)) {
    return(ggplot() +
      annotate("text", x = 0, y = 0, label = .treatment_empty_plot_label(te, probe_set)) +
      labs(title = title) + theme_void())
  }
  te <- te[keep, , drop = FALSE]

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
#' Detects strata, FU levels, treatment contrasts, and DNAm (filtered probes) automatically.
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

  # Use the cells actually present in each stratum; arms or visits can be absent.
  strata <- c("all", "male", "female")[!sapply(results[c("all", "male", "female")], is.null)]
  plot_cells <- lapply(strata, function(stratum) {
    te <- results[[stratum]]$treatment_effects
    if (is.null(te) || nrow(te) == 0L) return(NULL)
    if (!"CONTRAST" %in% names(te)) {
      # Validate legacy rows before assigning their binary label.
      te <- .select_treatment_effects(te)
    }
    if (!"TREATMENT" %in% names(te)) stop("Contrast results must contain TREATMENT labels")
    cells <- unique(te[, c("FU", "CONTRAST", "TREATMENT"), drop = FALSE])
    cells$STRATUM <- stratum
    cells
  })
  plot_cells <- do.call(rbind, plot_cells)
  if (is.null(plot_cells) || nrow(plot_cells) == 0L) stop("No treatment effects available to plot")
  has_filtered <- any(vapply(strata, function(stratum) {
    te <- results[[stratum]]$treatment_effects
    "IN_FILTERED_SET" %in% names(te) || any(.treatment_filtered_membership(te))
  }, logical(1)))
  probe_sets <- if (has_filtered) c("full", "filtered") else c("full")

  cat("Generating plots:\n")
  cat("  Strata:     ", paste(strata, collapse = ", "), "\n")
  cat("  FU levels:  ", paste(sort(unique(plot_cells$FU)), collapse = ", "), "\n")
  cat("  Treatments: ", paste(unique(plot_cells$TREATMENT), collapse = ", "), "\n")
  cat("  Probe sets: ", paste(probe_sets, collapse = ", "), "\n")

  for (ps in probe_sets) {
    # Unavailable inference must not remove a requested contrast/follow-up page.
    cells <- plot_cells

    # One PDF per probe set, one page per stratum x FU x contrast.
    qq_file <- file.path(figures_dir, paste0("qq_", ps, ".pdf"))
    pdf(qq_file, width = 7, height = 6)
    tryCatch({
      for (i in seq_len(nrow(cells))) {
        plot_qq(results, stratum = cells$STRATUM[i], fu = cells$FU[i],
                probe_set = ps, treatment = cells$TREATMENT[i])
      }
    }, finally = {
      dev.off()
    })
    cat("  Saved:", qq_file, "\n")

    volcano_file <- file.path(figures_dir, paste0("volcano_", ps, ".pdf"))
    pdf(volcano_file, width = 8, height = 6)
    tryCatch({
      for (i in seq_len(nrow(cells))) {
        print(plot_volcano(results, stratum = cells$STRATUM[i], fu = cells$FU[i],
                           probe_set = ps, treatment = cells$TREATMENT[i]))
      }
    }, finally = {
      dev.off()
    })
    cat("  Saved:", volcano_file, "\n")
  }

  cat("Done.\n")
}
