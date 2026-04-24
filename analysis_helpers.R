.apply_multiple_testing_correction <- function(results_df, p_value_col = "P_VALUE",
                                                group_col = "COEFFICIENT") {

  # Edge case: empty results
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }

  # Apply BH correction separately for each group (coefficient or FU level)
  results_df$BH_P_VALUE <- NA_real_

  for (grp in unique(results_df[[group_col]])) {
    grp_idx <- which(results_df[[group_col]] == grp)
    p_values_grp <- results_df[[p_value_col]][grp_idx]
    results_df$BH_P_VALUE[grp_idx] <- p.adjust(p_values_grp, method = "BH")
  }

  return(results_df)
}


# Prepares data for analysis by matching samples between pheno and omics,
# and computing the mapping from FU samples to their baseline values.
# Used by .perform_lm_analysis() and .perform_lme4_analysis()
.prepare_analysis_data <- function(pheno_df, omics_df, pheno_baseline, omics_baseline) {

  # Get sample IDs from omics (exclude ANALYTE_NAME column)
  omics_sample_ids <- colnames(omics_df)[-which(colnames(omics_df) == "ANALYTE_NAME")]

  # Find shared samples between pheno and omics
  shared_samples <- intersect(pheno_df$SAMPLE_ID, omics_sample_ids)

  # Filter pheno to shared samples
  pheno_merged <- pheno_df[pheno_df$SAMPLE_ID %in% shared_samples, ]

  # Baseline subject IDs for lookup
  baseline_subject_ids <- pheno_baseline$SUBJECT_ID

  # Convert omics_baseline to matrix for efficient numeric indexing
  # (avoids issues with duplicate column names when indexing with vectors)
  omics_baseline_matrix <- as.matrix(omics_baseline)

  # Map each sample in pheno_merged to its baseline:
  # 1. Find which row in pheno_baseline has the same SUBJECT_ID
  # 2. Get that row's SAMPLE_ID
  # 3. Find which column in omics_baseline_matrix that corresponds to
  sample_subjects <- pheno_merged$SUBJECT_ID
  baseline_idx <- match(sample_subjects, baseline_subject_ids)
  baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))

  list(
    pheno_merged = pheno_merged,
    shared_samples = shared_samples,
    baseline_col_idx = baseline_col_idx,
    omics_baseline_matrix = omics_baseline_matrix
  )
}


.perform_lm_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL, response_type = c("change", "level"),
                                 checkpoint_dir = NULL, checkpoint_batch_size = 2000L) {

    # Linear regression for single follow-up timepoint only
    # Extracts ALL fixed effect coefficients (treatment, covariates)

    response_type <- match.arg(response_type)

    # Prepare data using shared helper
    prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
    pheno_merged <- prep$pheno_merged
    shared_samples <- prep$shared_samples
    baseline_col_idx <- prep$baseline_col_idx
    omics_baseline_matrix <- prep$omics_baseline_matrix

    # Model data template — each worker copies this and fills in analyte columns
    model_data <- data.frame(pheno_merged)
    model_data$analyte <- NA_real_
    model_data$analyte_baseline <- NA_real_

    # Build model formula
    # analyte ~ CONTROL_STATUS + FEMALE + baseline_analyte + covariates
    covariate_terms <- c("FEMALE", "analyte_baseline")
    if (!is.null(additional_covariates)) {
      covariate_terms <- c(covariate_terms, additional_covariates)
    }

    # Exclude FEMALE if it has only one level
    if (length(unique(model_data$FEMALE)) == 1) {
      covariate_terms <- setdiff(covariate_terms, "FEMALE")
    }

    formula_str <- "analyte ~ CONTROL_STATUS"
    if (length(covariate_terms) > 0) {
      formula_str <- paste(formula_str, paste(covariate_terms, collapse = " + "), sep = " + ")
    }

    analyte_names <- omics_df$ANALYTE_NAME

    # Get FU level for this analysis (should be single value in LM analysis)
    fu_level <- unique(pheno_merged$FU)
    if (length(fu_level) != 1) {
      warning("LM analysis expects single FU level, found:", length(fu_level))
      if (length(fu_level) > 1) {
        fu_level <- max(fu_level)
      } else {
        return(NULL)
      }
    }
    fu_level <- as.integer(as.character(fu_level))

    if (!is.null(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    }

    batches     <- split(seq_along(analyte_names),
                         ceiling(seq_along(analyte_names) / checkpoint_batch_size))
    all_results <- vector("list", length(analyte_names))

    for (b in seq_along(batches)) {
      batch      <- batches[[b]]
      batch_file <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, paste0("batch_", b, ".rds")) else NULL

      if (!is.null(batch_file) && file.exists(batch_file)) {
        all_results[batch] <- readRDS(batch_file)
        next
      }

      # Pre-extract per-analyte data so omics_df and omics_baseline_matrix are
      # not captured in the future_map closure. Without this, future exports the
      # full matrices (~GB for DNAm) to every worker rather than one row each.
      batch_items <- lapply(batch, function(i) list(
        analyte_name  = analyte_names[i],
        fu_values     = as.numeric(omics_df[i, shared_samples]),
        baseline_vals = omics_baseline_matrix[i, baseline_col_idx]
      ))

      batch_results <- furrr::future_map(batch_items, function(item) {
        tryCatch({
          analyte_name  <- item$analyte_name
          fu_values     <- item$fu_values
          baseline_vals <- item$baseline_vals

          md <- model_data
          if (response_type == "change") {
            md$analyte <- fu_values - baseline_vals
          } else {
            md$analyte <- fu_values
          }
          md$analyte_baseline <- baseline_vals

          fit        <- lm(as.formula(formula_str), data = md)
          n_obs      <- nrow(fit$model)
          coef_table <- summary(fit)$coefficients

          coefs <- data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT  = rownames(coef_table),
            N_OBS        = n_obs,
            EFFECT_SIZE  = coef_table[, "Estimate"],
            SE           = coef_table[, "Std. Error"],
            P_VALUE      = coef_table[, "Pr(>|t|)"],
            stringsAsFactors = FALSE,
            row.names = NULL
          )

          ctrl_idx <- grepl("^CONTROL_STATUS", rownames(coef_table))
          te <- data.frame(
            ANALYTE_NAME = analyte_name,
            FU           = fu_level,
            EFFECT_SIZE  = coef_table[ctrl_idx, "Estimate"],
            SE           = coef_table[ctrl_idx, "Std. Error"],
            P_VALUE      = coef_table[ctrl_idx, "Pr(>|t|)"],
            stringsAsFactors = FALSE,
            row.names = NULL
          )

          list(coefficients = coefs, treatment_effects = te)

        }, error = function(e) {
          warning("Error processing analyte '", item$analyte_name, "': ", e$message)
          NULL
        })
      }, .options = furrr::furrr_options(seed = TRUE))

      if (!is.null(batch_file)) {
        saveRDS(batch_results, paste0(batch_file, ".tmp"))
        file.rename(paste0(batch_file, ".tmp"), batch_file)
      }

      all_results[batch] <- batch_results
    }

    all_results <- Filter(Negate(is.null), all_results)

    if (length(all_results) == 0) return(NULL)

    coefficients      <- do.call(rbind, lapply(all_results, `[[`, "coefficients"))
    treatment_effects <- do.call(rbind, lapply(all_results, `[[`, "treatment_effects"))
    row.names(coefficients)      <- NULL
    row.names(treatment_effects) <- NULL

    return(list(
      coefficients = coefficients,
      treatment_effects = treatment_effects
    ))
}


.perform_lme4_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL, response_type = c("change", "level"),
                                   checkpoint_dir = NULL, checkpoint_batch_size = 2000L) {

    require(lme4)
    require(lmerTest)
    require(emmeans)

    response_type <- match.arg(response_type)

    # Prepare data using shared helper
    prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
    pheno_merged <- prep$pheno_merged
    shared_samples <- prep$shared_samples
    baseline_col_idx <- prep$baseline_col_idx
    omics_baseline_matrix <- prep$omics_baseline_matrix

    # Model data template — each worker copies this and fills in analyte columns
    model_data <- data.frame(pheno_merged)
    model_data$analyte <- NA_real_
    model_data$analyte_baseline <- NA_real_

    # Build model formula
    # analyte ~ CONTROL_STATUS * factor(FU) + FEMALE + baseline_analyte + covariates + (1|SUBJECT_ID)
    # Extracts ALL fixed effect coefficients from the model
    covariate_terms <- c("FEMALE", "analyte_baseline")
    if (!is.null(additional_covariates)) {
      covariate_terms <- c(covariate_terms, additional_covariates)
    }

    # Exclude FEMALE if it has only one level
    if (length(unique(model_data$FEMALE)) == 1) {
      covariate_terms <- setdiff(covariate_terms, "FEMALE")
    }

    formula_str <- "analyte ~ CONTROL_STATUS * FU"
    if (length(covariate_terms) > 0) {
      formula_str <- paste(formula_str, paste(covariate_terms, collapse = " + "), sep = " + ")
    }
    formula_str <- paste(formula_str, "+ (1|SUBJECT_ID)")

    analyte_names <- omics_df$ANALYTE_NAME

    if (!is.null(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    }

    batches     <- split(seq_along(analyte_names),
                         ceiling(seq_along(analyte_names) / checkpoint_batch_size))
    all_results <- vector("list", length(analyte_names))

    for (b in seq_along(batches)) {
      batch      <- batches[[b]]
      batch_file <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, paste0("batch_", b, ".rds")) else NULL

      if (!is.null(batch_file) && file.exists(batch_file)) {
        all_results[batch] <- readRDS(batch_file)
        next
      }

      # Pre-extract per-analyte data so omics_df and omics_baseline_matrix are
      # not captured in the future_map closure. Without this, future exports the
      # full matrices (~GB for DNAm) to every worker rather than one row each.
      batch_items <- lapply(batch, function(i) list(
        analyte_name  = analyte_names[i],
        fu_values     = as.numeric(omics_df[i, shared_samples]),
        baseline_vals = omics_baseline_matrix[i, baseline_col_idx]
      ))

      batch_results <- furrr::future_map(batch_items, function(item) {
        tryCatch({
          analyte_name  <- item$analyte_name
          fu_values     <- item$fu_values
          baseline_vals <- item$baseline_vals

          md <- model_data
          if (response_type == "change") {
            md$analyte <- fu_values - baseline_vals
          } else {
            md$analyte <- fu_values
          }
          md$analyte_baseline <- baseline_vals

          fit        <- lmerTest::lmer(as.formula(formula_str), data = md, REML = FALSE,
                                      control = lme4::lmerControl(calc.derivs = FALSE))
          n_obs      <- nobs(fit)
          coef_table <- summary(fit)$coefficients

          coefs <- data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT  = rownames(coef_table),
            N_OBS        = n_obs,
            EFFECT_SIZE  = coef_table[, "Estimate"],
            SE           = coef_table[, "Std. Error"],
            P_VALUE      = coef_table[, "Pr(>|t|)"],
            stringsAsFactors = FALSE,
            row.names = NULL
          )

          emm      <- emmeans::emmeans(fit, ~ CONTROL_STATUS | FU)
          contr_df <- as.data.frame(pairs(emm, reverse = TRUE))

          te <- data.frame(
            ANALYTE_NAME = analyte_name,
            FU           = as.integer(as.character(contr_df$FU)),
            EFFECT_SIZE  = contr_df$estimate,
            SE           = contr_df$SE,
            P_VALUE      = contr_df$p.value,
            stringsAsFactors = FALSE,
            row.names = NULL
          )

          list(coefficients = coefs, treatment_effects = te)

        }, error = function(e) {
          warning("Error processing analyte '", item$analyte_name, "': ", e$message)
          NULL
        })
      }, .options = furrr::furrr_options(seed = TRUE, packages = c("lme4", "lmerTest", "emmeans")))

      if (!is.null(batch_file)) {
        saveRDS(batch_results, paste0(batch_file, ".tmp"))
        file.rename(paste0(batch_file, ".tmp"), batch_file)
      }

      all_results[batch] <- batch_results
    }

    all_results <- Filter(Negate(is.null), all_results)

    if (length(all_results) == 0) return(NULL)

    coefficients      <- do.call(rbind, lapply(all_results, `[[`, "coefficients"))
    treatment_effects <- do.call(rbind, lapply(all_results, `[[`, "treatment_effects"))
    row.names(coefficients)      <- NULL
    row.names(treatment_effects) <- NULL

    return(list(
      coefficients = coefficients,
      treatment_effects = treatment_effects
    ))
}


.perform_analysis <- function(pheno_df, omics_df, omics_type, mixed_effects, additional_covariates = NULL, response_type = c("change", "level"),
                              checkpoint_dir = NULL, checkpoint_batch_size = 2000L) {

  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]
  baseline_sample_ids <- pheno_baseline$SAMPLE_ID
  omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]

  pheno_analysis <- pheno_df[pheno_df$FU != 0, ]
  omics_analysis <- omics_df

  max_fu <- max(as.numeric(as.character(pheno_analysis$FU)), na.rm = TRUE)
  if (max_fu == 1) {
    results <- .perform_lm_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates, response_type, checkpoint_dir, checkpoint_batch_size)
  } else {
    results <- .perform_lme4_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates, response_type, checkpoint_dir, checkpoint_batch_size)
  }

  if (!is.null(results)) {
    if (!is.null(results$coefficients) && nrow(results$coefficients) > 0) {
      results$coefficients <- .apply_multiple_testing_correction(
        results$coefficients, group_col = "COEFFICIENT"
      )
      results$coefficients <- results$coefficients[
        order(results$coefficients$ANALYTE_NAME, results$coefficients$COEFFICIENT),
      ]
    }

    if (!is.null(results$treatment_effects) && nrow(results$treatment_effects) > 0) {
      results$treatment_effects <- .apply_multiple_testing_correction(
        results$treatment_effects, group_col = "FU"
      )
      results$treatment_effects <- results$treatment_effects[
        order(results$treatment_effects$ANALYTE_NAME, results$treatment_effects$FU),
      ]
    }
  }

  return(results)
}


# Helper function to add BH_P_VALUE_FILTERED column to a results data frame
# Applies BH correction only to probes in filtered_probes, NA for others
.add_filtered_bh_column <- function(df, filtered_probes, group_col) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  df$BH_P_VALUE_FILTERED <- NA_real_

  for (grp in unique(df[[group_col]])) {
    # Subset to this group AND filtered probes
    idx <- which(df[[group_col]] == grp & df$ANALYTE_NAME %in% filtered_probes)
    if (length(idx) > 0) {
      df$BH_P_VALUE_FILTERED[idx] <- p.adjust(df$P_VALUE[idx], method = "BH")
    }
  }

  return(df)
}


# Helper function to add filtered BH correction to all strata
.add_filtered_bh_correction <- function(outputs, filtered_probes) {
  for (stratum in c("all", "male", "female")) {
    if (is.null(outputs[[stratum]])) next

    # Add column to coefficients
    outputs[[stratum]]$coefficients <- .add_filtered_bh_column(
      outputs[[stratum]]$coefficients, filtered_probes, group_col = "COEFFICIENT"
    )

    # Add column to treatment_effects
    outputs[[stratum]]$treatment_effects <- .add_filtered_bh_column(
      outputs[[stratum]]$treatment_effects, filtered_probes, group_col = "FU"
    )
  }
  return(outputs)
}


.run_stratified_analysis <- function(pheno_list, omics_list, omics_type,
                                     additional_covariates, response_type = c("change", "level"),
                                     filtered_probes = NULL,
                                     checkpoint_dir = NULL, checkpoint_batch_size = 2000L) {

  response_type <- match.arg(response_type)

  outputs <- list(all = NULL, male = NULL, female = NULL)

  for (dataset in c("all", "male", "female")) {

    if (is.null(pheno_list[[dataset]])) next

    stratum_checkpoint_dir <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, response_type, dataset)
    } else NULL

    analysis_results <- .perform_analysis(
      pheno_list[[dataset]],
      omics_list[[dataset]],
      omics_type,
      pheno_list$requires_mixed_effects,
      additional_covariates,
      response_type,
      stratum_checkpoint_dir,
      checkpoint_batch_size
    )

    outputs[[dataset]] <- list(
      coefficients      = analysis_results$coefficients,
      treatment_effects = analysis_results$treatment_effects
    )
  }

  if (!is.null(filtered_probes)) {
    outputs <- .add_filtered_bh_correction(outputs, filtered_probes)
  }

  return(outputs)
}