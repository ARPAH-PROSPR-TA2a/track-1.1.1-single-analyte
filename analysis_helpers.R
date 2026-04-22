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


.perform_lm_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL, response_type = c("change", "level")) {

    # Linear regression for single follow-up timepoint only
    # Extracts ALL fixed effect coefficients (treatment, covariates)

    # Initialize results - raw coefficients
    coefficients <- data.frame(
      ANALYTE_NAME = character(),
      COEFFICIENT = character(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )

    # Initialize results - treatment effects (single FU for LM)
    treatment_effects <- data.frame(
      ANALYTE_NAME = character(),
      FU = integer(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )

    # Prepare data using shared helper
    prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
    pheno_merged <- prep$pheno_merged
    shared_samples <- prep$shared_samples
    baseline_col_idx <- prep$baseline_col_idx
    omics_baseline_matrix <- prep$omics_baseline_matrix

    # Pre-allocate model data template (avoid copying in loop)
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

    # Pre-compute loop invariants
    analyte_names <- omics_df$ANALYTE_NAME

    # Get FU level for this analysis (should be single value in LM analysis)
    fu_level <- unique(pheno_merged$FU)
    if (length(fu_level) != 1) {
      warning("LM analysis expects single FU level, found:", length(fu_level))
      if (length(fu_level) > 1) {
        fu_level <- max(fu_level)  # Use max FU if multiple present
      } else {
        return(NULL)  # No data to analyze
      }
    }
    
    # Fit model for each analyte
    for (i in seq_along(analyte_names)) {
      tryCatch({
        analyte_name <- analyte_names[i]
        
        # Get raw FU analyte values for shared samples (convert data.frame row to numeric vector)
        fu_values <- as.numeric(omics_df[i, shared_samples])
        
        # Get baseline analyte values using column position indexing
        # (avoids issues with duplicate column names)
        baseline_vals <- omics_baseline_matrix[i, baseline_col_idx]

        # Update model data with current analyte values
        # response_type determines whether we model change or absolute level
        if (match.arg(response_type) == "change") {
          model_data$analyte <- fu_values - baseline_vals
        } else {
          model_data$analyte <- fu_values
        }
        model_data$analyte_baseline <- baseline_vals

         # Fit linear model
         fit <- lm(as.formula(formula_str), data = model_data)
         fit_summary <- summary(fit)
         
        # Extract all fixed effect coefficients
        coef_table <- fit_summary$coefficients

        # Loop through all coefficients
        for (coef_name in rownames(coef_table)) {
          # Extract coefficient info
          effect_size <- coef_table[coef_name, "Estimate"]
          se <- coef_table[coef_name, "Std. Error"]
          p_value <- coef_table[coef_name, "Pr(>|t|)"]

          # Add to coefficients
          coefficients <- rbind(coefficients, data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT = coef_name,
            EFFECT_SIZE = effect_size,
            SE = se,
            P_VALUE = p_value,
            stringsAsFactors = FALSE
          ))

          # Extract treatment effect (CONTROL_STATUS coefficient)
          if (grepl("^CONTROL_STATUS", coef_name)) {
            treatment_effects <- rbind(treatment_effects, data.frame(
              ANALYTE_NAME = analyte_name,
              FU = as.integer(as.character(fu_level)),
              EFFECT_SIZE = effect_size,
              SE = se,
              P_VALUE = p_value,
              stringsAsFactors = FALSE
            ))
          }
        }

    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }

  # Return NULL if no results
  if (nrow(coefficients) == 0) {
    return(NULL)
  }

  return(list(
    coefficients = coefficients,
    treatment_effects = treatment_effects
  ))
}


.perform_lme4_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL, response_type = c("change", "level")) {

    # Load required packages
    require(lme4)
    require(lmerTest)
    require(emmeans)

    # Initialize results - raw coefficients
    coefficients <- data.frame(
      ANALYTE_NAME = character(),
      COEFFICIENT = character(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )

    # Initialize results - treatment effects at each FU
    treatment_effects <- data.frame(
      ANALYTE_NAME = character(),
      FU = integer(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )

    # Prepare data using shared helper
    prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
    pheno_merged <- prep$pheno_merged
    shared_samples <- prep$shared_samples
    baseline_col_idx <- prep$baseline_col_idx
    omics_baseline_matrix <- prep$omics_baseline_matrix

    # Pre-allocate model data template (avoid copying in loop)
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

    # Pre-compute loop invariants
    analyte_names <- omics_df$ANALYTE_NAME

    # Fit model for each analyte
    for (i in seq_along(analyte_names)) {
      tryCatch({
        analyte_name <- analyte_names[i]
        
        # Get raw FU analyte values for shared samples (convert data.frame row to numeric vector)
        fu_values <- as.numeric(omics_df[i, shared_samples])
        
        # Get baseline analyte values using column position indexing
        # (avoids issues with duplicate column names)
        baseline_vals <- omics_baseline_matrix[i, baseline_col_idx]

        # Update model data with current analyte values
        # response_type determines whether we model change or absolute level
        if (match.arg(response_type) == "change") {
          model_data$analyte <- fu_values - baseline_vals
        } else {
          model_data$analyte <- fu_values
        }
        model_data$analyte_baseline <- baseline_vals

         # Fit lmer model
         fit <- lmerTest::lmer(as.formula(formula_str), data = model_data, REML = FALSE,
                               control = lmerControl(calc.derivs = FALSE))

        # Extract coefficients with Satterthwaite p-values
        fit_summary <- summary(fit)
        coef_table <- fit_summary$coefficients

        # Extract all fixed effect coefficients
        for (coef_name in rownames(coef_table)) {
          effect_size <- coef_table[coef_name, "Estimate"]
          se <- coef_table[coef_name, "Std. Error"]
          p_value <- coef_table[coef_name, "Pr(>|t|)"]

          coefficients <- rbind(coefficients, data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT = coef_name,
            EFFECT_SIZE = effect_size,
            SE = se,
            P_VALUE = p_value,
            stringsAsFactors = FALSE
          ))
        }

        # Extract treatment effects at each FU using emmeans
        emm <- emmeans(fit, ~ CONTROL_STATUS | FU)
        contr <- pairs(emm, reverse = TRUE)  # (treatment - control)
        contr_df <- as.data.frame(contr)

        for (j in seq_len(nrow(contr_df))) {
          treatment_effects <- rbind(treatment_effects, data.frame(
            ANALYTE_NAME = analyte_name,
            FU = as.integer(as.character(contr_df$FU[j])),
            EFFECT_SIZE = contr_df$estimate[j],
            SE = contr_df$SE[j],
            P_VALUE = contr_df$p.value[j],
            stringsAsFactors = FALSE
          ))
        }

    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }

  # Return NULL if no results
  if (nrow(coefficients) == 0) {
    return(NULL)
  }

  return(list(
    coefficients = coefficients,
    treatment_effects = treatment_effects
  ))
}


.perform_analysis <- function(pheno_df, omics_df, omics_type, mixed_effects, additional_covariates = NULL, response_type = c("change", "level")) {

  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]
  baseline_sample_ids <- pheno_baseline$SAMPLE_ID
  omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]

  pheno_analysis <- pheno_df[pheno_df$FU != 0, ]
  omics_analysis <- omics_df

  max_fu <- max(as.numeric(as.character(pheno_analysis$FU)), na.rm = TRUE)
  if (max_fu == 1) {
    results <- .perform_lm_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates, response_type)
  } else {
    results <- .perform_lme4_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates, response_type)
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
                                     filtered_probes = NULL) {

  response_type <- match.arg(response_type)

  outputs <- list(all = NULL, male = NULL, female = NULL)

  for (dataset in c("all", "male", "female")) {

    if (is.null(pheno_list[[dataset]])) next

    analysis_results <- .perform_analysis(
      pheno_list[[dataset]],
      omics_list[[dataset]],
      omics_type,
      pheno_list$requires_mixed_effects,
      additional_covariates,
      response_type
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