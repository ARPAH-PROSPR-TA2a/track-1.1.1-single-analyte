.apply_multiple_testing_correction <- function(results_df, p_value_col = "P_VALUE") {
  
  # Edge case: empty results
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }
  
  # Apply BH correction separately for each coefficient
  results_df$BH_P_VALUE <- NA_real_
  
  for (coef in unique(results_df$COEFFICIENT)) {
    coef_idx <- which(results_df$COEFFICIENT == coef)
    p_values_coef <- results_df[[p_value_col]][coef_idx]
    results_df$BH_P_VALUE[coef_idx] <- p.adjust(p_values_coef, method = "BH")
  }
  
  return(results_df)
}


.create_randomization_report <- function(pheno_df, omics_df) {
  
  # Filter to baseline (FU=0)
  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]
  
  # Filter omics to baseline samples
  baseline_sample_ids <- pheno_baseline$SAMPLE_ID
  omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]
  
  # Initialize results data.frame
  results <- data.frame(
    ANALYTE_NAME = character(),
    MEAN_DIFFERENCE = numeric(),
    COHENS_D = numeric(),
    SE = numeric(),
    P_VALUE = numeric(),
    stringsAsFactors = FALSE
  )
  
  # T-test per analyte
  analyte_names <- omics_df$ANALYTE_NAME
  
  for (i in seq_along(analyte_names)) {
    tryCatch({
      # Extract values for this analyte
      analyte_values <- as.numeric(omics_baseline[i, ])
      control_status <- pheno_baseline$CONTROL_STATUS
      
      # Split by group
      values_control <- analyte_values[control_status == 0]
      values_treatment <- analyte_values[control_status == 1]
      
      # Remove NAs
      values_control <- na.omit(values_control)
      values_treatment <- na.omit(values_treatment)
      
      # Skip if insufficient samples in either group
      if (length(values_control) < 2 || length(values_treatment) < 2) {
        return()
      }
      
      # Calculate statistics
      n_control <- length(values_control)
      n_treatment <- length(values_treatment)
      mean_control <- mean(values_control)
      mean_treatment <- mean(values_treatment)
      var_control <- var(values_control)
      var_treatment <- var(values_treatment)
      
      # Mean difference (Treatment - Control)
      mean_diff <- mean_treatment - mean_control
      
      # Standard Error (Welch's approach - doesn't assume equal variances)
      se <- sqrt((var_control / n_control) + (var_treatment / n_treatment))
      
      # Cohen's d (using pooled SD)
      pooled_sd <- sqrt(((n_control - 1) * var_control + (n_treatment - 1) * var_treatment) / 
                          (n_control + n_treatment - 2))
      cohens_d <- mean_diff / pooled_sd
      
      # Welch's t-test
      tt <- t.test(values_treatment, values_control, var.equal = FALSE)
      
      # Add to results
      results <- rbind(results, data.frame(
        ANALYTE_NAME = analyte_names[i],
        MEAN_DIFFERENCE = mean_diff,
        COHENS_D = cohens_d,
        SE = se,
        P_VALUE = tt$p.value,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }
  
   # Apply multiple testing correction (global across all analytes)
   if (nrow(results) > 0) {
     results$BH_P_VALUE <- p.adjust(results$P_VALUE, method = "BH")
     return(results)
  } else {
    return(NULL)
  }
}

.perform_lm_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL) {
   
    # Linear regression for single follow-up timepoint only
    # Extracts ALL fixed effect coefficients (treatment, covariates)
    
    # Initialize results
    results <- data.frame(
      ANALYTE_NAME = character(),
      COEFFICIENT = character(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )
   
   # Get sample IDs from omics (exclude ANALYTE_NAME column)
   omics_sample_ids <- colnames(omics_df)[-which(colnames(omics_df) == "ANALYTE_NAME")]
   
   # Filter pheno to shared samples
   shared_samples <- intersect(pheno_df$SAMPLE_ID, omics_sample_ids)
   pheno_merged <- pheno_df[pheno_df$SAMPLE_ID %in% shared_samples, ]

    # Baseline data for lookup
    baseline_subject_ids <- pheno_baseline$SUBJECT_ID
    
    # Convert omics_baseline to matrix for efficient numeric indexing
    # (avoids issues with duplicate column names when indexing with vectors)
    omics_baseline_matrix <- as.matrix(omics_baseline)
    
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
    sample_subjects <- pheno_merged$SUBJECT_ID
    baseline_idx <- match(sample_subjects, baseline_subject_ids)
    
    # Map baseline indices to column positions for omics_baseline_matrix indexing
    # Note: baseline_idx contains row positions in pheno_baseline;
    # we need to convert these to column positions in omics_baseline_matrix
    baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))
    
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
       
        # Compute change scores (on-the-fly)
        analyte_change <- fu_values - baseline_vals
        
        # Update model data with current analyte values
        model_data$analyte <- analyte_change
        model_data$analyte_baseline <- baseline_vals
        
         # Fit linear model
         fit <- lm(as.formula(formula_str), data = model_data)
         fit_summary <- summary(fit)
         
        # Extract all fixed effect coefficients (except intercept)
        coef_table <- fit_summary$coefficients
        
        # Loop through all coefficients
        for (coef_name in rownames(coef_table)) {
          # Extract coefficient info
          effect_size <- coef_table[coef_name, "Estimate"]
          se <- coef_table[coef_name, "Std. Error"]
          p_value <- coef_table[coef_name, "Pr(>|t|)"]
          
          # Add to results
          results <- rbind(results, data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT = coef_name,
            EFFECT_SIZE = effect_size,
            SE = se,
            P_VALUE = p_value,
            stringsAsFactors = FALSE
          ))
        }
      
    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }
  
  # Return NULL if no results
  if (nrow(results) == 0) {
    return(NULL)
  }
  
  return(results)
}


.perform_lme4_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL) {
   
    # Load required packages
    require(lme4)
   
    # Initialize results
    results <- data.frame(
      ANALYTE_NAME = character(),
      COEFFICIENT = character(),
      EFFECT_SIZE = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Get sample IDs from omics (exclude ANALYTE_NAME column)
   omics_sample_ids <- colnames(omics_df)[-which(colnames(omics_df) == "ANALYTE_NAME")]
   
   # Filter pheno to shared samples
   shared_samples <- intersect(pheno_df$SAMPLE_ID, omics_sample_ids)
   pheno_merged <- pheno_df[pheno_df$SAMPLE_ID %in% shared_samples, ]
   
   # Get FU levels present in data
   fu_levels <- sort(unique(pheno_merged$FU))
   
   # Baseline data for lookup
   baseline_subject_ids <- pheno_baseline$SUBJECT_ID
   
   # Convert omics_baseline to matrix for efficient numeric indexing
   # (avoids issues with duplicate column names when indexing with vectors)
   omics_baseline_matrix <- as.matrix(omics_baseline)
   
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
    sample_subjects <- pheno_merged$SUBJECT_ID
    baseline_idx <- match(sample_subjects, baseline_subject_ids)
    
    # Map baseline indices to column positions for omics_baseline_matrix indexing
    # Note: baseline_idx contains row positions in pheno_baseline;
    # we need to convert these to column positions in omics_baseline_matrix
    baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))
    
    # Fit model for each analyte
    for (i in seq_along(analyte_names)) {
      tryCatch({
        analyte_name <- analyte_names[i]
        
        # Get raw FU analyte values for shared samples (convert data.frame row to numeric vector)
        fu_values <- as.numeric(omics_df[i, shared_samples])
        
        # Get baseline analyte values using column position indexing
        # (avoids issues with duplicate column names)
        baseline_vals <- omics_baseline_matrix[i, baseline_col_idx]
       
        # Compute change scores (on-the-fly)
        analyte_change <- fu_values - baseline_vals
        
        # Update model data with current analyte values
        model_data$analyte <- analyte_change
        model_data$analyte_baseline <- baseline_vals
        
         # Fit lmer model
         fit <- lmer(as.formula(formula_str), data = model_data, REML = FALSE)
        
        # Extract coefficients directly from summary
        fit_summary <- summary(fit)
        coef_table <- fit_summary$coefficients
        
        # Compute degrees of freedom for p-values
        # Using approximation: df = nrow(model_data) - number of fixed effects
        n_fixed_effects <- nrow(coef_table)
        df_approx <- nrow(model_data) - n_fixed_effects
        
        # Extract all fixed effect coefficients
        for (coef_name in rownames(coef_table)) {
          # Extract coefficient info
          effect_size <- coef_table[coef_name, "Estimate"]
          se <- coef_table[coef_name, "Std. Error"]
          t_stat <- coef_table[coef_name, "t value"]
          # Compute p-value from t-statistic
          p_value <- 2 * pt(-abs(t_stat), df = df_approx)
          
          results <- rbind(results, data.frame(
            ANALYTE_NAME = analyte_name,
            COEFFICIENT = coef_name,
            EFFECT_SIZE = effect_size,
            SE = se,
            P_VALUE = p_value,
            stringsAsFactors = FALSE
          ))
        }
      
    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }
  
  # Return NULL if no results
  if (nrow(results) == 0) {
    return(NULL)
  }
  
  return(results)
}


.perform_limma_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL, requires_mixed_effects) {
   
   # Limma analysis for DNAm (DNA methylation) data
   # Vectorized approach: fits all FU levels simultaneously for speed with high-dimensional data
   # Uses empirical Bayes moderation for variance estimation (appropriate for 1M+ analytes)
   
   require(limma)
   
   # Get sample IDs from omics (exclude ANALYTE_NAME column)
   omics_sample_ids <- colnames(omics_df)[-which(colnames(omics_df) == "ANALYTE_NAME")]
   
   # Filter pheno to shared samples
   shared_samples <- intersect(pheno_df$SAMPLE_ID, omics_sample_ids)
   pheno_merged <- pheno_df[pheno_df$SAMPLE_ID %in% shared_samples, ]
   
   # Get FU levels present in data
   fu_levels <- sort(unique(pheno_merged$FU))
   
   # Baseline data for lookup
   baseline_subject_ids <- pheno_baseline$SUBJECT_ID
   
   # Convert omics_baseline to matrix for efficient numeric indexing
   omics_baseline_matrix <- as.matrix(omics_baseline)
   
   # Convert full omics data to matrix (analytes × samples)
   omics_values <- as.matrix(omics_df[, shared_samples])
   
   # Pre-compute baseline values for each sample using numeric indexing (O(1) vs O(S×N))
   sample_subjects <- pheno_merged$SUBJECT_ID
   baseline_idx <- match(sample_subjects, baseline_subject_ids)
   baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))
   
   # Get baseline values as matrix (analytes × samples)
   omics_baseline_merged <- omics_baseline_matrix[, baseline_col_idx, drop = FALSE]
   
    # Compute change scores vectorized: analyte_change is (n_analytes × n_samples) matrix
    analyte_change <- omics_values - omics_baseline_merged
    
     # Build design matrix based on FU structure
     if (requires_mixed_effects) {
       message("LIMMA: Multiple FU - using full model with FU effects")
       pheno_merged$FU_factor <- factor(pheno_merged$FU)
       design <- model.matrix(~ CONTROL_STATUS * FU_factor + FEMALE, data = pheno_merged)
     } else {
       message("LIMMA: Single FU - using simple model without FU effects")
       design <- model.matrix(~ CONTROL_STATUS + FEMALE, data = pheno_merged)
     }
     
     # Exclude FEMALE if it has only one level
     if ("FEMALE" %in% colnames(design) && length(unique(design[, "FEMALE"])) == 1) {
       design <- design[, colnames(design) != "FEMALE", drop = FALSE]
     }
     
     # Add additional covariates if provided
    if (!is.null(additional_covariates)) {
      for (cov in additional_covariates) {
        if (cov %in% colnames(pheno_merged)) {
          cov_vals <- pheno_merged[[cov]]
          if (!all(is.na(cov_vals))) {
            design <- cbind(design, cov_vals)
            colnames(design)[ncol(design)] <- cov
          }
        }
      }
    }
    
    # Fit model based on FU structure
    if (requires_mixed_effects) {
      # Multiple FU: estimate within-subject correlation
      cor <- duplicateCorrelation(analyte_change, design, block = pheno_merged$SUBJECT_ID)
      fit <- lmFit(analyte_change, design, 
                   block = pheno_merged$SUBJECT_ID, 
                   correlation = cor$consensus.correlation)
    } else {
      # Single FU: no repeated measures, no correlation estimation needed
      fit <- lmFit(analyte_change, design)
    }
    fit <- eBayes(fit)
    
    # Initialize results data frame with pre-allocated capacity (avoid rbind in loop)
    # Note: number of rows = n_analytes * number_of_coefficients
    n_analytes <- nrow(analyte_change)
    n_coefficients <- ncol(design)  # One column per coefficient (including intercept, which we'll skip)
    max_rows <- n_analytes * n_coefficients
    
    results <- data.frame(
      ANALYTE_NAME = character(max_rows),
      COEFFICIENT = character(max_rows),
      EFFECT_SIZE = numeric(max_rows),
      SE = numeric(max_rows),
      P_VALUE = numeric(max_rows),
      stringsAsFactors = FALSE
    )
    
    row_idx <- 0
    
    # Extract all fixed effect coefficients (except intercept)
    coef_names <- colnames(design)
    
    for (coef_name in coef_names) {
      # Find coefficient index in design matrix
      coef_idx <- which(colnames(design) == coef_name)
      
      if (length(coef_idx) > 0) {
        # Extract coefficients, SEs, and p-values (vectorized across analytes)
        effect_sizes <- fit$coefficients[, coef_idx]
        ses <- fit$stdev.unscaled[, coef_idx] * fit$sigma
        p_values <- fit$p.value[, coef_idx]
        
        # Add to results (vectorized append)
        for (j in seq_len(n_analytes)) {
          row_idx <- row_idx + 1
          results$ANALYTE_NAME[row_idx] <- omics_df$ANALYTE_NAME[j]
          results$COEFFICIENT[row_idx] <- coef_name
          results$EFFECT_SIZE[row_idx] <- effect_sizes[j]
          results$SE[row_idx] <- ses[j]
          results$P_VALUE[row_idx] <- p_values[j]
        }
      }
    }
    
    # Trim results to actual rows used
    results <- results[1:row_idx, ]
    
    # Return NULL if no results
    if (nrow(results) == 0) {
      return(NULL)
    }
    
    return(results)
}


.perform_analysis <- function(pheno_df, omics_df, omics_type, mixed_effects, additional_covariates = NULL) {
  
  # STEP 1: Extract baseline data for all analysis functions
  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]
  baseline_sample_ids <- pheno_baseline$SAMPLE_ID
  omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]
  
  # Filter to post-baseline analysis data (FU > 0)
  pheno_analysis <- pheno_df[pheno_df$FU != 0, ]
  omics_analysis <- omics_df
  
  # STEP 2: Dispatch to appropriate analysis function
  if (omics_type == "DNAm") {
    # Limma handles both single and multiple FU cases
    results <- .perform_limma_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates, mixed_effects)
  } else {
    # Proteomics/Metabolomics - check if LM or LME4
    max_fu <- max(as.numeric(as.character(pheno_analysis$FU)), na.rm = TRUE)
    if (max_fu == 1) {
      # Single follow-up: use linear regression
      results <- .perform_lm_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates)
    } else {
      # Multiple follow-ups: use mixed effects
      results <- .perform_lme4_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates)
    }
  }
  
   # STEP 3: Apply multiple testing correction
   if (!is.null(results) && nrow(results) > 0) {
     results <- .apply_multiple_testing_correction(results)
     
     # Sort by ANALYTE_NAME first, then COEFFICIENT for consistent ordering
     results <- results[order(results$ANALYTE_NAME, results$COEFFICIENT), ]
   }
   
   return(results)
}