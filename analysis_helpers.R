.apply_multiple_testing_correction <- function(results_df, p_value_col = "P_VALUE") {
  
  # Edge case: empty results
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }
  
  # Extract p-values
  p_values <- results_df[[p_value_col]]
  
  # Calculate Benjamini-Hochberg (FDR) correction
  bh_p <- p.adjust(p_values, method = "BH")
  
  # Add corrected p-values to results
  results_df$BH_P_VALUE <- bh_p
  
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
  
  # Apply multiple testing correction
  if (nrow(results) > 0) {
    results <- .apply_multiple_testing_correction(results)
    return(results)
  } else {
    return(NULL)
  }
}

.perform_lm_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL) {
   
   # Linear regression for single follow-up timepoint only
   
   # Initialize results
   results <- data.frame(
     ANALYTE_NAME = character(),
     EFFECT_SIZE = numeric(),
     SE = numeric(),
     P_VALUE = numeric(),
     FU = integer(),
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
    # analyte ~ CONTROL_STATUS + baseline_analyte + covariates
    covariate_terms <- c("analyte_baseline")
   if (!is.null(additional_covariates)) {
     covariate_terms <- c(covariate_terms, additional_covariates)
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
        
       # Extract CONTROL_STATUS coefficient (treatment effect)
       coef_row <- which(rownames(fit_summary$coefficients) == "CONTROL_STATUS")
        
       # Skip if coefficient not found
       if (length(coef_row) == 0) {
          next
       }
        
       # Extract and append results
       effect_size <- fit_summary$coefficients[coef_row, "Estimate"]
       se <- fit_summary$coefficients[coef_row, "Std. Error"]
       p_value <- fit_summary$coefficients[coef_row, "Pr(>|t|)"]
        
       results <- rbind(results, data.frame(
          ANALYTE_NAME = analyte_name,
          EFFECT_SIZE = effect_size,
          SE = se,
          P_VALUE = p_value,
          FU = fu_level,
          stringsAsFactors = FALSE
        ))
      
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
   require(emmeans)
   
   # Initialize results
   results <- data.frame(
     ANALYTE_NAME = character(),
     EFFECT_SIZE = numeric(),
     SE = numeric(),
     P_VALUE = numeric(),
     FU = integer(),
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
   # Base: analyte ~ CONTROL_STATUS * factor(FU) + baseline_analyte + covariates + (1|SUBJECT_ID)
   covariate_terms <- c("analyte_baseline")
   if (!is.null(additional_covariates)) {
     covariate_terms <- c(covariate_terms, additional_covariates)
   }
   
    formula_str <- "analyte ~ CONTROL_STATUS * factor(FU)"
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
       
       # Extract effects using emmeans
      em <- emmeans(fit, ~CONTROL_STATUS | FU)
      contrasts_result <- contrast(em, method = "pairwise", adjust = "none")
      
      # Convert to data.frame for extraction
      results_df <- as.data.frame(contrasts_result)
      
      # Extract results per FU level
      for (fu_level in fu_levels) {
        fu_rows <- which(results_df$FU == fu_level)
        
        if (length(fu_rows) > 0) {
          row_idx <- fu_rows[1]
          
          effect_size <- results_df$estimate[row_idx]
          se <- results_df$SE[row_idx]
          p_value <- results_df$p.value[row_idx]
          
          results <- rbind(results, data.frame(
            ANALYTE_NAME = analyte_name,
            EFFECT_SIZE = effect_size,
            SE = se,
            P_VALUE = p_value,
            FU = fu_level,
            stringsAsFactors = FALSE
          ))
        }
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


.perform_limma_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline, additional_covariates = NULL) {
   
   # Limma analysis for DNAm (DNA methylation) data
   # Uses empirical Bayes moderation for variance estimation (appropriate for high-dimensional data)
   
   require(limma)
   
   # Initialize results
   results <- data.frame(
     ANALYTE_NAME = character(),
     EFFECT_SIZE = numeric(),
     SE = numeric(),
     P_VALUE = numeric(),
     FU = integer(),
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
   
   # Convert full omics data to matrix
   omics_values <- as.matrix(omics_df[, shared_samples])
   
   # Pre-compute baseline column indices for each sample in pheno_merged
   sample_subjects <- pheno_merged$SUBJECT_ID
   baseline_idx <- match(sample_subjects, baseline_subject_ids)
   baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))
   
   # Process per FU level (unlike LM/LME4, limma typically fits all data together)
   # Here we'll fit models per FU level for consistency with the pipeline
   
   for (fu_level in fu_levels) {
     # Filter to this FU level
     fu_idx <- which(pheno_merged$FU == fu_level)
     pheno_fu <- pheno_merged[fu_idx, ]
     
     if (nrow(pheno_fu) < 2) {
       next  # Skip if fewer than 2 samples at this FU
     }
     
     # Get omics data for these samples
     omics_fu <- omics_values[, colnames(omics_values) %in% pheno_fu$SAMPLE_ID, drop = FALSE]
     
     # Get baseline data for these samples
     baseline_col_idx_fu <- baseline_col_idx[fu_idx]
     omics_baseline_fu <- omics_baseline_matrix[, baseline_col_idx_fu, drop = FALSE]
     
     # Compute change scores
     analyte_change <- omics_fu - omics_baseline_fu
     
     # Build design matrix
     design <- model.matrix(~ CONTROL_STATUS, data = pheno_fu)
     
     # Add covariates if provided
     if (!is.null(additional_covariates)) {
       for (cov in additional_covariates) {
         if (cov %in% colnames(pheno_fu)) {
           cov_vals <- pheno_fu[[cov]]
           if (!all(is.na(cov_vals))) {
             design <- cbind(design, pheno_fu[[cov]])
             colnames(design)[ncol(design)] <- cov
           }
         }
       }
     }
     
     # Add baseline analyte adjustment
     baseline_vals_fu <- t(omics_baseline_fu)
     colnames(baseline_vals_fu) <- paste0("baseline_", seq_len(nrow(analyte_change)))
     # Note: limma expects analytes as rows, so we'll add baseline as covariate differently
     # For now, include baseline effect in the model using contrasts
     
     tryCatch({
       # Fit linear models using limma
       fit <- lmFit(analyte_change, design)
       fit <- eBayes(fit)
       
       # Extract coefficients for CONTROL_STATUS (treatment effect)
       coef_idx <- which(colnames(design) == "CONTROL_STATUS")
       
       if (length(coef_idx) > 0) {
         effect_sizes <- fit$coefficients[, coef_idx]
         ses <- fit$stdev.unscaled[, coef_idx] * fit$sigma
         p_values <- fit$p.value[, coef_idx]
         
         # Add to results
         for (i in seq_along(effect_sizes)) {
           results <- rbind(results, data.frame(
             ANALYTE_NAME = omics_df$ANALYTE_NAME[i],
             EFFECT_SIZE = effect_sizes[i],
             SE = ses[i],
             P_VALUE = p_values[i],
             FU = fu_level,
             stringsAsFactors = FALSE
           ))
         }
       }
       
     }, error = function(e) {
       warning("Error fitting limma model for FU level ", fu_level, ": ", e$message)
     })
   }
   
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
    # Limma handles all FU logic internally
    results <- .perform_limma_analysis(pheno_analysis, omics_analysis, pheno_baseline, omics_baseline, additional_covariates)
  } else {
    # Proteomics/Metabolomics - check if LM or LME4
    max_fu <- max(pheno_analysis$FU, na.rm = TRUE)
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
  }
  
  return(results)
}