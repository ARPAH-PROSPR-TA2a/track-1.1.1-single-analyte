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
    # analyte ~ CONTROL_STATUS * factor(FU) + baseline_analyte + covariates + (1|SUBJECT_ID)
    # Extracts CONTROL_STATUS coefficient for FU=1, CONTROL_STATUS:factor(FU)2 for FU=2, etc.
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
        
        # Extract coefficients directly from summary
        fit_summary <- summary(fit)
        coef_table <- fit_summary$coefficients
        
        # Compute degrees of freedom for p-values
        # Using approximation: df = nrow(model_data) - number of fixed effects
        n_fixed_effects <- nrow(coef_table)
        df_approx <- nrow(model_data) - n_fixed_effects
        
        # Extract results per FU level
        for (fu_level in fu_levels) {
          if (fu_level == fu_levels[1]) {
            # First FU level: extract CONTROL_STATUS coefficient
            coef_name <- "CONTROL_STATUS"
          } else {
            # Subsequent FU levels: extract CONTROL_STATUS:factor(FU) interaction
            coef_name <- paste0("CONTROL_STATUS:factor(FU)", fu_level)
          }
          
          # Check if coefficient exists in model
          if (coef_name %in% rownames(coef_table)) {
            effect_size <- coef_table[coef_name, "Estimate"]
            se <- coef_table[coef_name, "Std. Error"]
            t_stat <- coef_table[coef_name, "t value"]
            # Compute p-value from t-statistic
            p_value <- 2 * pt(-abs(t_stat), df = df_approx)
            
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
   
    # Build design matrix for ALL samples with CONTROL_STATUS and FU interaction
    # Implicit baseline adjustment via change scores; repeated measures via duplicateCorrelation
    # Approximates LME4's (1|SUBJECT_ID) random intercept using LIMMA's block/correlation approach
    pheno_merged$FU_factor <- factor(pheno_merged$FU)
    design <- model.matrix(~ CONTROL_STATUS * FU_factor, data = pheno_merged)
   
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
   
    tryCatch({
      # Estimate within-subject correlation using duplicateCorrelation
      # This approximates random intercept structure (matches LME4's (1|SUBJECT_ID))
      cor <- duplicateCorrelation(analyte_change, design, block = pheno_merged$SUBJECT_ID)
      
      # Fit linear models for all analytes simultaneously using limma (vectorized)
      # Pass block and correlation to account for repeated measures within subjects
      fit <- lmFit(analyte_change, design, 
                   block = pheno_merged$SUBJECT_ID, 
                   correlation = cor$consensus.correlation)
      fit <- eBayes(fit)
     
     # Initialize results data frame with pre-allocated capacity (avoid rbind in loop)
     n_analytes <- nrow(analyte_change)
     n_fu_levels <- length(fu_levels)
     max_rows <- n_analytes * n_fu_levels
     
     results <- data.frame(
       ANALYTE_NAME = character(max_rows),
       EFFECT_SIZE = numeric(max_rows),
       SE = numeric(max_rows),
       P_VALUE = numeric(max_rows),
       FU = integer(max_rows),
       stringsAsFactors = FALSE
     )
     
     row_idx <- 0
     
      # Extract results per FU level
      # Each FU level gets its own direct coefficient from the model (no combining)
      for (i in seq_along(fu_levels)) {
        fu_level <- fu_levels[i]
        
        # Determine which coefficient to extract for this FU level
        if (fu_level == fu_levels[1]) {
          # First FU level: just CONTROL_STATUS coefficient
          coef_name <- "CONTROL_STATUS"
        } else {
          # Subsequent FU levels: extract CONTROL_STATUS:FU_factor interaction directly
          coef_name <- paste0("CONTROL_STATUS:FU_factor", fu_level)
        }
        
        # Find coefficient index
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
            results$EFFECT_SIZE[row_idx] <- effect_sizes[j]
            results$SE[row_idx] <- ses[j]
            results$P_VALUE[row_idx] <- p_values[j]
            results$FU[row_idx] <- fu_level
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
     
   }, error = function(e) {
     warning("Error fitting limma model: ", e$message)
     return(NULL)
   })
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