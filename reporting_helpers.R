.create_pheno_data_report <- function(pheno_df) {

  # One row per (FU, FEMALE) cell. Subject-level counts (N_SUBJECTS,
  # N_CONTROL, N_TREATMENT) dedupe by SUBJECT_ID so multiple samples for
  # the same person at the same FU don't inflate the totals. N_SAMPLES
  # is row-level so the difference between people present and samples
  # collected stays visible.
  groups <- unique(pheno_df[, c("FU", "FEMALE")])
  groups <- groups[order(groups$FU, groups$FEMALE), , drop = FALSE]

  report <- data.frame(
    FU          = groups$FU,
    FEMALE      = groups$FEMALE,
    N_SUBJECTS  = NA_integer_,
    N_CONTROL   = NA_integer_,
    N_TREATMENT = NA_integer_,
    N_SAMPLES   = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(report))) {
    cell <- pheno_df[pheno_df$FU == report$FU[i] &
                       pheno_df$FEMALE == report$FEMALE[i], ]
    report$N_SUBJECTS[i]  <- length(unique(cell$SUBJECT_ID))
    report$N_CONTROL[i]   <- length(unique(cell$SUBJECT_ID[cell$CONTROL_STATUS == 0]))
    report$N_TREATMENT[i] <- length(unique(cell$SUBJECT_ID[cell$CONTROL_STATUS == 1]))
    report$N_SAMPLES[i]   <- nrow(cell)
  }

  rownames(report) <- NULL
  report
}


.create_omics_data_report <- function(pheno_df, omics_df) {

  # Restrict to baseline (FU=0) samples so the report describes the
  # pre-treatment reference distribution for each analyte.
  baseline_sample_ids <- pheno_df$SAMPLE_ID[pheno_df$FU == 0]

  analyte_names <- omics_df$ANALYTE_NAME

  omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
  omics_numeric <- omics_numeric[, colnames(omics_numeric) %in% baseline_sample_ids, drop = FALSE]
  
  # Initialize results data.frame
  report <- data.frame(
    ANALYTE_NAME = analyte_names,
    N_NONMISSING = NA_integer_,
    MEAN = NA_real_,
    MEDIAN = NA_real_,
    SD = NA_real_,
    MIN = NA_real_,
    MAX = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Calculate per-analyte statistics
  for (i in seq_along(analyte_names)) {
    analyte_values <- as.numeric(omics_numeric[i, ])
    
    report$N_NONMISSING[i] <- sum(!is.na(analyte_values))
    report$MEAN[i] <- mean(analyte_values, na.rm = TRUE)
    report$MEDIAN[i] <- median(analyte_values, na.rm = TRUE)
    report$SD[i] <- sd(analyte_values, na.rm = TRUE)
    report$MIN[i] <- min(analyte_values, na.rm = TRUE)
    report$MAX[i] <- max(analyte_values, na.rm = TRUE)
  }
  
  return(report)
}


.create_addx_covariate_report <- function(pheno_df, covariate_names) {

   # Handle NULL or empty covariate_names
   if (is.null(covariate_names) || length(covariate_names) == 0) {
     return(NULL)
   }

   # Restrict to baseline (FU=0) samples so the report describes the
   # pre-treatment reference distribution for each covariate.
   pheno_df <- pheno_df[pheno_df$FU == 0, ]

   # Initialize results list to build then convert to data.frame
   results_list <- list(
     COVARIATE_NAME = character(),
     TYPE = character(),
     N_NA = integer(),
     SUMMARY = list()
   )
   
   # Process each covariate
   for (i in seq_along(covariate_names)) {
     covar_name <- covariate_names[i]
     covar_data <- pheno_df[[covar_name]]
     
     # Count NAs
     n_na <- sum(is.na(covar_data))
     
     # Determine type and generate summary
     if (is.numeric(covar_data)) {
       covar_type <- "numeric"
       summary_stats <- list(
         mean = mean(covar_data, na.rm = TRUE),
         median = median(covar_data, na.rm = TRUE),
         sd = sd(covar_data, na.rm = TRUE),
         min = min(covar_data, na.rm = TRUE),
         max = max(covar_data, na.rm = TRUE)
       )
     } else if (is.factor(covar_data)) {
       covar_type <- "factor"
       level_counts <- table(covar_data, useNA = "no")
       summary_stats <- list(
         n_levels = nlevels(covar_data),
         level_names = levels(covar_data),
         counts = as.numeric(level_counts)
       )
     } else if (is.logical(covar_data)) {
       covar_type <- "logical"
       summary_stats <- list(
         n_true = sum(covar_data == TRUE, na.rm = TRUE),
         n_false = sum(covar_data == FALSE, na.rm = TRUE)
       )
     } else {
       # Fallback for unexpected types
       covar_type <- "unknown"
       summary_stats <- list()
     }
     
     # Add to results list
     results_list$COVARIATE_NAME <- c(results_list$COVARIATE_NAME, covar_name)
     results_list$TYPE <- c(results_list$TYPE, covar_type)
     results_list$N_NA <- c(results_list$N_NA, n_na)
     results_list$SUMMARY[[i]] <- summary_stats
   }
   
   # Convert to data.frame
   report <- data.frame(
     COVARIATE_NAME = results_list$COVARIATE_NAME,
     TYPE = results_list$TYPE,
     N_NA = results_list$N_NA,
     stringsAsFactors = FALSE
   )
   
   # Add SUMMARY as a list column
   report$SUMMARY <- results_list$SUMMARY

   return(report)
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
