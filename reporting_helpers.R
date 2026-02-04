.create_pheno_data_report <- function(pheno_df) {
  
  # Basic counts
  report <- list(
    N_SAMPLES = nrow(pheno_df),
    N_SUBJECTS = length(unique(pheno_df$SUBJECT_ID)),
    N_FEMALE = sum(pheno_df$FEMALE == 1),
    N_MALE = sum(pheno_df$FEMALE == 0),
    N_CONTROL = sum(pheno_df$CONTROL_STATUS == 0),
    N_TREATMENT = sum(pheno_df$CONTROL_STATUS == 1)
  )
  
  # Follow-up counts (only for FU levels present in data)
  fu_levels <- sort(unique(pheno_df$FU))
  for (fu in fu_levels) {
    report[[paste0("N_SAMPLES_FU", fu)]] <- sum(pheno_df$FU == fu)
  }
  
  return(report)
}


.create_omics_data_report <- function(omics_df) {
  
  analyte_names <- omics_df$ANALYTE_NAME
  
  omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
  
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
