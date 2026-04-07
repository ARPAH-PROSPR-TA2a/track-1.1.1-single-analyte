source("validation_helpers.R")
source("reporting_helpers.R")
source("analysis_helpers.R")

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

# Helper function to subset omics data to a specific set of analytes
.subset_omics_list <- function(omics_list, analyte_subset) {
  if (is.null(analyte_subset)) {
    return(omics_list)
  }
  
  subsetted <- list()
  for (dataset in c("all", "male", "female")) {
    if (is.null(omics_list[[dataset]])) {
      subsetted[[dataset]] <- NULL
      next
    }
    
    omics_df <- omics_list[[dataset]]
    matching_analytes <- omics_df$ANALYTE_NAME %in% analyte_subset
    subsetted[[dataset]] <- omics_df[matching_analytes, ]
  }
  
  return(subsetted)
}

# Helper function to run analysis across all strata (all, male, female)
.run_stratified_analysis <- function(pheno_list, omics_list, omics_type,
                                     additional_covariates, analyte_subset = NULL) {
  
  # Subset omics data if analyte_subset provided
  omics_list <- .subset_omics_list(omics_list, analyte_subset)
  
  outputs <- list(all = NULL, male = NULL, female = NULL)
  
  for (dataset in c("all", "male", "female")) {
    
    if (is.null(pheno_list[[dataset]])) {
      next
    }
    
    # Generate data summary reports
    pheno_report <- .create_pheno_data_report(pheno_list[[dataset]])
    omics_report <- .create_omics_data_report(omics_list[[dataset]])
    
    if (!is.null(additional_covariates)) {
      covariates_report <- .create_addx_covariate_report(pheno_list[[dataset]], additional_covariates)
    } else {
      covariates_report <- NULL
    }
    
    # Run randomization analysis
    randomization_report <- .create_randomization_report(pheno_list[[dataset]], omics_list[[dataset]])
    
    # Perform omics-wide association analysis
    analysis_results <- .perform_analysis(
      pheno_list[[dataset]],
      omics_list[[dataset]],
      omics_type,
      pheno_list$requires_mixed_effects,
      additional_covariates
    )

    # analysis_results is now a list with coefficients and treatment_effects
    outputs[[dataset]] <- list(
      coefficients = analysis_results$coefficients,
      treatment_effects = analysis_results$treatment_effects,
      omics_summary = omics_report,
      pheno_summary = pheno_report,
      covariates_summary = covariates_report,
      randomization_summary = randomization_report
    )
  }
  
  return(outputs)
}

FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type = "Proteomics",
                           additional_covariates = NULL) {
  
  .validate_omics_type(omics_type)
  
  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)
  
  if (omics_type == "DNAm") {
    # Load probe lists for DNAm analysis
    full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")

    # Check probe coverage against available data
    available_probes <- omics_list$all$ANALYTE_NAME

    full_present <- sum(full_probes %in% available_probes)
    full_missing <- length(full_probes) - full_present
    if (full_present == 0) {
      stop("No probes from full probe list found in data")
    }
    if (full_missing > 0) {
      warning(sprintf("%d of %d probes from full probe list not found in data",
                      full_missing, length(full_probes)))
    }

    filtered_present <- sum(filtered_probes %in% available_probes)
    filtered_missing <- length(filtered_probes) - filtered_present
    if (filtered_present == 0) {
      stop("No probes from filtered probe list found in data")
    }
    if (filtered_missing > 0) {
      warning(sprintf("%d of %d probes from filtered probe list not found in data",
                      filtered_missing, length(filtered_probes)))
    }

    # Run analysis ONCE on full probes (LIMMA benefits from borrowing info across all probes)
    outputs <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                        additional_covariates, full_probes)

    # Add BH_P_VALUE_FILTERED column (BH correction on filtered probes only)
    outputs <- .add_filtered_bh_correction(outputs, filtered_probes)

  } else {
    # Proteomics/Metabolomics: single analysis with all analytes
    outputs <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                        additional_covariates)
  }
  
  return(outputs)
}
