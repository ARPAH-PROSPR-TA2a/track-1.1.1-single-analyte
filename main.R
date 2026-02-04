source("validation_helpers.R")
source("reporting_helpers.R")
source("analysis_helpers.R")

FAST_omics_WAS <- function(pheno, 
                           omics, 
                           omics_type = "Proteomics",
                           additional_covariates = NULL) {
  
  # Validate omics_type
  .validate_omics_type(omics_type)
  
  # Validate incoming data and create stratified datasets (all, male, female)
  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)
  
  # Initialize output structure
  outputs <- list(all = NULL, male = NULL, female = NULL)
  
  # Iterate over datasets
  for (dataset in c("all", "male", "female")) {
    
    # Skip if dataset is NULL
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
    
    # Store results for this dataset
    outputs[[dataset]] <- list(
      results = analysis_results,
      omics_summary = omics_report,
      pheno_summary = pheno_report,
      covariates_summary = covariates_report,
      randomization_summary = randomization_report
    )
  }
  
  return(outputs)
}