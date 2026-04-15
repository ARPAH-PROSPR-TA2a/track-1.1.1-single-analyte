source("validation_helpers.R")
source("reporting_helpers.R")
source("analysis_helpers.R")

FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type = "Proteomics",
                           additional_covariates = NULL) {

  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)

  outputs_change <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                             additional_covariates, "change")
  outputs_level <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                            additional_covariates, "level")

  return(list(change = outputs_change, level = outputs_level))
}
