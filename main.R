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

  filtered_probes <- NULL
  if (omics_type == "DNAm") {
    full_probes     <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, filtered_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, full_probes)
  }

  reports        <- .generate_reports(pheno_list, omics_list, additional_covariates)
  analysis_change <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "change", filtered_probes)
  analysis_level  <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "level", filtered_probes)

  return(list(
    analysis_change = analysis_change,
    analysis_level  = analysis_level,
    reports         = reports
  ))
}
