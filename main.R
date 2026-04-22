source("validation_helpers.R")
source("reporting_helpers.R")
source("analysis_helpers.R")

FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type = "Proteomics",
                           additional_covariates = NULL,
                           n_cores = NULL,
                           checkpoint_dir = NULL,
                           checkpoint_batch_size = 2000L) {

  # Auto-detect cores if not specified, leaving one free for the OS.
  # Note: detectCores() may overcount in HPC/container environments — set
  # n_cores explicitly if running on a cluster with allocated core limits.
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }

  old_plan <- future::plan("list")
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_cores > 1L) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

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

  analysis_change <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "change", filtered_probes,
                                              checkpoint_dir, checkpoint_batch_size)
  analysis_level  <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "level", filtered_probes,
                                              checkpoint_dir, checkpoint_batch_size)

  return(list(
    analysis_change = analysis_change,
    analysis_level  = analysis_level
  ))
}

FAST_omics_WAS_reports <- function(pheno,
                                   omics,
                                   omics_type = "Proteomics",
                                   additional_covariates = NULL) {

  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)

  if (omics_type == "DNAm") {
    full_probes     <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, filtered_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, full_probes)
  }

  return(.generate_reports(pheno_list, omics_list, additional_covariates))
}
