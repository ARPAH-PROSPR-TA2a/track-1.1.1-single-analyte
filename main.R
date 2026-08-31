.log_111 <- function(verbose, text) {
  if (isTRUE(verbose)) {
    line <- paste0("[1.1.1] ", text)
    message(line)

    progress_log <- getOption("track111.progress_log")
    if (is.character(progress_log) && length(progress_log) == 1L &&
        !is.na(progress_log) && nzchar(progress_log)) {
      cat(
        "[1.1.1] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        " | ", text, "\n",
        file = progress_log, append = TRUE, sep = ""
      )
    }
  }
  invisible(NULL)
}

source("validation_helpers.R")
source("reporting_helpers.R")
source("analysis_helpers.R")

FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type = "Proteomics",
                           additional_covariates = NULL,
                           n_cores = NULL,
                           checkpoint_dir = NULL,
                           checkpoint_batch_size = 2000L,
                           verbose = FALSE) {

  analysis_started <- proc.time()[["elapsed"]]

  # Auto-detect cores if not specified, leaving one free for the OS.
  # Note: detectCores() may overcount in HPC/container environments — set
  # n_cores explicitly if running on a cluster with allocated core limits.
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }

  # future snapshots certain R options when a plan is set, including `exact`
  # (controls [[ partial matching). In a fresh session `exact` is unset (NULL).
  # future stores NULL, then tries to restore it via options(exact = NULL),
  # which fails in newer R because `exact` is a formal arg to base::options(),
  # not a plain option name. Set it to its documented default so future
  # captures a valid value instead of NULL.
  
  if ("exact" %in% names(.Options)) {
    .Options$exact <- NULL
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (n_cores > 1L) {
    future::plan(future::multicore, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

  .log_111(verbose, "analysis validating inputs and preparing data (serial)")
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

  backend <- if (n_cores > 1L) "multicore" else "sequential"
  checkpoint_status <- if (is.null(checkpoint_dir)) {
    "disabled"
  } else {
    paste0(checkpoint_dir, " (batch size=", checkpoint_batch_size, ")")
  }

  .log_111(
    verbose,
    paste0(
      "analysis starting: ", nrow(omics_list$all), " analytes, ",
      nrow(pheno_list$all), " phenotype rows; backend=", backend,
      ", effective workers=", future::nbrOfWorkers(),
      "; checkpoints=", checkpoint_status
    )
  )

  analysis_change <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "change", filtered_probes,
                                              checkpoint_dir, checkpoint_batch_size,
                                              verbose = verbose)
  analysis_level  <- .run_stratified_analysis(pheno_list, omics_list, omics_type,
                                              additional_covariates, "level", filtered_probes,
                                              checkpoint_dir, checkpoint_batch_size,
                                              verbose = verbose)

  .log_111(
    verbose,
    sprintf("analysis complete (%.1f sec)", proc.time()[["elapsed"]] - analysis_started)
  )

  return(list(
    analysis_change = analysis_change,
    analysis_level  = analysis_level
  ))
}

FAST_omics_WAS_reports <- function(pheno,
                                   omics,
                                   omics_type = "Proteomics",
                                   additional_covariates = NULL,
                                   verbose = FALSE) {

  reports_started <- proc.time()[["elapsed"]]
  .log_111(verbose, "reports validating inputs (serial)")

  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)

  if (omics_type == "DNAm") {
    full_probes     <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, filtered_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, full_probes)
  }

  reports <- .generate_reports(
    pheno_list,
    omics_list,
    additional_covariates,
    verbose = verbose
  )

  .log_111(
    verbose,
    sprintf("reports complete (%.1f sec)", proc.time()[["elapsed"]] - reports_started)
  )

  return(reports)
}
