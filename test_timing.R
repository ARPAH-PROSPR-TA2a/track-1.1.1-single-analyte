# ===== TIMING TEST FOR DNAm PIPELINES =====
#
# Runs the DNAm pipeline once with LM (single FU) and once with LME4 (multiple FU),
# timing each internal step to identify performance bottlenecks.
#
# Steps timed:
#   [Setup]   Data loading (readRDS)
#   [Setup]   Data preparation / subsetting
#   [Both]    .validate_pheno / .validate_omics / DNAm probe validation
#   [Both]    .generate_reports
#   [Both]    .run_stratified_analysis  (change model)
#   [Both]    .run_stratified_analysis  (level model)

source("main.R")
require(lme4)

# ===== TIMING HELPERS =====

.timings <- list()

#' Time an expression and store the elapsed seconds
#' @param label Character label for this step
#' @param expr Expression to evaluate
#' @return Result of the expression (invisibly)
time_step <- function(label, expr) {
  cat("  Running: ", label, "...\n", sep = "")
  t <- system.time(result <- expr)
  elapsed <- unname(t["elapsed"])
  .timings[[label]] <<- elapsed
  cat("    -> ", round(elapsed, 2), "s\n", sep = "")
  invisible(result)
}

#' Print a formatted timing summary table
print_timing_summary <- function() {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("TIMING SUMMARY\n")
  cat(strrep("=", 60), "\n")

  labels  <- names(.timings)
  seconds <- unlist(.timings)
  total   <- sum(seconds)

  max_label <- max(nchar(labels))
  fmt <- paste0("  %-", max_label, "s  %7.2fs  (%5.1f%%)\n")

  for (i in seq_along(labels)) {
    cat(sprintf(fmt, labels[i], seconds[i], 100 * seconds[i] / total))
  }

  cat(strrep("-", 60), "\n")
  cat(sprintf(paste0("  %-", max_label, "s  %7.2fs  (100.0%%)\n"), "TOTAL", total))
  cat(strrep("=", 60), "\n")
}


# =============================================================================
# DATA LOADING
# =============================================================================

cat("\n[SETUP] Loading data\n")

pheno_raw <- time_step("load pheno_raw", readRDS("PracticeData/pheno_example.rds"))
omics_raw <- time_step("load omics_raw", readRDS("PracticeData/synth_small_betas.rds"))


# =============================================================================
# DATA PREPARATION
# =============================================================================

cat("\n[SETUP] Preparing data\n")

prep_result <- time_step("prepare pheno/omics", {
  # Convert haven_labelled columns if present
  if (any(sapply(pheno_raw, function(x) inherits(x, "haven_labelled")))) {
    for (col in names(pheno_raw)) {
      if (inherits(pheno_raw[[col]], "haven_labelled")) {
        raw_values <- as.vector(pheno_raw[[col]])
        numeric_attempt <- tryCatch(as.numeric(raw_values), error = function(e) NA)
        if (!all(is.na(numeric_attempt))) {
          pheno_raw[[col]] <- numeric_attempt
        } else {
          pheno_raw[[col]] <- as.character(raw_values)
        }
      }
    }
  }

  pheno <- pheno_raw[!duplicated(pheno_raw$SAMPLE_ID), ]
  pheno <- pheno[complete.cases(pheno[, c("SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "CONTROL_STATUS")]), ]

  analyte_names <- rownames(omics_raw)
  sample_names  <- colnames(omics_raw)
  omics_full    <- as.data.frame(omics_raw)
  colnames(omics_full) <- sample_names
  omics_full <- cbind(ANALYTE_NAME = analyte_names, omics_full, stringsAsFactors = FALSE)

  full_probes     <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
  filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")

  in_full     <- omics_full$ANALYTE_NAME %in% full_probes
  in_filtered <- omics_full$ANALYTE_NAME %in% filtered_probes

  filtered_available <- omics_full[in_filtered, ]
  full_only          <- omics_full[in_full & !in_filtered, ]
  n_filtered_use     <- min(100, nrow(filtered_available))
  n_full_only_use    <- min(400, nrow(full_only))
  omics <- rbind(filtered_available[seq_len(n_filtered_use), ],
                 full_only[seq_len(n_full_only_use), ])

  pheno_single_fu <- pheno[pheno$FU %in% c(0, 1), ]
  pheno_single_fu <- pheno_single_fu[!duplicated(pheno_single_fu$SAMPLE_ID), ]
  pheno_multi_fu  <- pheno

  list(pheno = pheno,
       omics = omics,
       pheno_single_fu = pheno_single_fu,
       pheno_multi_fu  = pheno_multi_fu,
       full_probes     = full_probes,
       filtered_probes = filtered_probes)
})

pheno           <- prep_result$pheno
omics           <- prep_result$omics
pheno_single_fu <- prep_result$pheno_single_fu
pheno_multi_fu  <- prep_result$pheno_multi_fu
full_probes     <- prep_result$full_probes
filtered_probes <- prep_result$filtered_probes
additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")

cat("\n  Data summary:\n")
cat("    Pheno (single FU): ", nrow(pheno_single_fu), " samples, FU: ",
    paste(sort(unique(pheno_single_fu$FU)), collapse = ", "), "\n", sep = "")
cat("    Pheno (multi  FU): ", nrow(pheno_multi_fu),  " samples, FU: ",
    paste(sort(unique(pheno_multi_fu$FU)),  collapse = ", "), "\n", sep = "")
cat("    Omics:             ", nrow(omics), " analytes\n", sep = "")


# =============================================================================
# DNAm + SINGLE FU  (LM)
# =============================================================================

cat("\n[LM] DNAm + Single FU\n")

lm_pheno_list <- time_step("LM: validate_pheno",
  .validate_pheno(pheno_single_fu, additional_covariates))

lm_omics_list_raw <- time_step("LM: validate_omics",
  .validate_omics(omics, lm_pheno_list))

lm_omics_list <- time_step("LM: DNAm probe validation + subset", {
  .validate_dnam_probe_coverage(full_probes, filtered_probes, lm_omics_list_raw$all$ANALYTE_NAME)
  .subset_omics_list(lm_omics_list_raw, full_probes)
})

lm_reports <- time_step("LM: generate_reports",
  .generate_reports(lm_pheno_list, lm_omics_list, additional_covariates))

lm_change <- time_step("LM: run_stratified_analysis (change)",
  .run_stratified_analysis(lm_pheno_list, lm_omics_list, "DNAm",
                           additional_covariates, "change", filtered_probes))

lm_level <- time_step("LM: run_stratified_analysis (level)",
  .run_stratified_analysis(lm_pheno_list, lm_omics_list, "DNAm",
                           additional_covariates, "level", filtered_probes))

cat("  LM change coefficients: ", nrow(lm_change$all$coefficients), " rows, ",
    length(unique(lm_change$all$coefficients$ANALYTE_NAME)), " analytes\n", sep = "")
cat("  LM level  coefficients: ", nrow(lm_level$all$coefficients),  " rows, ",
    length(unique(lm_level$all$coefficients$ANALYTE_NAME)),  " analytes\n", sep = "")


# =============================================================================
# DNAm + MULTIPLE FU  (LME4)
# =============================================================================

cat("\n[LME4] DNAm + Multiple FU\n")

lme4_pheno_list <- time_step("LME4: validate_pheno",
  .validate_pheno(pheno_multi_fu, additional_covariates))

lme4_omics_list_raw <- time_step("LME4: validate_omics",
  .validate_omics(omics, lme4_pheno_list))

lme4_omics_list <- time_step("LME4: DNAm probe validation + subset", {
  .validate_dnam_probe_coverage(full_probes, filtered_probes, lme4_omics_list_raw$all$ANALYTE_NAME)
  .subset_omics_list(lme4_omics_list_raw, full_probes)
})

lme4_reports <- time_step("LME4: generate_reports",
  .generate_reports(lme4_pheno_list, lme4_omics_list, additional_covariates))

lme4_change <- time_step("LME4: run_stratified_analysis (change)",
  .run_stratified_analysis(lme4_pheno_list, lme4_omics_list, "DNAm",
                           additional_covariates, "change", filtered_probes))

lme4_level <- time_step("LME4: run_stratified_analysis (level)",
  .run_stratified_analysis(lme4_pheno_list, lme4_omics_list, "DNAm",
                           additional_covariates, "level", filtered_probes))

cat("  LME4 change coefficients: ", nrow(lme4_change$all$coefficients), " rows, ",
    length(unique(lme4_change$all$coefficients$ANALYTE_NAME)), " analytes\n", sep = "")
cat("  LME4 level  coefficients: ", nrow(lme4_level$all$coefficients),  " rows, ",
    length(unique(lme4_level$all$coefficients$ANALYTE_NAME)),  " analytes\n", sep = "")


# =============================================================================
# TIMING SUMMARY
# =============================================================================

print_timing_summary()
