# ===== TEST FOR MAIN.R =====
#
# Tests the 2x2 matrix of analysis scenarios:
#
#                    |  Single FU (LM)  |  Multiple FU (LME4)
# -------------------|------------------|---------------------
# Non-DNAm           |  Test 1           |  Test 2
# (Proteomics/Metab) |  LM regression    |  LME4 mixed effects
# -------------------|------------------|---------------------
# DNAm               |  Test 3           |  Test 4
#                    |  LM               |  LME4
#
# Each test runs both change and level models (change = analyte change from baseline,
# level = absolute analyte level at follow-up). Analysis results are nested as:
#   list(analysis_change = list(all, male, female),
#        analysis_level  = list(all, male, female))
# Reports are generated separately via FAST_omics_WAS_reports, returning:
#   list(
#     pheno_summary         = data.frame(...),   # study-level
#     variable_summaries    = list(              # per stratum x FU x Tx cell
#       all    = list(omics_FU0_Tx0, omics_FU0_Tx1, ..., covariates_FU0_Tx0, ...),
#       male   = list(...),
#       female = list(...)
#     ),
#     randomization_reports = list(              # study-level
#       omics_report     = data.frame(...),
#       covariate_report = data.frame(...)
#     )
#   )
#
# DNAm results have the same structure as non-DNAm, with an additional
# BH_P_VALUE_FILTERED column that contains BH-corrected p-values for
# filtered probes only (NA for non-filtered probes).

source("main.R")
require(lme4)


# ===== TESTING HELPERS =====

#' Run a list of validation checks and print results
#' @param checks Named list of logical values
#' @return TRUE if all checks pass, FALSE otherwise
run_checks <- function(checks) {
  all_pass <- TRUE
  for (name in names(checks)) {
    status <- if (checks[[name]]) "PASS" else "FAIL"
    cat(status, " ", name, "\n", sep = "")
    if (!checks[[name]]) all_pass <- FALSE
  }
  all_pass
}

#' Print summary of analysis results
#' @param results Results object (with coefficients and treatment_effects)
#' @param label Label for output
print_results_summary <- function(results, label = "") {
  if (label != "") cat(label, "\n", sep = "")
  if (!is.null(results$coefficients)) {
    res <- results$coefficients
    cat("  Coefficients rows: ", nrow(res), "\n")
    cat("  Analytes:          ", length(unique(res$ANALYTE_NAME)), "\n")
    cat("  Columns:           ", paste(colnames(res), collapse = ", "), "\n")
    fu_terms <- unique(res$COEFFICIENT[grepl("(^|:)FU", res$COEFFICIENT)])
    if (length(fu_terms)) cat("  FU terms:          ", paste(fu_terms, collapse = ", "), "\n")
  }
  if (!is.null(results$treatment_effects)) {
    te <- results$treatment_effects
    cat("  Treatment effects: ", nrow(te), " rows\n")
    cat("  FU levels:         ", paste(sort(unique(te$FU)), collapse = ", "), "\n")
    if ("FU_LABEL" %in% names(te)) {
      cat("  FU labels:         ", paste(unique(paste0(te$FU, " = ", te$FU_LABEL)), collapse = ", "), "\n")
    }
  }
}

#' Validate treatment effects for single FU analysis
#' @param results Results object with $coefficients and $treatment_effects
#' @return List of check results
validate_treatment_effects_single_fu <- function(results) {
  coefs <- results$coefficients
  te <- results$treatment_effects

  sample_analyte <- te$ANALYTE_NAME[1]
  te_row <- te[te$ANALYTE_NAME == sample_analyte, ]
  coef_row <- coefs[coefs$ANALYTE_NAME == sample_analyte &
                      grepl("^TREATMENT_GROUP", coefs$COEFFICIENT), ]

  cat("  Sample analyte: ", sample_analyte, "\n")
  cat("  Treatment effect: ", round(te_row$EFFECT_SIZE, 6),
      " TREATMENT_GROUP coef: ", round(coef_row$EFFECT_SIZE, 6), "\n\n")

  list(
    "Effect size matches TREATMENT_GROUP" = abs(te_row$EFFECT_SIZE - coef_row$EFFECT_SIZE) < 1e-10,
    "SE matches TREATMENT_GROUP" = abs(te_row$SE - coef_row$SE) < 1e-10,
    "P-value matches TREATMENT_GROUP" = abs(te_row$P_VALUE - coef_row$P_VALUE) < 1e-10
  )
}

#' Validate treatment effects for multiple FU analysis
#' @param results Results object with $coefficients and $treatment_effects
#' @param method Label for output ("LME4", etc.)
#' @return List of check results
validate_treatment_effects_multi_fu <- function(results, method = "LME4") {
  coefs <- results$coefficients
  te <- results$treatment_effects

  sample_analyte <- te$ANALYTE_NAME[1]
  te_analyte <- te[te$ANALYTE_NAME == sample_analyte, ]
  coefs_analyte <- coefs[coefs$ANALYTE_NAME == sample_analyte, ]

  # Get relevant coefficients
  ctrl_coef <- coefs_analyte[grepl("^TREATMENT_GROUP[0-9]*$", coefs_analyte$COEFFICIENT), ]
  interaction_coef <- coefs_analyte[grepl("TREATMENT_GROUP.*:.*FU", coefs_analyte$COEFFICIENT), ]

  # Treatment effect at FU=1 should match TREATMENT_GROUP coefficient
  te_fu1 <- te_analyte[te_analyte$FU == 1, ]
  effect_fu1_expected <- ctrl_coef$EFFECT_SIZE

  # Treatment effect at FU=2 should equal TREATMENT_GROUP + interaction
  te_fu2 <- te_analyte[te_analyte$FU == 2, ]
  effect_fu2_expected <- ctrl_coef$EFFECT_SIZE + interaction_coef$EFFECT_SIZE

  # Naive SE (assuming independence) for comparison
  naive_se_fu2 <- sqrt(ctrl_coef$SE^2 + interaction_coef$SE^2)

  cat("  Sample analyte: ", sample_analyte, "\n")
  cat("  FU=1: ", method, "=", round(te_fu1$EFFECT_SIZE, 6),
      " expected=", round(effect_fu1_expected, 6), "\n")
  cat("  FU=2: ", method, "=", round(te_fu2$EFFECT_SIZE, 6),
      " expected=", round(effect_fu2_expected, 6), "\n")
  cat("  FU=2 SE: ", method, "=", round(te_fu2$SE, 6),
      " naive=", round(naive_se_fu2, 6), "\n\n")

  list(
    "FU=1 effect matches TREATMENT_GROUP" = abs(te_fu1$EFFECT_SIZE - effect_fu1_expected) < 1e-6,
    "FU=2 effect matches coefficient sum" = abs(te_fu2$EFFECT_SIZE - effect_fu2_expected) < 1e-6,
    "FU=2 SE differs from naive (uses covariance)" = abs(te_fu2$SE - naive_se_fu2) > 1e-10
  )
}

#' Standard structural checks for non-DNAm analysis results (one response type)
#' @param analysis Analysis sub-object (e.g. results$analysis_change) with all/male/female strata
#' @param expected_fu_levels Expected FU levels in treatment effects
#' @return Named list of check results
check_non_dnam_structure <- function(analysis, expected_fu_levels) {
  list(
    "Coefficients exist" = !is.null(analysis$all$coefficients),
    "Coefficients have rows" = nrow(analysis$all$coefficients) > 0,
    "Coefficients have N_OBS" = "N_OBS" %in% colnames(analysis$all$coefficients),
    "Coefficients have BH correction" = "BH_P_VALUE" %in% colnames(analysis$all$coefficients),
    "Treatment effects exist" = !is.null(analysis$all$treatment_effects),
    "Treatment effects have rows" = nrow(analysis$all$treatment_effects) > 0,
    "Treatment effects FU levels correct" = all(sort(unique(analysis$all$treatment_effects$FU)) %in% expected_fu_levels),
    "Treatment effects have BH correction" = "BH_P_VALUE" %in% colnames(analysis$all$treatment_effects),
    "Sex stratification works" = !is.null(analysis$male$coefficients) && !is.null(analysis$female$coefficients),
    "Output structure correct" = all(c("all", "male", "female") %in% names(analysis))
  )
}

#' Standard structural checks for DNAm analysis results (flat structure with BH_P_VALUE_FILTERED column)
#' @param analysis Analysis sub-object (e.g. results$analysis_change) with all/male/female strata
#' @param expected_fu_levels Expected FU levels in treatment effects
#' @param filtered_probes Vector of filtered probe names for validation
#' @return Named list of check results
check_dnam_structure <- function(analysis, expected_fu_levels, filtered_probes) {
  coefs <- analysis$all$coefficients
  te <- analysis$all$treatment_effects

  filtered_in_results <- coefs$ANALYTE_NAME %in% filtered_probes
  filtered_have_bh <- all(!is.na(coefs$BH_P_VALUE_FILTERED[filtered_in_results]))
  non_filtered_have_na <- all(is.na(coefs$BH_P_VALUE_FILTERED[!filtered_in_results]))

  list(
    "Has all/male/female strata" = all(c("all", "male", "female") %in% names(analysis)),
    "Coefficients exist" = !is.null(coefs),
    "Coefficients have rows" = nrow(coefs) > 0,
    "Has effect sizes" = all(!is.na(coefs$EFFECT_SIZE)),
    "Coef has N_OBS" = "N_OBS" %in% colnames(coefs),
    "Coef has BH_P_VALUE" = "BH_P_VALUE" %in% colnames(coefs),
    "Coef has BH_P_VALUE_FILTERED" = "BH_P_VALUE_FILTERED" %in% colnames(coefs),
    "Filtered probes have BH_P_VALUE_FILTERED" = filtered_have_bh,
    "Non-filtered probes have NA BH_P_VALUE_FILTERED" = non_filtered_have_na,
    "Treatment effects exist" = !is.null(te),
    "Treatment effects have rows" = nrow(te) > 0,
    "TE has BH_P_VALUE" = "BH_P_VALUE" %in% colnames(te),
    "TE has BH_P_VALUE_FILTERED" = "BH_P_VALUE_FILTERED" %in% colnames(te),
    "TE FU levels correct" = all(sort(unique(te$FU)) %in% expected_fu_levels),
    "Sex stratification works" = !is.null(analysis$male$coefficients) && !is.null(analysis$female$coefficients)
  )
}

#' Check top-level structure: has analysis_change and analysis_level
#' @param results Full results object from FAST_omics_WAS
#' @return Named list of check results
check_top_level_structure <- function(results) {
  list(
    "Has analysis_change and analysis_level" = all(c("analysis_change", "analysis_level") %in% names(results)),
    "Does not contain reports" = !"reports" %in% names(results),
    "analysis_change has all/male/female" = all(c("all", "male", "female") %in% names(results$analysis_change)),
    "analysis_level has all/male/female" = all(c("all", "male", "female") %in% names(results$analysis_level))
  )
}

#' Check reports structure: pheno_summary, variable_summaries, randomization_reports
#' @param reports Reports object from FAST_omics_WAS_reports
#' @param has_covariates Whether additional_covariates were passed (default TRUE)
#' @return Named list of check results
check_reports_structure <- function(reports, has_covariates = TRUE) {
  vs <- reports$variable_summaries
  rr <- reports$randomization_reports

  all_keys   <- names(vs$all)
  omics_keys <- all_keys[startsWith(all_keys, "omics_")]
  covar_keys <- all_keys[startsWith(all_keys, "covariates_")]

  list(
    "Has pheno_summary"                  = !is.null(reports$pheno_summary),
    "Pheno summary has rows"             = !is.null(reports$pheno_summary) && nrow(reports$pheno_summary) > 0,
    "Has variable_summaries"             = !is.null(vs),
    "variable_summaries has all stratum" = !is.null(vs$all),
    "Baseline omics entries present"     = all(c("omics_FU0_Tx0", "omics_FU0_Tx1") %in% omics_keys),
    "Omics entries have rows"            = length(omics_keys) > 0 && all(sapply(omics_keys, function(k) nrow(vs$all[[k]]) > 0)),
    "Covariate entries present"          = !has_covariates || length(covar_keys) > 0,
    "Covariate entries have rows"        = !has_covariates || all(sapply(covar_keys, function(k) nrow(vs$all[[k]]) > 0)),
    "Sex stratification works"           = !is.null(vs$male) && !is.null(vs$female),
    "Has randomization_reports"          = !is.null(rr),
    "Has omics_report"                   = !is.null(rr$omics_report),
    "Omics report has rows"              = !is.null(rr$omics_report) && nrow(rr$omics_report) > 0,
    "Has covariate_report"               = !has_covariates || !is.null(rr$covariate_report),
    "Covariate report has rows"          = !has_covariates || (!is.null(rr$covariate_report) && nrow(rr$covariate_report) > 0)
  )
}

# ===== DATA SETUP =====

cat("Test Suite for main.R\n")
cat("==========================\n\n")

# Small, deterministic fixtures with both filtered and full-only DNAm probes.
source("tests/synthetic_data.R")
fixture <- .make_synthetic_test_data(n_analytes = 16L)
pheno <- fixture$pheno
omics <- fixture$omics
full_probes <- fixture$full_probes
filtered_probes <- fixture$filtered_probes

cat("Synthetic probe coverage:\n")
cat("  In full list:     ", sum(omics$ANALYTE_NAME %in% full_probes), "\n")
cat("  In filtered list: ", sum(omics$ANALYTE_NAME %in% filtered_probes), "\n\n")

additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")

cat("Data Loaded\n")
cat("  Pheno:     ", nrow(pheno), "samples, FU:", paste(sort(unique(pheno$FU)), collapse = ", "), "\n")
cat("  Omics:     ", nrow(omics), "analytes\n")
cat("  Covariates:", paste(additional_covariates, collapse = ", "), "\n")
cat("  (Omics subset ensures overlap with both DNAm probe lists)\n\n")

# Create single-FU and multi-FU pheno subsets
pheno_single_fu <- pheno[pheno$FU %in% c(0, 1), ]
pheno_single_fu <- pheno_single_fu[!duplicated(pheno_single_fu$SAMPLE_ID), ]
pheno_multi_fu <- pheno  # All FU levels
fu_labels_single <- c(`0` = "Baseline", `1` = "12mo")
fu_labels_multi <- c(`0` = "Baseline", `1` = "12mo", `2` = "24mo")

# Reuse the four labeled fits below for unlabeled/relabeling comparisons. R's session tempdir
# owns these files; the label section removes them after successful checks.
label_checkpoints <- tempfile("track111_comprehensive_")

# =============================================================================
# TEST 1: Non-DNAm + Single FU (uses LM)
# =============================================================================

cat("TEST 1: Non-DNAm + Single FU\n")
cat("============================\n")
cat("Method: Linear regression (LM)\n\n")

cat("Running analysis...\n")
results_test1 <- FAST_omics_WAS(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 4,
  fu_labels = fu_labels_single,
  checkpoint_dir = file.path(label_checkpoints, "1")
)
reports_test1 <- FAST_omics_WAS_reports(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  fu_labels = fu_labels_single
)

cat("Top-Level Structure Checks\n")
test1_top_pass <- run_checks(check_top_level_structure(results_test1))
cat("\n")

print_results_summary(results_test1$analysis_change$all, "Change Results")
print_results_summary(results_test1$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test1_change_struct_pass <- run_checks(check_non_dnam_structure(results_test1$analysis_change, expected_fu_levels = 1))
cat("\n")

cat("Level Structural Checks\n")
test1_level_struct_pass <- run_checks(check_non_dnam_structure(results_test1$analysis_level, expected_fu_levels = 1))
cat("\n")

cat("Reports Structural Checks\n")
test1_reports_pass <- run_checks(check_reports_structure(reports_test1))
cat("\n")

cat("Change Treatment Effect Validation\n")
test1_change_te_pass <- run_checks(validate_treatment_effects_single_fu(results_test1$analysis_change$all))
cat("\n")

cat("Level Treatment Effect Validation\n")
test1_level_te_pass <- run_checks(validate_treatment_effects_single_fu(results_test1$analysis_level$all))
cat("\n")

test1_pass <- test1_top_pass && test1_change_struct_pass && test1_level_struct_pass &&
              test1_reports_pass && test1_change_te_pass && test1_level_te_pass

# =============================================================================
# TEST 2: Non-DNAm + Multiple FU (uses LME4)
# =============================================================================

cat("TEST 2: Non-DNAm + Multiple FU\n")
cat("==============================\n")
cat("Method: Linear mixed effects (LME4)\n\n")

cat("Running analysis...\n")
results_test2 <- FAST_omics_WAS(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 4,
  fu_labels = fu_labels_multi,
  checkpoint_dir = file.path(label_checkpoints, "2")
)
reports_test2 <- FAST_omics_WAS_reports(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  fu_labels = fu_labels_multi
)

cat("Top-Level Structure Checks\n")
test2_top_pass <- run_checks(check_top_level_structure(results_test2))
cat("\n")

print_results_summary(results_test2$analysis_change$all, "Change Results")
print_results_summary(results_test2$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test2_change_struct_pass <- run_checks(check_non_dnam_structure(results_test2$analysis_change, expected_fu_levels = c(1, 2)))
cat("\n")

cat("Level Structural Checks\n")
test2_level_struct_pass <- run_checks(check_non_dnam_structure(results_test2$analysis_level, expected_fu_levels = c(1, 2)))
cat("\n")

cat("Reports Structural Checks\n")
test2_reports_pass <- run_checks(check_reports_structure(reports_test2))
cat("\n")

cat("Change Treatment Effect Validation\n")
test2_change_te_pass <- run_checks(validate_treatment_effects_multi_fu(results_test2$analysis_change$all, method = "emmeans"))
cat("\n")

cat("Level Treatment Effect Validation\n")
test2_level_te_pass <- run_checks(validate_treatment_effects_multi_fu(results_test2$analysis_level$all, method = "emmeans"))
cat("\n")

test2_pass <- test2_top_pass && test2_change_struct_pass && test2_level_struct_pass &&
              test2_reports_pass && test2_change_te_pass && test2_level_te_pass

# =============================================================================
# TEST 3: DNAm + Single FU (uses LM, with BH_P_VALUE_FILTERED column)
# =============================================================================

cat("TEST 3: DNAm + Single FU\n")
cat("========================\n")
cat("Method: Linear regression (LM) with BH_P_VALUE_FILTERED for filtered probes\n\n")

cat("Running analysis (expects warnings about missing probes)...\n")
results_test3 <- FAST_omics_WAS(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  n_cores = 4,
  fu_labels = fu_labels_single,
  checkpoint_dir = file.path(label_checkpoints, "3")
)
reports_test3 <- FAST_omics_WAS_reports(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  fu_labels = fu_labels_single
)

cat("\nResults Structure\n")
cat("  Top-level keys: ", paste(names(results_test3), collapse = ", "), "\n")

cat("\nChange Results\n")
print_results_summary(results_test3$analysis_change$all, "Change: All probes (with BH_P_VALUE_FILTERED column)")
n_filtered_bh <- sum(!is.na(results_test3$analysis_change$all$coefficients$BH_P_VALUE_FILTERED))
n_total <- length(unique(results_test3$analysis_change$all$coefficients$ANALYTE_NAME))
cat("  Probes with BH_P_VALUE_FILTERED: ", n_filtered_bh, " / ", n_total, "\n")

cat("\nLevel Results\n")
print_results_summary(results_test3$analysis_level$all, "Level: All probes (with BH_P_VALUE_FILTERED column)")
n_filtered_bh <- sum(!is.na(results_test3$analysis_level$all$coefficients$BH_P_VALUE_FILTERED))
n_total <- length(unique(results_test3$analysis_level$all$coefficients$ANALYTE_NAME))
cat("  Probes with BH_P_VALUE_FILTERED: ", n_filtered_bh, " / ", n_total, "\n")
cat("\n")

cat("Top-Level Structure Checks\n")
test3_top_pass <- run_checks(check_top_level_structure(results_test3))
cat("\n")

cat("Change Structural Checks\n")
test3_change_struct_pass <- run_checks(check_dnam_structure(results_test3$analysis_change, expected_fu_levels = 1, filtered_probes))
cat("\n")

cat("Level Structural Checks\n")
test3_level_struct_pass <- run_checks(check_dnam_structure(results_test3$analysis_level, expected_fu_levels = 1, filtered_probes))
cat("\n")

cat("Reports Structural Checks\n")
test3_reports_pass <- run_checks(check_reports_structure(reports_test3))
cat("\n")

cat("Change Treatment Effect Validation\n")
test3_change_te_pass <- run_checks(validate_treatment_effects_single_fu(results_test3$analysis_change$all))
cat("\n")

cat("Level Treatment Effect Validation\n")
test3_level_te_pass <- run_checks(validate_treatment_effects_single_fu(results_test3$analysis_level$all))
cat("\n")

test3_pass <- test3_top_pass && test3_change_struct_pass && test3_level_struct_pass &&
              test3_reports_pass && test3_change_te_pass && test3_level_te_pass

# =============================================================================
# TEST 4: DNAm + Multiple FU (uses LME4, with BH_P_VALUE_FILTERED column)
# =============================================================================

cat("TEST 4: DNAm + Multiple FU\n")
cat("==========================\n")
cat("Method: Linear mixed effects (LME4) with BH_P_VALUE_FILTERED for filtered probes\n\n")

cat("Running analysis (expects warnings about missing probes)...\n")
results_test4 <- FAST_omics_WAS(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  n_cores = 4,
  fu_labels = fu_labels_multi,
  checkpoint_dir = file.path(label_checkpoints, "4")
)
reports_test4 <- FAST_omics_WAS_reports(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  fu_labels = fu_labels_multi
)

cat("\nResults Structure\n")
cat("  Top-level keys: ", paste(names(results_test4), collapse = ", "), "\n")

cat("\nChange Results\n")
print_results_summary(results_test4$analysis_change$all, "Change: All probes (with BH_P_VALUE_FILTERED column)")
n_filtered_bh <- sum(!is.na(results_test4$analysis_change$all$coefficients$BH_P_VALUE_FILTERED))
n_total <- length(unique(results_test4$analysis_change$all$coefficients$ANALYTE_NAME))
cat("  Probes with BH_P_VALUE_FILTERED: ", n_filtered_bh, " / ", n_total, "\n")

cat("\nLevel Results\n")
print_results_summary(results_test4$analysis_level$all, "Level: All probes (with BH_P_VALUE_FILTERED column)")
n_filtered_bh <- sum(!is.na(results_test4$analysis_level$all$coefficients$BH_P_VALUE_FILTERED))
n_total <- length(unique(results_test4$analysis_level$all$coefficients$ANALYTE_NAME))
cat("  Probes with BH_P_VALUE_FILTERED: ", n_filtered_bh, " / ", n_total, "\n")
cat("\n")

cat("Top-Level Structure Checks\n")
test4_top_pass <- run_checks(check_top_level_structure(results_test4))
cat("\n")

cat("Change Structural Checks\n")
test4_change_struct_pass <- run_checks(check_dnam_structure(results_test4$analysis_change, expected_fu_levels = c(1, 2), filtered_probes))
cat("\n")

cat("Level Structural Checks\n")
test4_level_struct_pass <- run_checks(check_dnam_structure(results_test4$analysis_level, expected_fu_levels = c(1, 2), filtered_probes))
cat("\n")

cat("Reports Structural Checks\n")
test4_reports_pass <- run_checks(check_reports_structure(reports_test4))
cat("\n")

cat("Change Treatment Effect Validation\n")
test4_change_te_pass <- run_checks(validate_treatment_effects_multi_fu(results_test4$analysis_change$all, method = "emmeans"))
cat("\n")

cat("Level Treatment Effect Validation\n")
test4_level_te_pass <- run_checks(validate_treatment_effects_multi_fu(results_test4$analysis_level$all, method = "emmeans"))
cat("\n")

test4_pass <- test4_top_pass && test4_change_struct_pass && test4_level_struct_pass &&
              test4_reports_pass && test4_change_te_pass && test4_level_te_pass

# =============================================================================
# FOLLOW-UP LABELS: public API behavior across the same four scenarios
# =============================================================================

cat("Follow-up Labels\n")
test_labels_pass <- local({
  labels <- fu_labels_multi
  labeled_results <- list(results_test1, results_test2, results_test3, results_test4)
  reports <- list(reports_test1, reports_test2, reports_test3, reports_test4)
  checks <- list()
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  for (i in seq_along(labeled_results)) {
    mixed <- i %% 2L == 0L
    input <- if (mixed) pheno_multi_fu else pheno_single_fu
    type <- if (i <= 2L) "Proteomics" else "DNAm"
    visits <- if (mixed) labels else labels[1:2]
    checkpoint <- file.path(label_checkpoints, as.character(i))
    batches <- list.files(checkpoint, recursive = TRUE, full.names = TRUE)
    checksums <- tools::md5sum(batches)
    labeled <- labeled_results[[i]]
    unlabeled <- quiet(FAST_omics_WAS(input, omics, omics_type = type,
      additional_covariates = additional_covariates, n_cores = 1L,
      checkpoint_dir = checkpoint))
    reordered <- quiet(FAST_omics_WAS(input, omics, omics_type = type,
      additional_covariates = additional_covariates, n_cores = 1L,
      fu_labels = rev(visits), checkpoint_dir = checkpoint))
    prefix <- paste("Test", i, "labels:")
    checks[[paste(prefix, "metadata and map order")]] <-
      identical(labeled$metadata$fu_labels, visits) && is.null(unlabeled$metadata$fu_labels)
    checks[[paste(prefix, "fresh and cached labeled results match")]] <- identical(reordered, labeled)
    checks[[paste(prefix, "checkpoint payloads unchanged")]] <-
      length(batches) > 0L && identical(checksums, tools::md5sum(batches))

    # Strip only the expected presentation changes, then compare the entire
    # result, including estimates, SEs, status, full BH, and filtered BH.
    normalized <- labeled
    normalized$metadata$fu_labels <- NULL
    for (branch in c("analysis_change", "analysis_level")) {
      for (stratum in c("all", "male", "female")) {
        before <- unlabeled[[branch]][[stratum]]
        after <- labeled[[branch]][[stratum]]
        expected_terms <- before$coefficients$COEFFICIENT
        term_map <- c(FU2 = "FU24mo", `TREATMENT_GROUP1:FU2` = "TREATMENT_GROUP1:FU24mo")
        index <- match(expected_terms, names(term_map))
        expected_terms[!is.na(index)] <- unname(term_map[index[!is.na(index)]])
        checks[[paste(prefix, branch, stratum, "visit/interaction names")]] <-
          identical(after$coefficients$COEFFICIENT, expected_terms) &&
          identical(after$treatment_effects$FU_LABEL,
                    unname(visits[as.character(after$treatment_effects$FU)])) &&
          identical(names(after$treatment_effects)[2:3], c("FU", "FU_LABEL"))
        normalized[[branch]][[stratum]]$coefficients$COEFFICIENT <- before$coefficients$COEFFICIENT
        normalized[[branch]][[stratum]]$treatment_effects$FU_LABEL <- NULL
      }
    }
    checks[[paste(prefix, "all numerical results unchanged")]] <- identical(normalized, unlabeled)

    labeled_report <- reports[[i]]
    unlabeled_report <- quiet(FAST_omics_WAS_reports(input, omics, omics_type = type,
      additional_covariates = additional_covariates))
    checks[[paste(prefix, "report metadata")]] <- identical(labeled_report$metadata$fu_labels, visits)
    for (table in c("pheno_summary", "arm_counts")) {
      checks[[paste(prefix, table)]] <-
        identical(names(labeled_report[[table]])[1:2], c("FU", "FU_LABEL")) &&
        identical(labeled_report[[table]]$FU_LABEL,
                  unname(visits[as.character(labeled_report[[table]]$FU)]))
      labeled_report[[table]]$FU_LABEL <- NULL
    }
    labeled_report$metadata$fu_labels <- NULL
    checks[[paste(prefix, "report values and keys unchanged")]] <- identical(labeled_report, unlabeled_report)

    if (i == 2L) {
      fresh_unlabeled <- quiet(FAST_omics_WAS(input, omics, omics_type = type,
        additional_covariates = additional_covariates, n_cores = 1L))
      checks[["Independent unlabeled fit preserves numerical results"]] <-
        isTRUE(all.equal(fresh_unlabeled, unlabeled, tolerance = 1e-8))
      revised <- visits
      revised["2"] <- "Month 24 (final)"
      updated <- quiet(FAST_omics_WAS(input, omics, omics_type = type,
        additional_covariates = additional_covariates, n_cores = 1L,
        fu_labels = revised, checkpoint_dir = checkpoint))
      checks[["Changed labels reuse cached fits"]] <-
        identical(checksums, tools::md5sum(batches)) &&
        identical(updated$metadata$fu_labels, revised) &&
        "TREATMENT_GROUP1:FUMonth 24 (final)" %in% updated$analysis_change$all$coefficients$COEFFICIENT

      # Keep the plotting functions and QQ stub local to this test section.
      source("plotting_helpers.R", local = TRUE)
      checks[["Plot titles use labels and retain legacy fallback"]] <-
        grepl("24mo", plot_volcano(labeled, fu = 2)$labels$title, fixed = TRUE) &&
        grepl("FU=2", plot_volcano(unlabeled, fu = 2)$labels$title, fixed = TRUE) &&
        identical(plot_volcano(labeled, fu = 2, title = "Custom")$labels$title, "Custom")
      received_title <- NULL
      qq <- function(pvector, main, ...) received_title <<- main
      plot_qq(labeled, fu = 2)
      checks[["QQ title uses the visit label"]] <- grepl("24mo", received_title, fixed = TRUE)
    }
  }

  invalid_maps <- list(
    unnamed = unname(labels), missing_baseline = labels[-1], extra_code = c(labels, `3` = "36mo"),
    duplicate_code = setNames(unname(labels), c("0", "1", "1")),
    duplicate_label = c(`0` = "Baseline", `1` = "12mo", `2` = "12mo"),
    blank = c(`0` = "Baseline", `1` = "12mo", `2` = " "),
    missing = c(`0` = "Baseline", `1` = "12mo", `2` = NA_character_)
  )
  for (name in names(invalid_maps)) {
    error <- tryCatch({
      quiet(.validate_pheno(pheno, fu_labels = invalid_maps[[name]])); NULL
    }, error = identity)
    checks[[paste("Reject invalid follow-up map:", name)]] <-
      inherits(error, "error") && grepl("fu_labels", conditionMessage(error))
  }
  collision <- pheno
  collision$FU24mo <- collision$agebl
  error <- tryCatch({
    .label_fu_coefficients(unlabeled$analysis_change$all$coefficients,
                           collision, labels, "FU24mo"); NULL
  }, error = identity)
  checks[["Reject coefficient-name collisions"]] <-
    inherits(error, "error") && grepl("ambiguous", conditionMessage(error))
  run_checks(checks)
})
unlink(label_checkpoints, recursive = TRUE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("Summary\n")
cat("=======\n\n")

cat("Test Matrix Coverage:\n")
cat("                     |  Single FU  |  Multiple FU\n")
cat("---------------------|-------------|-------------\n")
cat("Non-DNAm (Prot/Met)  |  Test 1:", if (test1_pass) "PASS" else "FAIL",
    " |  Test 2:", if (test2_pass) "PASS" else "FAIL", "\n")
cat("DNAm (dual probes)   |  Test 3:", if (test3_pass) "PASS" else "FAIL",
    " |  Test 4:", if (test4_pass) "PASS" else "FAIL", "\n\n")

cat("Follow-up labels: ", if (test_labels_pass) "PASS" else "FAIL", "\n\n")
all_pass <- test1_pass && test2_pass && test3_pass && test4_pass && test_labels_pass
if (all_pass) {
  cat("ALL TESTS PASSED\n")
} else {
  cat("SOME TESTS FAILED\n")
}

stopifnot(all_pass)
