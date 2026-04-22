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
#   list(all, male, female)
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
  }
  if (!is.null(results$treatment_effects)) {
    te <- results$treatment_effects
    cat("  Treatment effects: ", nrow(te), " rows\n")
    cat("  FU levels:         ", paste(sort(unique(te$FU)), collapse = ", "), "\n")
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
                      grepl("^CONTROL_STATUS", coefs$COEFFICIENT), ]

  cat("  Sample analyte: ", sample_analyte, "\n")
  cat("  Treatment effect: ", round(te_row$EFFECT_SIZE, 6),
      " CONTROL_STATUS coef: ", round(coef_row$EFFECT_SIZE, 6), "\n\n")

  list(
    "Effect size matches CONTROL_STATUS" = abs(te_row$EFFECT_SIZE - coef_row$EFFECT_SIZE) < 1e-10,
    "SE matches CONTROL_STATUS" = abs(te_row$SE - coef_row$SE) < 1e-10,
    "P-value matches CONTROL_STATUS" = abs(te_row$P_VALUE - coef_row$P_VALUE) < 1e-10
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
  ctrl_coef <- coefs_analyte[grepl("^CONTROL_STATUS[0-9]*$", coefs_analyte$COEFFICIENT), ]
  interaction_coef <- coefs_analyte[grepl("CONTROL_STATUS.*:.*FU", coefs_analyte$COEFFICIENT), ]

  # Treatment effect at FU=1 should match CONTROL_STATUS coefficient
  te_fu1 <- te_analyte[te_analyte$FU == 1, ]
  effect_fu1_expected <- ctrl_coef$EFFECT_SIZE

  # Treatment effect at FU=2 should equal CONTROL_STATUS + interaction
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
    "FU=1 effect matches CONTROL_STATUS" = abs(te_fu1$EFFECT_SIZE - effect_fu1_expected) < 1e-6,
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

#' Check reports structure: pheno, omics, covariates, randomization summaries
#' @param reports Reports object from FAST_omics_WAS_reports
#' @return Named list of check results
check_reports_structure <- function(reports) {
  report <- reports$all
  list(
    "Reports all stratum exists" = !is.null(report),
    "Has pheno_summary" = !is.null(report$pheno_summary),
    "Has omics_summary" = !is.null(report$omics_summary),
    "Has randomization_summary" = !is.null(report$randomization_summary),
    "Pheno summary has rows" = nrow(report$pheno_summary) > 0,
    "Omics summary has rows" = nrow(report$omics_summary) > 0,
    "Randomization summary has rows" = nrow(report$randomization_summary) > 0,
    "Reports sex stratification works" = !is.null(reports$male) && !is.null(reports$female)
  )
}

# ===== DATA SETUP =====

cat("Test Suite for main.R\n")
cat("==========================\n\n")

# Load raw data
pheno_raw <- readRDS("PracticeData/pheno_example.rds")
omics_raw <- readRDS("PracticeData/synth_small_betas.rds")

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

# Prepare pheno data
pheno <- pheno_raw[!duplicated(pheno_raw$SAMPLE_ID), ]
pheno <- pheno[complete.cases(pheno[, c("SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "CONTROL_STATUS")]), ]

# Prepare omics data - ensure overlap with both probe lists for DNAm tests
analyte_names <- rownames(omics_raw)
sample_names <- colnames(omics_raw)
omics_full <- as.data.frame(omics_raw)
colnames(omics_full) <- sample_names
omics_full <- cbind(ANALYTE_NAME = analyte_names, omics_full, stringsAsFactors = FALSE)

# Load probe lists to ensure test data has overlap with both
full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")

# Find probes in our data that are in each list
in_full <- omics_full$ANALYTE_NAME %in% full_probes
in_filtered <- omics_full$ANALYTE_NAME %in% filtered_probes

cat("Probe overlap in raw data:\n")
cat("  In full list:     ", sum(in_full), "\n")
cat("  In filtered list: ", sum(in_filtered), "\n\n")

# Create test subset: 100 from filtered (if available) + 400 from full-only
filtered_probes_available <- omics_full[in_filtered, ]
full_only_probes <- omics_full[in_full & !in_filtered, ]

n_filtered <- min(100, nrow(filtered_probes_available))
n_full_only <- min(400, nrow(full_only_probes))

omics <- rbind(
  filtered_probes_available[1:n_filtered, ],
  full_only_probes[1:n_full_only, ]
)

cat("Test subset created:\n")
cat("  From filtered list: ", n_filtered, "\n")
cat("  From full-only:     ", n_full_only, "\n")
cat("  Total:              ", nrow(omics), "\n\n")

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
  n_cores = 4
)
reports_test1 <- FAST_omics_WAS_reports(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
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
  n_cores = 4
)
reports_test2 <- FAST_omics_WAS_reports(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
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
  n_cores = 4
)
reports_test3 <- FAST_omics_WAS_reports(
  pheno = pheno_single_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
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
  n_cores = 4
)
reports_test4 <- FAST_omics_WAS_reports(
  pheno = pheno_multi_fu,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
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

all_pass <- test1_pass && test2_pass && test3_pass && test4_pass
if (all_pass) {
  cat("ALL TESTS PASSED\n")
} else {
  cat("SOME TESTS FAILED\n")
}
