# ===== COMPREHENSIVE PIPELINE TEST =====
# 
# Tests all three analysis methods (LM, LME4, Limma) with the same data
# LM:    Single FU (filters to FU 0 and 1)
# LME4:  Multiple FU (keeps FU 0, 1, 2)
# Limma: Multiple FU with DNAm omics type

source("main.R")
require(lme4)
require(emmeans)
require(limma)

# ===== SETUP =====

cat("Comprehensive Pipeline Test\n")
cat("===========================\n\n")

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

# Prepare omics data
analyte_names <- rownames(omics_raw)
sample_names <- colnames(omics_raw)
omics <- as.data.frame(omics_raw)
colnames(omics) <- sample_names
omics <- cbind(ANALYTE_NAME = analyte_names, omics, stringsAsFactors = FALSE)
omics <- omics[1:50, ]  # Subsample for testing

additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")

cat("Data Loaded\n")
cat("  Pheno:     ", nrow(pheno), "samples, FU:", paste(sort(unique(pheno$FU)), collapse=", "), "\n")
cat("  Omics:     ", nrow(omics), "analytes\n")
cat("  Covariates:", paste(additional_covariates, collapse=", "), "\n\n")

# ===== TEST 1: LINEAR REGRESSION (LM) =====

cat("TEST 1: Linear Regression (LM, Single FU)\n")
cat("==========================================\n\n")

# Filter to single FU (0 and 1 only)
pheno_lm <- pheno[pheno$FU %in% c(0, 1), ]
pheno_lm <- pheno_lm[!duplicated(pheno_lm$SAMPLE_ID), ]

cat("Running LM Analysis\n")
results_lm <- FAST_omics_WAS(
  pheno = pheno_lm,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Results (LM)\n")
if (!is.null(results_lm$all$results)) {
  res <- results_lm$all$results
  cat("  Rows:     ", nrow(res), "\n")
  cat("  Analytes: ", length(unique(res$ANALYTE_NAME)), "\n")
  cat("  FU:       ", paste(sort(unique(res$FU)), collapse=", "), "\n")
  cat("  Columns:  ", paste(colnames(res), collapse=", "), "\n\n")
  
  cat("Sample Results (first 5 rows)\n")
  print(head(res, 5))
  cat("\n")
}

cat("Validation Checks (LM)\n")
checks_lm <- list(
  "Results exist" = !is.null(results_lm$all$results),
  "Has expected rows" = nrow(results_lm$all$results) > 0,
  "FU level is 1" = all(results_lm$all$results$FU == 1),
  "Has effect sizes" = all(!is.na(results_lm$all$results$EFFECT_SIZE)),
  "Has SE values" = all(!is.na(results_lm$all$results$SE)),
  "Has p-values" = all(!is.na(results_lm$all$results$P_VALUE)),
  "Has BH correction" = "BH_P_VALUE" %in% colnames(results_lm$all$results),
  "Sex stratification works" = !is.null(results_lm$male$results) && !is.null(results_lm$female$results)
)

lm_pass <- TRUE
for (name in names(checks_lm)) {
  status <- if (checks_lm[[name]]) "✓" else "✗"
  cat(status, " ", name, "\n", sep="")
  if (!checks_lm[[name]]) lm_pass <- FALSE
}
cat("\n")

# ===== TEST 2: LINEAR MIXED EFFECTS (LME4) =====

cat("TEST 2: Linear Mixed Effects (LME4, Multiple FU)\n")
cat("==============================================\n\n")

# Use full pheno data (all FU levels, including baseline)
pheno_lme4 <- pheno

cat("Running LME4 Analysis\n")
results_lme4 <- FAST_omics_WAS(
  pheno = pheno_lme4,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Results (LME4)\n")
if (!is.null(results_lme4$all$results)) {
  res <- results_lme4$all$results
  cat("  Rows:     ", nrow(res), "\n")
  cat("  Analytes: ", length(unique(res$ANALYTE_NAME)), "\n")
  cat("  FU:       ", paste(sort(unique(res$FU)), collapse=", "), "\n")
  cat("  Columns:  ", paste(colnames(res), collapse=", "), "\n\n")
  
  cat("Sample Results (first 5 rows)\n")
  print(head(res, 5))
  cat("\n")
}

cat("Validation Checks (LME4)\n")
checks_lme4 <- list(
  "Results exist" = !is.null(results_lme4$all$results),
  "Has expected rows" = nrow(results_lme4$all$results) > 0,
  "FU levels correct" = all(sort(unique(results_lme4$all$results$FU)) %in% c(1, 2)),
  "Has effect sizes" = all(!is.na(results_lme4$all$results$EFFECT_SIZE)),
  "Has SE values" = all(!is.na(results_lme4$all$results$SE)),
  "Has p-values" = all(!is.na(results_lme4$all$results$P_VALUE)),
  "Has BH correction" = "BH_P_VALUE" %in% colnames(results_lme4$all$results),
  "Sex stratification works" = !is.null(results_lme4$male$results) && !is.null(results_lme4$female$results)
)

lme4_pass <- TRUE
for (name in names(checks_lme4)) {
  status <- if (checks_lme4[[name]]) "✓" else "✗"
  cat(status, " ", name, "\n", sep="")
  if (!checks_lme4[[name]]) lme4_pass <- FALSE
}
cat("\n")

# ===== TEST 3: LIMMA (DNAm) =====

cat("TEST 3: Limma Analysis (DNAm, Multiple FU)\n")
cat("========================================\n\n")

# Use full pheno data (all FU levels, including baseline)
pheno_limma <- pheno

cat("Running Limma Analysis\n")
results_limma <- FAST_omics_WAS(
  pheno = pheno_limma,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

cat("Results (Limma)\n")
if (!is.null(results_limma$all$results)) {
  res <- results_limma$all$results
  cat("  Rows:     ", nrow(res), "\n")
  cat("  Analytes: ", length(unique(res$ANALYTE_NAME)), "\n")
  cat("  FU:       ", paste(sort(unique(res$FU)), collapse=", "), "\n")
  cat("  Columns:  ", paste(colnames(res), collapse=", "), "\n\n")
  
  cat("Sample Results (first 5 rows)\n")
  print(head(res, 5))
  cat("\n")
}

cat("Validation Checks (Limma)\n")
checks_limma <- list(
  "Results exist" = !is.null(results_limma$all$results),
  "Has expected rows" = nrow(results_limma$all$results) > 0,
  "FU levels correct" = all(sort(unique(results_limma$all$results$FU)) %in% c(1, 2)),
  "Has effect sizes" = all(!is.na(results_limma$all$results$EFFECT_SIZE)),
  "Has SE values" = all(!is.na(results_limma$all$results$SE)),
  "Has p-values" = all(!is.na(results_limma$all$results$P_VALUE)),
  "Has BH correction" = "BH_P_VALUE" %in% colnames(results_limma$all$results),
  "Sex stratification works" = !is.null(results_limma$male$results) && !is.null(results_limma$female$results)
)

limma_pass <- TRUE
for (name in names(checks_limma)) {
  status <- if (checks_limma[[name]]) "✓" else "✗"
  cat(status, " ", name, "\n", sep="")
  if (!checks_limma[[name]]) limma_pass <- FALSE
}
cat("\n")

# ===== TEST 4: LIMMA WITH SINGLE FU (EDGE CASE) =====

cat("TEST 4: Limma Analysis with Single FU (Edge Case)\n")
cat("================================================\n\n")

# Use LM pheno data (single FU: 0 and 1 only) with LIMMA
pheno_limma_single <- pheno_lm

cat("Running Limma Analysis with Single FU\n")
results_limma_single <- FAST_omics_WAS(
  pheno = pheno_limma_single,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

cat("Results (Limma with Single FU)\n")
if (!is.null(results_limma_single$all$results)) {
  res <- results_limma_single$all$results
  cat("  Rows:     ", nrow(res), "\n")
  cat("  Analytes: ", length(unique(res$ANALYTE_NAME)), "\n")
  cat("  FU:       ", paste(sort(unique(res$FU)), collapse=", "), "\n")
  cat("  Columns:  ", paste(colnames(res), collapse=", "), "\n\n")
  
  cat("Sample Results (first 5 rows)\n")
  print(head(res, 5))
  cat("\n")
}

cat("Validation Checks (Limma with Single FU)\n")
checks_limma_single <- list(
  "Results exist" = !is.null(results_limma_single$all$results),
  "Has expected rows" = nrow(results_limma_single$all$results) > 0,
  "FU level is 1" = all(results_limma_single$all$results$FU == 1),
  "Has effect sizes" = all(!is.na(results_limma_single$all$results$EFFECT_SIZE)),
  "Has SE values" = all(!is.na(results_limma_single$all$results$SE)),
  "Has p-values" = all(!is.na(results_limma_single$all$results$P_VALUE)),
  "Has BH correction" = "BH_P_VALUE" %in% colnames(results_limma_single$all$results),
  "Sex stratification works" = !is.null(results_limma_single$male$results) && !is.null(results_limma_single$female$results)
)

limma_single_pass <- TRUE
for (name in names(checks_limma_single)) {
  status <- if (checks_limma_single[[name]]) "✓" else "✗"
  cat(status, " ", name, "\n", sep="")
  if (!checks_limma_single[[name]]) limma_single_pass <- FALSE
}
cat("\n")

# ===== SUMMARY =====

cat("Summary\n")
cat("=======\n\n")

cat("Test Results:\n")
cat("  LM (Single FU):              ", if (lm_pass) "✓ PASS" else "✗ FAIL", "\n")
cat("  LME4 (Multiple FU):          ", if (lme4_pass) "✓ PASS" else "✗ FAIL", "\n")
cat("  Limma (DNAm, Multiple FU):   ", if (limma_pass) "✓ PASS" else "✗ FAIL", "\n")
cat("  Limma (DNAm, Single FU):     ", if (limma_single_pass) "✓ PASS" else "✗ FAIL", "\n\n")

if (lm_pass && lme4_pass && limma_pass && limma_single_pass) {
  cat("✓ ALL TESTS PASSED\n")
} else {
  cat("✗ SOME TESTS FAILED\n")
}
