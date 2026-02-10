# 
# Compare p-values and effect sizes across analysis methods
# Borrows test data setup from test_comprehensive.R
#

source("main.R")
require(lme4)
require(emmeans)
require(limma)

# ===== SETUP =====

cat("Method Comparison Analysis\n")
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

# Prepare omics data
analyte_names <- rownames(omics_raw)
sample_names <- colnames(omics_raw)
omics <- as.data.frame(omics_raw)
colnames(omics) <- sample_names
omics <- cbind(ANALYTE_NAME = analyte_names, omics, stringsAsFactors = FALSE)
omics <- omics[1:50, ]  # Subsample for testing

additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")

cat("Data Loaded\n")
cat("  Pheno:      ", nrow(pheno), " samples\n")
cat("  Omics:      ", nrow(omics), " analytes\n\n")

# ===== RUN ANALYSES =====

cat("Running analyses...\n\n")

# TEST 1: LM (Single FU)
pheno_lm <- pheno[pheno$FU %in% c(0, 1), ]
cat("TEST 1: LM (Single FU, N=", nrow(pheno_lm), ")\n")
results_lm <- FAST_omics_WAS(pheno = pheno_lm, omics = omics, omics_type = "Proteomics")
cat("  Results: ", nrow(results_lm$all$results), " rows\n\n")

# TEST 2: LME4 (Multiple FU)
cat("TEST 2: LME4 (Multiple FU, N=", nrow(pheno), ")\n")
results_lme4 <- FAST_omics_WAS(pheno = pheno, omics = omics, omics_type = "Proteomics")
cat("  Results: ", nrow(results_lme4$all$results), " rows\n\n")

# TEST 3: LIMMA (Multiple FU, DNAm)
cat("TEST 3: LIMMA (Multiple FU, N=", nrow(pheno), ")\n")
results_limma <- FAST_omics_WAS(pheno = pheno, omics = omics, omics_type = "DNAm")
cat("  Results: ", nrow(results_limma$all$results), " rows\n\n")

# ===== COMPARISON =====

cat("P-value Correlations for Shared Coefficients\n")
cat("============================================\n\n")

# Extract CONTROL_STATUS results from each method
lm_control <- subset(results_lm$all$results, COEFFICIENT == "CONTROL_STATUS")
lme4_control <- subset(results_lme4$all$results, COEFFICIENT == "CONTROL_STATUS")
limma_control <- subset(results_limma$all$results, COEFFICIENT == "CONTROL_STATUS")

# Merge on ANALYTE_NAME to ensure same analytes
comparison_df <- merge(lm_control[, c("ANALYTE_NAME", "P_VALUE", "EFFECT_SIZE")],
                       lme4_control[, c("ANALYTE_NAME", "P_VALUE", "EFFECT_SIZE")],
                       by = "ANALYTE_NAME", suffixes = c("_lm", "_lme4"))

comparison_df <- merge(comparison_df,
                       limma_control[, c("ANALYTE_NAME", "P_VALUE", "EFFECT_SIZE")],
                       by = "ANALYTE_NAME")
colnames(comparison_df)[colnames(comparison_df) == "P_VALUE"] <- "P_VALUE_limma"
colnames(comparison_df)[colnames(comparison_df) == "EFFECT_SIZE"] <- "EFFECT_SIZE_limma"

cat("CONTROL_STATUS Coefficient:\n")
cat("  N analytes with data: ", nrow(comparison_df), "\n\n")

# P-value correlations
cat("P-value Correlations:\n")
cor_lm_lme4_p <- cor(comparison_df$P_VALUE_lm, comparison_df$P_VALUE_lme4, use = "complete.obs")
cor_lm_limma_p <- cor(comparison_df$P_VALUE_lm, comparison_df$P_VALUE_limma, use = "complete.obs")
cor_lme4_limma_p <- cor(comparison_df$P_VALUE_lme4, comparison_df$P_VALUE_limma, use = "complete.obs")

cat("  LM vs LME4:   ", sprintf("%.4f", cor_lm_lme4_p), "\n")
cat("  LM vs LIMMA:  ", sprintf("%.4f", cor_lm_limma_p), "\n")
cat("  LME4 vs LIMMA:", sprintf("%.4f", cor_lme4_limma_p), "\n\n")

# Effect size correlations
cat("Effect Size Correlations:\n")
cor_lm_lme4_es <- cor(comparison_df$EFFECT_SIZE_lm, comparison_df$EFFECT_SIZE_lme4, use = "complete.obs")
cor_lm_limma_es <- cor(comparison_df$EFFECT_SIZE_lm, comparison_df$EFFECT_SIZE_limma, use = "complete.obs")
cor_lme4_limma_es <- cor(comparison_df$EFFECT_SIZE_lme4, comparison_df$EFFECT_SIZE_limma, use = "complete.obs")

cat("  LM vs LME4:   ", sprintf("%.4f", cor_lm_lme4_es), "\n")
cat("  LM vs LIMMA:  ", sprintf("%.4f", cor_lm_limma_es), "\n")
cat("  LME4 vs LIMMA:", sprintf("%.4f", cor_lme4_limma_es), "\n\n")

# Sample comparison (first 5 analytes)
cat("Sample Results (first 5 analytes, CONTROL_STATUS):\n")
sample_comparison <- comparison_df[1:5, c("ANALYTE_NAME", "P_VALUE_lm", "P_VALUE_lme4", "P_VALUE_limma",
                                           "EFFECT_SIZE_lm", "EFFECT_SIZE_lme4", "EFFECT_SIZE_limma")]
print(sample_comparison)
cat("\n")

# Range of differences
cat("P-value Difference Statistics:\n")
diff_lm_lme4 <- abs(comparison_df$P_VALUE_lm - comparison_df$P_VALUE_lme4)
diff_lm_limma <- abs(comparison_df$P_VALUE_lm - comparison_df$P_VALUE_limma)
diff_lme4_limma <- abs(comparison_df$P_VALUE_lme4 - comparison_df$P_VALUE_limma)

cat("  LM vs LME4:    Mean diff = ", sprintf("%.6f", mean(diff_lm_lme4, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_lm_lme4, na.rm = TRUE)), "\n")
cat("  LM vs LIMMA:   Mean diff = ", sprintf("%.6f", mean(diff_lm_limma, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_lm_limma, na.rm = TRUE)), "\n")
cat("  LME4 vs LIMMA: Mean diff = ", sprintf("%.6f", mean(diff_lme4_limma, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_lme4_limma, na.rm = TRUE)), "\n\n")

# Effect size difference statistics
cat("Effect Size Difference Statistics:\n")
diff_es_lm_lme4 <- abs(comparison_df$EFFECT_SIZE_lm - comparison_df$EFFECT_SIZE_lme4)
diff_es_lm_limma <- abs(comparison_df$EFFECT_SIZE_lm - comparison_df$EFFECT_SIZE_limma)
diff_es_lme4_limma <- abs(comparison_df$EFFECT_SIZE_lme4 - comparison_df$EFFECT_SIZE_limma)

cat("  LM vs LME4:    Mean diff = ", sprintf("%.6f", mean(diff_es_lm_lme4, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_es_lm_lme4, na.rm = TRUE)), "\n")
cat("  LM vs LIMMA:   Mean diff = ", sprintf("%.6f", mean(diff_es_lm_limma, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_es_lm_limma, na.rm = TRUE)), "\n")
cat("  LME4 vs LIMMA: Mean diff = ", sprintf("%.6f", mean(diff_es_lme4_limma, na.rm = TRUE)),
    ", Max diff = ", sprintf("%.6f", max(diff_es_lme4_limma, na.rm = TRUE)), "\n\n")

cat("Comparison Complete\n")
