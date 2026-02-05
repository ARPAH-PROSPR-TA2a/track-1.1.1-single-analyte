# ===== LME4 vs LIMMA COMPARISON =====
# Compare numerical results between LME4 and LIMMA methods
# Both using change scores with repeated measures handling
# LME4: (1|SUBJECT_ID) random intercept
# LIMMA: duplicateCorrelation() with block=SUBJECT_ID

source("main.R")
require(lme4)
require(emmeans)
require(limma)

cat("LME4 vs LIMMA Numerical Comparison\n")
cat("===================================\n\n")

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

# Run both analyses
cat("Running LME4 Analysis...\n")
results_lme4 <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Running LIMMA Analysis...\n")
results_limma <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

# Extract results
lme4_res <- results_lme4$all$results
limma_res <- results_limma$all$results

cat("\n\nComparison Results\n")
cat("==================\n\n")

# Merge results by ANALYTE_NAME and FU
comparison <- merge(
  lme4_res[, c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE")],
  limma_res[, c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE")],
  by = c("ANALYTE_NAME", "FU"),
  suffixes = c("_lme4", "_limma")
)

cat("Sample Comparisons (first 10 rows)\n")
print(head(comparison, 10))
cat("\n")

# Calculate correlations
cat("Correlation Analysis\n")
cat("--------------------\n")

# Effect size correlation
es_cor <- cor(comparison$EFFECT_SIZE_lme4, comparison$EFFECT_SIZE_limma, use = "complete.obs")
cat("Effect Size Correlation:  ", round(es_cor, 4), "\n")

# SE correlation
se_cor <- cor(comparison$SE_lme4, comparison$SE_limma, use = "complete.obs")
cat("SE Correlation:           ", round(se_cor, 4), "\n")

# P-value correlation (on -log scale for better comparison)
pval_cor <- cor(-log10(comparison$P_VALUE_lme4 + 1e-300), 
                 -log10(comparison$P_VALUE_limma + 1e-300), 
                 use = "complete.obs")
cat("P-value Correlation:      ", round(pval_cor, 4), "\n\n")

# Calculate differences
cat("Difference Analysis\n")
cat("-------------------\n")

# Effect size differences
es_diff <- abs(comparison$EFFECT_SIZE_lme4 - comparison$EFFECT_SIZE_limma)
cat("Effect Size Absolute Differences:\n")
cat("  Mean:   ", round(mean(es_diff, na.rm = TRUE), 6), "\n")
cat("  Median: ", round(median(es_diff, na.rm = TRUE), 6), "\n")
cat("  Max:    ", round(max(es_diff, na.rm = TRUE), 6), "\n")
cat("  Min:    ", round(min(es_diff, na.rm = TRUE), 6), "\n\n")

# SE differences
se_diff <- abs(comparison$SE_lme4 - comparison$SE_limma)
cat("SE Absolute Differences:\n")
cat("  Mean:   ", round(mean(se_diff, na.rm = TRUE), 6), "\n")
cat("  Median: ", round(median(se_diff, na.rm = TRUE), 6), "\n")
cat("  Max:    ", round(max(se_diff, na.rm = TRUE), 6), "\n")
cat("  Min:    ", round(min(se_diff, na.rm = TRUE), 6), "\n\n")

# Relative percent differences
cat("Relative Percent Differences (Effect Size)\n")
rel_diff <- abs((comparison$EFFECT_SIZE_lme4 - comparison$EFFECT_SIZE_limma) / 
                 (abs(comparison$EFFECT_SIZE_lme4) + 1e-10)) * 100
cat("  Mean:   ", round(mean(rel_diff, na.rm = TRUE), 2), "%\n")
cat("  Median: ", round(median(rel_diff, na.rm = TRUE), 2), "%\n\n")

# Summary
cat("Summary\n")
cat("=======\n")
if (es_cor > 0.95) {
  cat("✓ Effect sizes are highly correlated (r > 0.95)\n")
} else if (es_cor > 0.90) {
  cat("✓ Effect sizes are well correlated (r > 0.90)\n")
} else if (es_cor > 0.80) {
  cat("⚠ Effect sizes are moderately correlated (r > 0.80)\n")
} else {
  cat("✗ Effect sizes have low correlation (r < 0.80)\n")
}

if (mean(es_diff, na.rm = TRUE) < 0.01) {
  cat("✓ Effect sizes differ by < 0.01 on average\n")
} else if (mean(es_diff, na.rm = TRUE) < 0.05) {
  cat("✓ Effect sizes differ by < 0.05 on average\n")
} else {
  cat("⚠ Effect sizes differ by > 0.05 on average\n")
}

if (se_cor > 0.95) {
  cat("✓ SEs are highly correlated (r > 0.95)\n")
} else if (se_cor > 0.90) {
  cat("✓ SEs are well correlated (r > 0.90)\n")
} else {
  cat("⚠ SEs have lower correlation but may still be acceptable\n")
}

cat("\nModels are ", if(es_cor > 0.95 && se_cor > 0.90) "WELL ALIGNED" else "REASONABLY ALIGNED", "\n")
