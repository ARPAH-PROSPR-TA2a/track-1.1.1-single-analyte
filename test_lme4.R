# ===== TEST: LME4 PIPELINE =====
# 
# Comprehensive test of the LME4 (mixed effects) analysis pipeline
# Tests the complete workflow from raw data to statistical results

source("main.R")
require(lme4)
require(emmeans)

# ===== SETUP =====

cat("LME4 Pipeline Test\n")
cat("==================\n\n")

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

# ===== PROCESSING =====

cat("Running LME4 Analysis\n")
cat("====================\n\n")

results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

# ===== RESULTS =====

cat("Results\n")
cat("=======\n\n")

for (dataset in c("all", "male", "female")) {
  if (!is.null(results[[dataset]]$results)) {
    res <- results[[dataset]]$results
    cat(toupper(dataset), "\n")
    cat("  Rows:     ", nrow(res), "\n")
    cat("  Analytes: ", length(unique(res$ANALYTE_NAME)), "\n")
    cat("  FU:       ", paste(sort(unique(res$FU)), collapse=", "), "\n")
    cat("  Columns:  ", paste(colnames(res), collapse=", "), "\n\n")
  }
}

cat("Sample Results (first 10 rows, 'all' dataset)\n")
cat("============================================\n\n")
print(head(results$all$results, 10))

cat("\n\nValidation Checks\n")
cat("==================\n\n")

checks <- list(
  "Results exist" = !is.null(results$all$results),
  "Has expected rows" = nrow(results$all$results) > 0,
  "FU levels correct" = all(sort(unique(results$all$results$FU)) %in% c(1, 2)),
  "Has effect sizes" = all(!is.na(results$all$results$EFFECT_SIZE)),
  "Has SE values" = all(!is.na(results$all$results$SE)),
  "Has p-values" = all(!is.na(results$all$results$P_VALUE)),
  "Has BH correction" = "BH_P_VALUE" %in% colnames(results$all$results),
  "Sex stratification works" = !is.null(results$male$results) && !is.null(results$female$results)
)

all_pass <- TRUE
for (name in names(checks)) {
  status <- if (checks[[name]]) "✓" else "✗"
  cat(status, " ", name, "\n", sep="")
  if (!checks[[name]]) all_pass <- FALSE
}

cat("\n")
if (all_pass) {
  cat("✓ ALL CHECKS PASSED\n")
} else {
  cat("✗ SOME CHECKS FAILED\n")
}
