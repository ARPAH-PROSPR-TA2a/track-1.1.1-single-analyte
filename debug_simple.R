# ===== SIMPLE DEBUG: Check what the formulas are actually predicting =====

source("main.R")
require(lme4)
require(emmeans)
require(limma)

# Load data
pheno_raw <- readRDS("PracticeData/pheno_example.rds")
omics_raw <- readRDS("PracticeData/synth_small_betas.rds")

# Convert haven_labelled
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

sample_names <- colnames(omics_raw)
analyte_names <- rownames(omics_raw)
omics <- as.data.frame(omics_raw)
colnames(omics) <- sample_names
omics <- cbind(ANALYTE_NAME = analyte_names, omics, stringsAsFactors = FALSE)
omics <- omics[1:5, ]  # Just 5 analytes

cat("===== CHECKING CONTROL_STATUS ENCODING =====\n\n")
cat("CONTROL_STATUS unique values:\n")
print(table(pheno$CONTROL_STATUS))
cat("\n")
cat("CONTROL_STATUS as factor:\n")
print(contrasts(factor(pheno$CONTROL_STATUS)))
cat("\n")

# Show sample data
cat("Sample of pheno data (first 10 rows):\n")
print(pheno[1:10, c("SAMPLE_ID", "SUBJECT_ID", "FU", "CONTROL_STATUS")])
cat("\n")

# Run analyses
cat("Running full analyses...\n\n")

results_lme4 <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = "agebl"
)

results_limma <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = "agebl"
)

cat("===== RESULTS COMPARISON =====\n\n")

lme4_res <- results_lme4$all$results[1:10, ]
limma_res <- results_limma$all$results[1:10, ]

cat("LME4 Results (first 10):\n")
print(lme4_res)
cat("\n")

cat("LIMMA Results (first 10):\n")
print(limma_res)
cat("\n")

# Check sign of effects
cat("Effect size signs:\n")
cat("LME4 FU=1 (first 5 analytes): ", lme4_res$EFFECT_SIZE[lme4_res$FU == 1][1:5], "\n")
cat("LIMMA FU=1 (first 5 analytes): ", limma_res$EFFECT_SIZE[limma_res$FU == 1][1:5], "\n")
