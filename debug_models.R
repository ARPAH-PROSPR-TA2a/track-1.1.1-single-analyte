# ===== DEBUG: Check model coefficients =====

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

analyte_names <- rownames(omics_raw)
sample_names <- colnames(omics_raw)
omics <- as.data.frame(omics_raw)
colnames(omics) <- sample_names
omics <- cbind(ANALYTE_NAME = analyte_names, omics, stringsAsFactors = FALSE)
omics <- omics[1:3, ]  # Just 3 analytes for debugging

additional_covariates <- c("agebl")

cat("===== DEBUG: Single Analyte Comparison =====\n\n")

# Get one analyte to debug
omics_sample_ids <- colnames(omics)[-which(colnames(omics) == "ANALYTE_NAME")]
shared_samples <- intersect(pheno$SAMPLE_ID, omics_sample_ids)
pheno_merged <- pheno[pheno$SAMPLE_ID %in% shared_samples, ]

omics_baseline_matrix <- as.matrix(omics[, -1])  # Remove ANALYTE_NAME column
colnames(omics_baseline_matrix) <- omics_sample_ids

fu_values <- as.numeric(omics_baseline_matrix[1, match(shared_samples, colnames(omics_baseline_matrix))])

baseline_sample_ids <- pheno$SAMPLE_ID[pheno$FU == 0]
baseline_indices <- match(baseline_sample_ids, colnames(omics_baseline_matrix))
baseline_vals <- omics_baseline_matrix[1, baseline_indices]

cat("First Analyte: ", omics$ANALYTE_NAME[1], "\n")
cat("Number of samples: ", length(shared_samples), "\n")
cat("First 5 FU values: ", head(fu_values, 5), "\n")
cat("First 5 baseline values: ", head(baseline_vals, 5), "\n\n")

# === LME4 ANALYSIS ===
cat("===== LME4 ANALYSIS =====\n\n")

model_data <- data.frame(pheno_merged)
model_data$analyte <- NA_real_
model_data$analyte_baseline <- NA_real_

fu_levels <- sort(unique(pheno_merged$FU))
baseline_subject_ids <- pheno$SUBJECT_ID[pheno$FU == 0]
baseline_idx <- match(pheno_merged$SUBJECT_ID, baseline_subject_ids)
baseline_col_idx <- match(pheno$SAMPLE_ID[match(pheno_merged$SUBJECT_ID[which(!is.na(baseline_idx))], baseline_subject_ids)], 
                          colnames(omics_baseline_matrix))
# Simpler: for each subject in pheno_merged, get their baseline value
baseline_vals_all <- rep(NA_real_, nrow(pheno_merged))
for (i in 1:nrow(pheno_merged)) {
  subj <- pheno_merged$SUBJECT_ID[i]
  baseline_samp <- pheno$SAMPLE_ID[pheno$SUBJECT_ID == subj & pheno$FU == 0]
  if (length(baseline_samp) > 0) {
    baseline_col_idx <- match(baseline_samp[1], colnames(omics_baseline_matrix))
    if (!is.na(baseline_col_idx)) {
      baseline_vals_all[i] <- omics_baseline_matrix[1, baseline_col_idx]
    }
  }
}
analyte_change <- fu_values - baseline_vals_all

model_data$analyte <- analyte_change
model_data$analyte_baseline <- baseline_vals_all

formula_str <- "analyte ~ CONTROL_STATUS * factor(FU) + agebl + (1|SUBJECT_ID)"
cat("Formula: ", formula_str, "\n")
cat("Data summary:\n")
cat("  Rows:", nrow(model_data), "\n")
cat("  Change score range:", range(model_data$analyte), "\n")
cat("  Baseline range:", range(model_data$analyte_baseline), "\n\n")

fit_lme4 <- lmer(as.formula(formula_str), data = model_data, REML = FALSE)
cat("LME4 Fit Summary:\n")
print(summary(fit_lme4))
cat("\n")

em <- emmeans(fit_lme4, ~CONTROL_STATUS | FU)
cat("Emmeans object:\n")
print(em)
cat("\n")

contrasts_result <- contrast(em, method = "pairwise", adjust = "none")
cat("Contrasts result:\n")
print(as.data.frame(contrasts_result))
cat("\n")

# === LIMMA ANALYSIS ===
cat("===== LIMMA ANALYSIS =====\n\n")

# Prepare data for LIMMA
pheno_limma <- pheno_merged
pheno_limma$FU_factor <- factor(pheno_limma$FU)
pheno_limma$agebl_cov <- pheno_limma$agebl

design <- model.matrix(~ CONTROL_STATUS * FU_factor + agebl_cov, data = pheno_limma)
cat("Design matrix columns:\n")
cat(colnames(design), "\n\n")

# Create matrix with just first analyte
analyte_matrix <- matrix(analyte_change, nrow = 1)

cat("Analyte matrix shape: ", nrow(analyte_matrix), "x", ncol(analyte_matrix), "\n\n")

cor <- duplicateCorrelation(analyte_matrix, design, block = pheno_limma$SUBJECT_ID)
cat("Consensus correlation: ", cor$consensus.correlation, "\n\n")

fit_limma <- lmFit(analyte_matrix, design,
                   block = pheno_limma$SUBJECT_ID,
                   correlation = cor$consensus.correlation)
fit_limma <- eBayes(fit_limma)

cat("LIMMA Coefficients:\n")
cat("Column names: ", colnames(design), "\n")
print(fit_limma$coefficients)
cat("\n")

cat("LIMMA Coefficients (first row only):\n")
print(fit_limma$coefficients[1, ])
cat("\n")

cat("LIMMA P-values:\n")
print(fit_limma$p.value[1, ])
cat("\n")

# Extract CONTROL_STATUS coefficient for FU=1
coef_idx_cs <- which(colnames(design) == "CONTROL_STATUS")
cat("CONTROL_STATUS coefficient (FU=1): ", fit_limma$coefficients[1, coef_idx_cs], "\n")
cat("CONTROL_STATUS p-value: ", fit_limma$p.value[1, coef_idx_cs], "\n\n")

# Extract combined effect for FU=2
coef_int <- "CONTROL_STATUS:FU_factor2"
coef_idx_base <- which(colnames(design) == "CONTROL_STATUS")
coef_idx_int <- which(colnames(design) == coef_int)

if (length(coef_idx_int) > 0) {
  combined_effect <- fit_limma$coefficients[1, coef_idx_base] + fit_limma$coefficients[1, coef_idx_int]
  cat("CONTROL_STATUS effect at FU=2 (combined): ", combined_effect, "\n")
  cat("  Base: ", fit_limma$coefficients[1, coef_idx_base], "\n")
  cat("  Interaction: ", fit_limma$coefficients[1, coef_idx_int], "\n")
}
