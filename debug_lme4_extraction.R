# Debug LME4 extraction issue

source("main.R")
require(lme4)

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
omics <- omics[1:3, ]

cat("===== DEBUG LME4 EXTRACTION =====\n\n")

# Simulate what happens in the function
omics_sample_ids <- colnames(omics)[-which(colnames(omics) == "ANALYTE_NAME")]
shared_samples <- intersect(pheno$SAMPLE_ID, omics_sample_ids)
pheno_merged <- pheno[pheno$SAMPLE_ID %in% shared_samples, ]

fu_levels <- sort(unique(pheno_merged$FU))
cat("FU levels:", fu_levels, "\n\n")

baseline_subject_ids <- pheno$SUBJECT_ID[pheno$FU == 0]
omics_baseline_matrix <- as.matrix(omics[, -1])
# Only set colnames if dimensions match
if (ncol(omics_baseline_matrix) == length(omics_sample_ids)) {
  colnames(omics_baseline_matrix) <- omics_sample_ids
}

model_data <- data.frame(pheno_merged)
model_data$analyte <- NA_real_
model_data$analyte_baseline <- NA_real_

formula_str <- "analyte ~ CONTROL_STATUS * factor(FU) + agebl + (1|SUBJECT_ID)"

# Test with first analyte
analyte_name <- omics$ANALYTE_NAME[1]
fu_values <- as.numeric(omics[1, shared_samples])

baseline_idx <- match(pheno_merged$SUBJECT_ID, baseline_subject_ids)
baseline_col_idx <- match(pheno$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))
baseline_vals <- omics_baseline_matrix[1, baseline_col_idx]

analyte_change <- fu_values - baseline_vals
model_data$analyte <- analyte_change
model_data$analyte_baseline <- baseline_vals

cat("Fitting model...\n")
fit <- lmer(as.formula(formula_str), data = model_data, REML = FALSE)
cat("Model fitted successfully\n\n")

cat("Getting summary...\n")
fit_summary <- summary(fit)
cat("Summary obtained\n\n")

cat("Coefficient table:\n")
print(fit_summary$coefficients)
cat("\n")

cat("Row names of coefficient table:\n")
print(rownames(fit_summary$coefficients))
cat("\n")

cat("Checking for coefficient names:\n")
for (fu_level in fu_levels) {
  if (fu_level == fu_levels[1]) {
    coef_name <- "CONTROL_STATUS"
  } else {
    coef_name <- paste0("CONTROL_STATUS:factor(FU)", fu_level)
  }
  cat("Looking for:", coef_name, "... ")
  if (coef_name %in% rownames(fit_summary$coefficients)) {
    cat("FOUND\n")
  } else {
    cat("NOT FOUND\n")
  }
}
