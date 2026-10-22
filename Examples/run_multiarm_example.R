# Open this saved script in RStudio and click Source, or run each section in order.
# Running sections individually uses rstudioapi to locate the open script.
# From an R console, source("/path/to/Track1.1.1/Examples/run_multiarm_example.R").
# All data below are synthetic. Nothing is read from participant data files.
# The input data and results remain in your R workspace for inspection.

# 1. Load the pipeline ----
# Locate this script without changing your working directory.
example_dir <- local({
  source_frames <- Filter(function(frame) !is.null(frame$ofile), sys.frames())
  if (length(source_frames)) {
    frame <- tail(source_frames, 1L)[[1L]]
    return(if (isTRUE(frame$chdir)) getwd() else dirname(normalizePath(frame$ofile)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getSourceEditorContext()$path
  } else {
    script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  }
  if (length(script_path) != 1L || is.na(script_path) || !nzchar(script_path)) {
    stop("Could not locate this example. Open the saved script in RStudio and click Source, or source() it by its full path. Running sections individually requires rstudioapi.")
  }
  dirname(normalizePath(script_path, mustWork = TRUE))
})
source(file.path(example_dir, "..", "main.R"), local = TRUE)

# 2. Create synthetic data ----
set.seed(20261012)
n_by_arm <- c(`0` = 40L, `1` = 32L, `2` = 24L)
labels <- c(`0` = "Usual care", `1` = "Diet", `2` = "Diet + exercise")
arm <- rep(names(n_by_arm), n_by_arm)
n_subjects <- length(arm)
n_analytes <- 6L
followups <- 0:2

subjects <- data.frame(
  SUBJECT_ID = sprintf("synthetic_%03d", seq_len(n_subjects)),
  TREATMENT_GROUP = factor(arm, levels = names(n_by_arm)),
  FEMALE = factor(rep(0:1, length.out = n_subjects), levels = 0:1)
)
pheno <- subjects[rep(seq_len(n_subjects), each = length(followups)), ]
pheno$FU <- factor(rep(followups, times = n_subjects), levels = followups)
pheno$SAMPLE_ID <- paste0(pheno$SUBJECT_ID, "_FU", pheno$FU)
row.names(pheno) <- NULL

# Active arms have distinct effects against the shared usual-care control.
effects_by_arm <- list(
  `0` = rep(0, n_analytes),
  `1` = c(0.8, -0.5, 0.4, 0, 0, 0),
  `2` = c(-0.4, 0.7, 0.4, 0, 0, 0)
)
baseline <- matrix(rnorm(n_analytes * n_subjects), nrow = n_analytes)
subject_change <- matrix(rnorm(n_analytes * n_subjects, sd = 0.4),
                         nrow = n_analytes)
values <- matrix(NA_real_, nrow = n_analytes, ncol = nrow(pheno))
for (i in seq_len(n_subjects)) {
  arm_effect <- effects_by_arm[[arm[i]]]
  for (fu in followups) {
    j <- (i - 1L) * length(followups) + fu + 1L
    values[, j] <- if (fu == 0L) {
      baseline[, i]
    } else {
      baseline[, i] + subject_change[, i] + fu * (0.1 + arm_effect) +
        rnorm(n_analytes, sd = 0.5)
    }
  }
}
colnames(values) <- pheno$SAMPLE_ID
omics <- data.frame(
  ANALYTE_NAME = paste0("synthetic_analyte_", seq_len(n_analytes)),
  values,
  check.names = FALSE
)

# 3. Fit the models and generate reports ----
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "other",
  treatment_labels = labels,
  n_cores = 1
)
reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "other",
  treatment_labels = labels
)

# 4. Inspect the results ----
# One joint fit; two separate treatment-effect data frames for each stratum.
effects <- split(results$analysis_change$all$treatment_effects,
                 results$analysis_change$all$treatment_effects$TREATMENT)
effects_diet <- effects[["Diet"]]
effects_diet_exercise <- effects[["Diet + exercise"]]

cat("Synthetic example: 40 usual care, 32 diet, 24 diet + exercise; two follow-ups.\n")
print(results$metadata$treatments)
print(effects_diet)
print(effects_diet_exercise)
cat("Per-arm report counts:\n")
print(reports$arm_counts)
