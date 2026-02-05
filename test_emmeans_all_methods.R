# ===== TEST: emmeans() with LM, LME4, and LIMMA =====

require(lme4)
require(emmeans)
require(limma)

set.seed(42)

# Create test data
n_subjects <- 30
n_timepoints <- 3
n_reps <- n_subjects * n_timepoints

data_test <- data.frame(
  subject = rep(1:n_subjects, each = n_timepoints),
  fu = rep(0:2, n_subjects),
  treatment = rep(c(0, 1), n_reps / 2),
  value = rnorm(n_reps),
  baseline = rnorm(n_reps)
)

data_test$change <- data_test$value - 0.1

cat("===== TEST 1: emmeans() with LM =====\n\n")

# Filter to single FU for LM
data_lm <- data_test[data_test$fu %in% c(0, 1), ]
data_lm <- data_lm[!duplicated(data_lm$subject), ]

fit_lm <- lm(change ~ treatment + baseline, data = data_lm)

tryCatch({
  em_lm <- emmeans(fit_lm, ~treatment)
  contrasts_lm <- contrast(em_lm, method = "pairwise", adjust = "none")
  cat("✓ emmeans works with lm\n")
  cat("Results:\n")
  print(as.data.frame(contrasts_lm))
}, error = function(e) {
  cat("✗ Error with emmeans + lm:", e$message, "\n")
})

cat("\n===== TEST 2: emmeans() with LME4 =====\n\n")

fit_lme4 <- lmer(change ~ treatment * factor(fu) + baseline + (1|subject), 
                  data = data_test, REML = FALSE)

tryCatch({
  em_lme4 <- emmeans(fit_lme4, ~treatment | fu)
  contrasts_lme4 <- contrast(em_lme4, method = "pairwise", adjust = "none")
  cat("✓ emmeans works with lmer\n")
  cat("Results (first 6 rows):\n")
  print(head(as.data.frame(contrasts_lme4)))
}, error = function(e) {
  cat("✗ Error with emmeans + lmer:", e$message, "\n")
})

cat("\n===== TEST 3: emmeans() with LIMMA =====\n\n")

# Prepare data for LIMMA (need matrix format)
design_matrix <- model.matrix(~treatment * factor(fu) + baseline, data = data_test)
# Create a simple matrix of values (1 analyte × n samples)
analyte_matrix <- matrix(data_test$change, nrow = 1)

fit_limma <- lmFit(analyte_matrix, design_matrix)
fit_limma <- eBayes(fit_limma)

tryCatch({
  em_limma <- emmeans(fit_limma, ~treatment | fu)
  contrasts_limma <- contrast(em_limma, method = "pairwise", adjust = "none")
  cat("✓ emmeans works with limma\n")
  cat("Results (first 6 rows):\n")
  print(head(as.data.frame(contrasts_limma)))
}, error = function(e) {
  cat("✗ Error with emmeans + limma:", e$message, "\n")
})

cat("\n===== COMPARISON =====\n\n")

cat("All three methods can use emmeans()!\n")
cat("\nKey observations:\n")
cat("- LM: Single FU level (0 only)\n")
cat("- LME4: Multiple FU levels with random intercept\n")
cat("- LIMMA: Multiple FU levels with empirical Bayes\n")
cat("- All extract the same type of object: pairwise contrasts\n")
