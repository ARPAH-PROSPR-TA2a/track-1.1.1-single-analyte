# Test: Extract coefficients directly from LME4 (like LIMMA does) instead of using emmeans

require(lme4)
require(emmeans)

# Simulated data
set.seed(42)
n_subjects <- 20
n_timepoints <- 3
n_reps <- n_subjects * n_timepoints

data_test <- data.frame(
  subject = rep(1:n_subjects, each = n_timepoints),
  fu = rep(0:2, n_subjects),
  treatment = rep(c(0, 1), n_reps / 2),
  value = rnorm(n_reps),
  baseline = rnorm(n_reps)
)

data_test$change <- data_test$value - 0.1  # Simple difference

cat("===== COMPARING COEFFICIENT EXTRACTION METHODS =====\n\n")

# Model
formula_str <- "change ~ treatment * factor(fu) + baseline + (1|subject)"
fit <- lmer(as.formula(formula_str), data = data_test, REML = FALSE)

cat("Model formula: ", formula_str, "\n\n")

# Method 1: emmeans pairwise contrasts
cat("===== METHOD 1: emmeans() with pairwise contrasts =====\n")
em <- emmeans(fit, ~treatment | fu)
contrasts_result <- contrast(em, method = "pairwise", adjust = "none")
print(as.data.frame(contrasts_result))

# Method 2: Direct coefficient extraction (like LIMMA)
cat("\n===== METHOD 2: Direct coefficient extraction =====\n")
fixed_effects <- fixef(fit)
cat("Fixed effects:\n")
print(fixed_effects)
cat("\n")

cat("Extracting CONTROL_STATUS (treatment) coefficient at each FU:\n")
for (fu_val in c(0, 1, 2)) {
  fu_factor <- paste0("factor(fu)", fu_val)
  
  if (fu_val == 0) {
    # First level: just treatment coefficient
    treatment_effect <- fixed_effects["treatment"]
    cat(sprintf("FU=%d: treatment = %.6f (direct)\n", fu_val, treatment_effect))
  } else {
    # Subsequent levels: treatment + interaction
    treatment_effect <- fixed_effects["treatment"]
    interaction_name <- paste0("treatment:factor(fu)", fu_val)
    if (interaction_name %in% names(fixed_effects)) {
      combined_effect <- treatment_effect + fixed_effects[interaction_name]
      cat(sprintf("FU=%d: treatment + interaction = %.6f + %.6f = %.6f\n", 
                  fu_val, treatment_effect, fixed_effects[interaction_name], combined_effect))
    }
  }
}

cat("\n===== COMPARISON =====\n")
cat("emmeans pairwise contrasts give adjusted marginal means differences\n")
cat("Direct coefficients give model regression coefficients\n")
cat("These can differ in magnitude/interpretation!\n")
