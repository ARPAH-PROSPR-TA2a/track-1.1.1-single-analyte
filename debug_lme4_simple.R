# Simple LME4 extraction debug

source("main.R")
require(lme4)

# Just test the extraction logic directly
cat("===== TESTING COEFFICIENT EXTRACTION =====\n\n")

# Create simple test data
set.seed(42)
test_data <- data.frame(
  subject = rep(1:10, each=3),
  fu = rep(0:2, 10),
  treatment = rep(c(0,1), 15),
  change = rnorm(30),
  baseline = rnorm(30),
  age = rep(rnorm(10), each=3)
)

fu_levels <- sort(unique(test_data$fu))
cat("FU levels:", fu_levels, "\n\n")

# Fit model
fit <- lmer(change ~ treatment * factor(fu) + baseline + age + (1|subject), 
            data = test_data, REML = FALSE)

cat("Model fitted. Getting summary...\n")
fit_summary <- summary(fit)

cat("\nCoefficient table:\n")
coef_table <- fit_summary$coefficients
print(coef_table)

cat("\n\nRow names:", rownames(coef_table), "\n")

cat("\n\nExtracting for each FU level:\n")
for (fu_level in fu_levels) {
  if (fu_level == fu_levels[1]) {
    coef_name <- "treatment"
  } else {
    coef_name <- paste0("treatment:factor(fu)", fu_level)
  }
  
  cat("\nFU=", fu_level, " -> Looking for: '", coef_name, "'\n", sep="")
  
  if (coef_name %in% rownames(coef_table)) {
    effect_size <- coef_table[coef_name, "Estimate"]
    se <- coef_table[coef_name, "Std. Error"]
    p_value <- coef_table[coef_name, "Pr(>|t|)"]
    cat("  Found! ES=", round(effect_size, 4), ", SE=", round(se, 4), ", p=", round(p_value, 4), "\n")
  } else {
    cat("  NOT FOUND\n")
  }
}
