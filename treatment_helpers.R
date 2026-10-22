# Codes identify arms throughout the models; labels describe them in the outputs.
.prepare_treatments <- function(pheno, treatment_labels = NULL) {
  codes <- as.character(pheno$TREATMENT_GROUP)
  numbers <- suppressWarnings(as.integer(codes))
  if (anyNA(numbers) || any(numbers < 0L) || any(as.character(numbers) != codes)) {
    stop("TREATMENT_GROUP must contain non-missing non-negative integer codes: 0 = control; 1, 2, ... = active treatments.")
  }
  arms <- as.character(sort(unique(numbers)))
  if (!"0" %in% arms || length(arms) < 2L) {
    stop("pheno must contain control (TREATMENT_GROUP == 0) and at least one active treatment.")
  }

  # Validate before dropping duplicate visits or incomplete covariates.
  by_subject <- split(codes, pheno$SUBJECT_ID, drop = TRUE)
  if (any(vapply(by_subject, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("TREATMENT_GROUP must be constant within SUBJECT_ID.")
  }
  if (is.null(treatment_labels)) {
    if (length(arms) > 2L) {
      stop('Multiple active treatments require treatment_labels covering every code, including control 0; e.g., c(`0` = "Control", `1` = "Diet", `2` = "Exercise").')
    }
    treatment_labels <- setNames(c("Control", "Treatment"), arms)
  }
  if (!is.character(treatment_labels) || !is.null(dim(treatment_labels)) ||
      is.null(names(treatment_labels)) || anyNA(names(treatment_labels)) ||
      anyDuplicated(names(treatment_labels)) || !setequal(names(treatment_labels), arms)) {
    stop("treatment_labels must be a named character vector with unique code names matching every observed TREATMENT_GROUP code exactly, including 0.")
  }
  if (anyNA(treatment_labels) || any(!nzchar(trimws(treatment_labels))) ||
      anyDuplicated(treatment_labels)) {
    stop("treatment_labels must contain unique, non-missing, non-blank labels.")
  }
  treatment_labels <- treatment_labels[arms]
  active_arms <- setdiff(arms, "0")
  treatments <- data.frame(
    ARM = active_arms, TREATMENT = unname(treatment_labels[active_arms]),
    REFERENCE = unname(treatment_labels["0"]),
    CONTRAST = paste0(active_arms, "_vs_0"), stringsAsFactors = FALSE
  )
  pheno$TREATMENT_GROUP <- factor(codes, levels = arms)
  list(pheno = pheno, treatments = treatments, treatment_labels = treatment_labels)
}

.treatment_metadata <- function(treatments, treatment_labels) {
  list(
    schema_version = "2.0", treatment_labels = treatment_labels,
    treatments = treatments,
    p_adjust = "BH within analysis, stratum, FU and contrast"
  )
}

# Keep requested comparisons visible even when an arm has no usable observations.
.empty_treatment_effects <- function(analyte_name, fu_levels, arm_levels,
                                     status = "not_estimable") {
  grid <- expand.grid(FU = as.integer(as.character(fu_levels)),
                      TREATMENT = setdiff(arm_levels, "0"),
                      stringsAsFactors = FALSE)
  data.frame(
    ANALYTE_NAME = analyte_name, FU = grid$FU,
    CONTRAST = paste0(grid$TREATMENT, "_vs_0"),
    TREATMENT = grid$TREATMENT, REFERENCE = "0",
    EFFECT_SIZE = NA_real_, SE = NA_real_, P_VALUE = NA_real_,
    STATUS = status, stringsAsFactors = FALSE
  )
}

# The ordinary single-FU model already has the requested contrasts as coefficients.
# Use emmeans for mixed models and for unusual/aliased LM parameterizations.
.extract_treatment_effects <- function(fit, analyte_name, fu_levels, arm_levels,
                                       mixed = FALSE, coefficient_table = NULL) {
  out <- .empty_treatment_effects(analyte_name, fu_levels, arm_levels)
  direct_coefficients <- !mixed &&
    identical(fit$xlevels$TREATMENT_GROUP[1L], "0") &&
    identical(fit$contrasts$TREATMENT_GROUP, "contr.treatment") &&
    !anyNA(stats::coef(fit))
  if (direct_coefficients) {
    if (is.null(coefficient_table)) coefficient_table <- summary(fit)$coefficients
    idx <- match(paste0("TREATMENT_GROUP", out$TREATMENT), rownames(coefficient_table))
    out$EFFECT_SIZE <- coefficient_table[idx, "Estimate"]
    out$SE <- coefficient_table[idx, "Std. Error"]
    out$P_VALUE <- coefficient_table[idx, "Pr(>|t|)"]
  } else {
    emm <- emmeans::emmeans(fit, if (mixed) ~ TREATMENT_GROUP | FU else ~ TREATMENT_GROUP,
                            nesting = NULL)
    observed_arms <- as.character(unique(as.data.frame(emm)$TREATMENT_GROUP))
    if (!"0" %in% observed_arms) return(out)
    active_arms <- intersect(setdiff(arm_levels, "0"), observed_arms)
    if (!length(active_arms)) return(out)
    contrasts <- lapply(active_arms, function(arm) {
      as.numeric(observed_arms == arm) - as.numeric(observed_arms == "0")
    })
    names(contrasts) <- paste0(active_arms, "_vs_0")
    estimates <- as.data.frame(emmeans::contrast(emm, method = contrasts, adjust = "none"))
    if (!mixed) estimates$FU <- fu_levels[[1L]]
    idx <- match(paste(out$CONTRAST, out$FU), paste(estimates$contrast, estimates$FU))
    out$EFFECT_SIZE <- estimates$estimate[idx]
    out$SE <- estimates$SE[idx]
    out$P_VALUE <- estimates$p.value[idx]
  }
  ok <- is.finite(out$EFFECT_SIZE) & is.finite(out$SE) & is.finite(out$P_VALUE)
  out$STATUS[ok] <- "ok"
  out[!ok, c("EFFECT_SIZE", "SE", "P_VALUE")] <- NA_real_
  out
}
