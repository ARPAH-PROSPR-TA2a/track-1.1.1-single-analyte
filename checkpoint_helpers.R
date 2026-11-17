# Checkpoints contain model results before labels and multiple-testing correction.
.read_analysis_checkpoint <- function(batch_file, analyte_names, fu_levels,
                                      arm_levels, checkpoint_state = NULL) {
  if (is.null(batch_file) || !file.exists(batch_file)) return(NULL)

  invalid <- function(reason) {
    stop("Cannot reuse checkpoint '", batch_file, "': ", reason,
         ". Use a fresh checkpoint_dir for this analysis.", call. = FALSE)
  }
  payload <- tryCatch(readRDS(batch_file), error = function(e) {
    invalid(paste("could not read the file -", conditionMessage(e)))
  })

  coefficient_columns <- c("ANALYTE_NAME", "COEFFICIENT", "N_OBS",
                           "EFFECT_SIZE", "SE", "P_VALUE")
  legacy_effect_columns <- c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE")
  effect_columns <- c("ANALYTE_NAME", "FU", "CONTRAST", "TREATMENT", "REFERENCE",
                      "EFFECT_SIZE", "SE", "P_VALUE", "STATUS")
  numeric_columns <- c("EFFECT_SIZE", "SE", "P_VALUE")

  has_columns <- function(df, columns) {
    is.data.frame(df) && identical(sort(names(df)), sort(columns))
  }
  has_numeric_columns <- function(df, columns) {
    all(vapply(df[columns], function(x) is.numeric(x) && is.null(dim(x)), logical(1)))
  }
  has_identifiers <- function(x) {
    is.atomic(x) && is.null(dim(x)) && !anyNA(x)
  }
  is_result_pair <- function(x) {
    is.list(x) && !is.data.frame(x) && length(x) == 2L &&
      setequal(names(x), c("coefficients", "treatment_effects"))
  }
  valid_coefficients <- function(df, legacy = FALSE) {
    columns <- coefficient_columns
    if (legacy && !"N_OBS" %in% names(df)) columns <- setdiff(columns, "N_OBS")
    has_columns(df, columns) && nrow(df) > 0L &&
      has_identifiers(df$ANALYTE_NAME) && is.character(df$COEFFICIENT) &&
      is.null(dim(df$COEFFICIENT)) &&
      !anyNA(df$COEFFICIENT) &&
      has_numeric_columns(df, setdiff(columns, c("ANALYTE_NAME", "COEFFICIENT")))
  }
  is_legacy_result <- function(x) {
    if (!is_result_pair(x)) return(FALSE)
    if (is.null(x$coefficients) && is.null(x$treatment_effects)) return(TRUE)
    effects <- x$treatment_effects
    valid_coefficients(x$coefficients, legacy = TRUE) &&
      has_columns(effects, legacy_effect_columns) && nrow(effects) > 0L &&
      has_identifiers(effects$ANALYTE_NAME) &&
      has_numeric_columns(effects, c("FU", numeric_columns))
  }

  # The first checkpoint format saved one result (or NULL) per analyte.
  legacy_list <- is.list(payload) && !is.data.frame(payload) &&
    is.null(names(payload)) && length(payload) == length(analyte_names) &&
    length(payload) > 0L && all(vapply(payload, function(x) {
      is.null(x) || (is_legacy_result(x) && !is.null(x$coefficients))
    }, logical(1)))
  if (is_legacy_result(payload) || legacy_list) {
    if (is.null(checkpoint_state) || !isTRUE(checkpoint_state$notified)) {
      directory <- if (is.null(checkpoint_state)) dirname(batch_file) else checkpoint_state$directory
      message("Older-format checkpoints detected in '", directory,
              "'. FAST will recompute affected batches and replace them as they complete.")
      if (!is.null(checkpoint_state)) checkpoint_state$notified <- TRUE
    }
    return(NULL)
  }

  if (!is_result_pair(payload)) invalid("unrecognized result structure")
  effects <- payload$treatment_effects
  if (!has_columns(effects, effect_columns) || nrow(effects) == 0L) {
    invalid("unrecognized treatment-effect columns or empty treatment effects")
  }
  text_columns <- c("CONTRAST", "TREATMENT", "REFERENCE", "STATUS")
  if (!has_identifiers(effects$ANALYTE_NAME) ||
      !all(vapply(effects[text_columns], function(x) {
        is.character(x) && is.null(dim(x))
      }, logical(1))) ||
      anyNA(effects[text_columns]) ||
      !has_numeric_columns(effects, c("FU", numeric_columns))) {
    invalid("invalid treatment-effect column types or identifiers")
  }
  if (!all(effects$STATUS %in% c("ok", "not_estimable", "fit_error"))) {
    invalid("unrecognized treatment-effect status")
  }

  active_arms <- setdiff(as.character(arm_levels), "0")
  expected_keys <- expand.grid(
    ANALYTE_NAME = as.character(analyte_names),
    FU = as.numeric(as.character(fu_levels)),
    CONTRAST = paste0(active_arms, "_vs_0"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  actual_keys <- data.frame(
    ANALYTE_NAME = as.character(effects$ANALYTE_NAME),
    FU = as.numeric(effects$FU), CONTRAST = effects$CONTRAST,
    stringsAsFactors = FALSE
  )
  if (anyNA(actual_keys) || anyDuplicated(actual_keys) ||
      nrow(actual_keys) != nrow(expected_keys)) {
    invalid("missing or duplicated analyte, follow-up, or contrast rows")
  }
  actual_keys <- actual_keys[do.call(order, actual_keys), , drop = FALSE]
  expected_keys <- expected_keys[do.call(order, expected_keys), , drop = FALSE]
  rownames(actual_keys) <- rownames(expected_keys) <- NULL
  if (!identical(actual_keys, expected_keys)) {
    invalid("analyte, follow-up, or contrast rows do not match this batch")
  }
  if (!all(effects$TREATMENT %in% active_arms) ||
      !all(effects$REFERENCE == "0") ||
      !all(effects$CONTRAST == paste0(effects$TREATMENT, "_vs_0"))) {
    invalid("treatment or reference codes do not match the contrasts")
  }

  coefficients <- payload$coefficients
  if (is.null(coefficients)) {
    if (any(effects$STATUS == "ok")) {
      invalid("successful treatment effects have no coefficients")
    }
  } else {
    if (!valid_coefficients(coefficients) ||
        !all(as.character(coefficients$ANALYTE_NAME) %in% as.character(analyte_names)) ||
        anyDuplicated(coefficients[c("ANALYTE_NAME", "COEFFICIENT")])) {
      invalid("invalid coefficient columns or analyte rows")
    }
    successful_analytes <- effects$ANALYTE_NAME[effects$STATUS == "ok"]
    if (!all(successful_analytes %in% coefficients$ANALYTE_NAME)) {
      invalid("successful treatment effects are missing coefficients")
    }
  }
  payload
}


.save_analysis_checkpoint <- function(batch_payload, batch_file) {
  if (is.null(batch_file)) return(invisible(NULL))

  temporary_file <- paste0(batch_file, ".tmp")
  on.exit(unlink(temporary_file), add = TRUE)
  tryCatch(saveRDS(batch_payload, temporary_file), error = function(e) {
    stop("Could not write checkpoint '", batch_file, "': ", conditionMessage(e),
         call. = FALSE)
  })
  if (!file.rename(temporary_file, batch_file)) {
    stop("Could not replace checkpoint '", batch_file,
         "'. The previous checkpoint was retained.", call. = FALSE)
  }
  invisible(NULL)
}
