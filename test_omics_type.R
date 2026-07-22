source("validation_helpers.R")

conditions <- character()
result <- withCallingHandlers(
  .validate_omics_type("other"),
  message = function(condition) {
    conditions <<- c(conditions, conditionMessage(condition))
    invokeRestart("muffleMessage")
  },
  warning = function(condition) {
    conditions <<- c(conditions, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)

stopifnot(is.null(result), length(conditions) == 0L)

invalid_error <- tryCatch(
  .validate_omics_type("unsupported"),
  error = identity
)

stopifnot(
  inherits(invalid_error, "error"),
  grepl("other", conditionMessage(invalid_error), fixed = TRUE)
)

cat("omics_type validation tests passed\n")
