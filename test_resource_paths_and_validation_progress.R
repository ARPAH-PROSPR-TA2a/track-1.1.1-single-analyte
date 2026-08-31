source("main.R")

run_resource_and_progress_tests <- function() {
  original_wd <- getwd()
  original_root <- .track111_root
  original_progress_log <- getOption("track111.progress_log")
  on.exit(setwd(original_wd), add = TRUE)
  on.exit(assign(".track111_root", original_root, envir = .GlobalEnv), add = TRUE)
  on.exit(options(track111.progress_log = original_progress_log), add = TRUE)
  options(track111.progress_log = NULL)

  setwd(tempdir())
  probes <- .load_dnam_probe_lists()
  stopifnot(length(probes$full) > 0L, length(probes$filtered) > 0L)
  setwd(original_wd)

  assign(".track111_root", tempfile("missing-track111-"), envir = .GlobalEnv)
  missing_resource_error <- tryCatch(
    FAST_omics_WAS(
      pheno = data.frame(),
      omics = data.frame(),
      omics_type = "DNAm",
      n_cores = 1L
    ),
    error = conditionMessage
  )
  stopifnot(grepl(
    "Required DNAm probe-list file",
    missing_resource_error,
    fixed = TRUE
  ))
  assign(".track111_root", original_root, envir = .GlobalEnv)

  pheno_list <- list(
    all = data.frame(SAMPLE_ID = c("s1", "s2", "s3")),
    male = NULL,
    female = NULL
  )
  omics <- data.frame(
    ANALYTE_NAME = c(
      "variable",
      "constant",
      "partial_variable",
      "partial_constant",
      "all_na"
    ),
    s1 = c(1, 4, NA, NA, NA),
    s2 = c(2, 4, 1, 5, NA),
    s3 = c(3, 4, 2, 5, NA),
    check.names = FALSE
  )

  validate <- function(verbose) {
    warnings_seen <- character()
    result <- NULL
    messages_seen <- capture.output(
      result <- withCallingHandlers(
        .validate_omics(
          omics,
          pheno_list,
          verbose = verbose,
          progress_every = 2L
        ),
        warning = function(warning_condition) {
          warnings_seen <<- c(
            warnings_seen,
            conditionMessage(warning_condition)
          )
          invokeRestart("muffleWarning")
        }
      ),
      type = "message"
    )
    list(result = result, warnings = warnings_seen, messages = messages_seen)
  }

  quiet <- validate(FALSE)
  verbose <- validate(TRUE)
  quality_lines <- verbose$messages[
    startsWith(verbose$messages, "[1.1.1] omics quality checks")
  ]

  stopifnot(
    identical(quiet$result, verbose$result),
    identical(
      quiet$warnings,
      c(
        "3 analytes contain NA values",
        "2 analytes have near-zero variance"
      )
    ),
    identical(quiet$warnings, verbose$warnings),
    !any(startsWith(quiet$messages, "[1.1.1] omics quality checks")),
    length(quality_lines) == 4L,
    grepl("0/5", quality_lines[[1L]], fixed = TRUE),
    grepl("2/5", quality_lines[[2L]], fixed = TRUE),
    grepl("4/5", quality_lines[[3L]], fixed = TRUE),
    grepl("5/5 analytes complete", quality_lines[[4L]], fixed = TRUE),
    grepl(
      "3 with NA, 2 near-zero variance",
      quality_lines[[4L]],
      fixed = TRUE
    ),
    identical(verbose$result$all, omics)
  )

  cat("Resource path and validation progress tests passed.\n")
}

run_resource_and_progress_tests()
