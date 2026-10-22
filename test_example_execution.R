# Run from the repository root: Rscript --vanilla test_example_execution.R
# Exercise example loading and workspace behavior; statistical tests live elsewhere.
local({
  example <- normalizePath("Examples/run_multiarm_example.R", mustWork = TRUE)
  repository <- dirname(dirname(example))
  original_wd <- getwd()
  original_globals <- ls(.GlobalEnv, all.names = TRUE)
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  original_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  temporary <- tempfile("track111-example-execution-")
  dir.create(temporary)
  temporary <- normalizePath(temporary, mustWork = TRUE)
  on.exit({
    setwd(original_wd)
    if (had_seed) {
      assign(".Random.seed", original_seed, .GlobalEnv)
    } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    unlink(temporary, recursive = TRUE)
  }, add = TRUE)

  check <- function(condition, label) {
    if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
    cat("PASS: ", label, "\n", sep = "")
  }
  expected_objects <- c("pheno", "omics", "results", "reports",
                         "effects_diet", "effects_diet_exercise")
  helpers <- c("FAST_omics_WAS", "FAST_omics_WAS_reports", ".log_111",
               ".validate_pheno", ".validate_omics", ".prepare_treatments",
               ".perform_lm_analysis", ".generate_reports")

  setwd(temporary)
  environments <- lapply(c(FALSE, TRUE), function(change_directory) {
    workspace <- new.env(parent = .GlobalEnv)
    invisible(capture.output(suppressWarnings(suppressMessages(
      source(example, local = workspace, chdir = change_directory)
    ))))
    mode <- paste0("source(chdir = ", change_directory, ")")
    check(identical(getwd(), temporary), paste(mode, "preserves the caller's working directory"))
    check(all(vapply(expected_objects, exists, logical(1), envir = workspace,
                     inherits = FALSE)) &&
            nrow(workspace$effects_diet) > 0L && nrow(workspace$effects_diet_exercise) > 0L,
          paste(mode, "leaves inputs and both treatment effects inspectable"))
    check(identical(workspace$.track111_root, repository) &&
            all(vapply(helpers, function(name) {
              exists(name, workspace, inherits = FALSE) &&
                identical(environment(get(name, workspace)), workspace)
            }, logical(1))),
          paste(mode, "loads pipeline helpers into the requested environment"))
    workspace
  })
  check(isTRUE(all.equal(environments[[1L]]$results, environments[[2L]]$results,
                        tolerance = 1e-8)) &&
          isTRUE(all.equal(environments[[1L]]$reports, environments[[2L]]$reports,
                           tolerance = 1e-8)),
        "both source modes produce equivalent analysis and report objects")
  check(!length(setdiff(ls(.GlobalEnv, all.names = TRUE),
                         c(original_globals, ".Random.seed"))),
        "example inputs, outputs, and pipeline helpers do not leak into the global workspace")

  # Simulate RStudio's saved-editor path in an isolated R process. This checks
  # section-by-section evaluation without claiming to operate an actual RStudio UI.
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    script <- file.path(temporary, "simulated_editor.R")
    writeLines(c(
      "local({",
      "  example <- normalizePath(commandArgs(TRUE)[1L], mustWork = TRUE)",
      "  ns <- asNamespace('rstudioapi')",
      "  available <- get('isAvailable', ns)",
      "  editor <- get('getSourceEditorContext', ns)",
      "  on.exit({",
      "    assignInNamespace('isAvailable', available, ns = 'rstudioapi')",
      "    assignInNamespace('getSourceEditorContext', editor, ns = 'rstudioapi')",
      "  }, add = TRUE)",
      "  assignInNamespace('isAvailable', function(...) TRUE, ns = 'rstudioapi')",
      "  assignInNamespace('getSourceEditorContext', function(...) list(path = example),",
      "                    ns = 'rstudioapi')",
      "  setwd(tempdir())",
      "  caller_wd <- getwd()",
      "  workspace <- new.env(parent = globalenv())",
      "  for (expr in parse(example)) {",
      "    if (is.call(expr) && identical(expr[[1L]], as.name('<-')) &&",
      "        identical(expr[[2L]], as.name('results'))) break",
      "    eval(expr, envir = workspace)",
      "  }",
      "  stopifnot(identical(getwd(), caller_wd),",
      "            identical(workspace$.track111_root, dirname(dirname(example))),",
      "            exists('pheno', workspace, inherits = FALSE),",
      "            exists('omics', workspace, inherits = FALSE),",
      "            identical(environment(workspace$FAST_omics_WAS), workspace),",
      "            identical(environment(workspace$.validate_omics), workspace),",
      "            !exists('FAST_omics_WAS', globalenv(), inherits = FALSE))",
      "})"
    ), script)
    log <- file.path(temporary, "simulated_editor.log")
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(script), shQuote(example)),
                      stdout = log, stderr = log)
    if (status != 0L) cat(readLines(log), sep = "\n")
    check(status == 0L, "simulated saved RStudio editor supports setup and data sections")
  } else {
    cat("SKIP: simulated RStudio editor check (rstudioapi is not installed).\n")
  }
  cat("Example execution tests passed.\n")
})
