# Run from the repository root: Rscript --vanilla test_dnam_filtered_membership.R
source("main.R")
source("tests/synthetic_data.R")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS: ", label, "\n", sep = "")
}

local({
  fixture <- .make_synthetic_test_data(n_analytes = 6L)
  pheno <- fixture$pheno
  omics <- fixture$omics
  probes <- omics$ANALYTE_NAME
  # Included probe 1 is estimable; probe 2 lacks control observations;
  # probe 3 has only female observations, so its all-sex fit fails.
  omics[2L, pheno$SAMPLE_ID[pheno$TREATMENT_GROUP == "0"]] <- NA_real_
  omics[3L, pheno$SAMPLE_ID[pheno$FEMALE == "0"]] <- NA_real_
  directory <- tempfile("track111_membership_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  for (multiarm in c(FALSE, TRUE)) {
    input <- pheno
    labels <- NULL
    if (multiarm) {
      active <- unique(input$SUBJECT_ID[input$TREATMENT_GROUP == "1"])
      second_arm <- active[rep(c(FALSE, FALSE, TRUE, TRUE), length.out = length(active))]
      codes <- as.character(input$TREATMENT_GROUP)
      codes[input$SUBJECT_ID %in% second_arm] <- "2"
      input$TREATMENT_GROUP <- factor(codes, levels = 0:2)
      labels <- c(`0` = "Control", `1` = "Diet", `2` = "Exercise")
    }
    for (mixed in c(FALSE, TRUE)) {
      p <- if (mixed) input else input[input$FU != "2", ]
      name <- paste(if (multiarm) "multiarm" else "binary", if (mixed) "LME4" else "LM")
      checkpoint <- file.path(directory, gsub(" ", "_", name))
      run <- function(type, checkpoint_dir = NULL) quiet(FAST_omics_WAS(
        p, omics, omics_type = type, additional_covariates = "agebl",
        treatment_labels = labels, n_cores = 1L, checkpoint_dir = checkpoint_dir,
        checkpoint_batch_size = 3L
      ))
      reference <- run("other")
      fresh <- run("DNAm", checkpoint)
      without_dnam_columns <- fresh

      for (analysis in c("analysis_change", "analysis_level")) {
        for (stratum in c("all", "male", "female")) {
          te <- fresh[[analysis]][[stratum]]$treatment_effects
          check(identical(te$IN_FILTERED_SET, te$ANALYTE_NAME %in% fixture$filtered_probes),
                paste(name, analysis, stratum, "has logical membership on every contrast row"))
          check(!"IN_FILTERED_SET" %in% names(fresh[[analysis]][[stratum]]$coefficients) &&
                  !"IN_FILTERED_SET" %in% names(reference[[analysis]][[stratum]]$treatment_effects),
                paste(name, analysis, stratum, "adds membership only to DNAm treatment effects"))

          expected_bh <- rep(NA_real_, nrow(te))
          for (idx in split(seq_len(nrow(te)), interaction(te$FU, te$CONTRAST, drop = TRUE))) {
            included <- idx[te$ANALYTE_NAME[idx] %in% fixture$filtered_probes]
            expected_bh[included] <- p.adjust(te$P_VALUE[included], "BH")
          }
          check(isTRUE(all.equal(te$BH_P_VALUE_FILTERED, expected_bh)),
                paste(name, analysis, stratum, "preserves filtered BH calculation"))

          for (slot in c("coefficients", "treatment_effects")) {
            table <- without_dnam_columns[[analysis]][[stratum]][[slot]]
            table$BH_P_VALUE_FILTERED <- NULL
            table$IN_FILTERED_SET <- NULL
            without_dnam_columns[[analysis]][[stratum]][[slot]] <- table
          }
        }
        te <- fresh[[analysis]]$all$treatment_effects
        valid <- te$ANALYTE_NAME == probes[1L]
        absent <- te$ANALYTE_NAME == probes[2L]
        failed <- te$ANALYTE_NAME == probes[3L]
        excluded <- te$ANALYTE_NAME %in% probes[4:6]
        check(all(te$STATUS[valid] == "ok") && all(te$IN_FILTERED_SET[valid]) &&
                all(is.finite(te$BH_P_VALUE_FILTERED[valid])) &&
                all(te$STATUS[absent] == "not_estimable") && all(te$IN_FILTERED_SET[absent]) &&
                all(is.na(te$BH_P_VALUE_FILTERED[absent])) &&
                all(te$STATUS[failed] == "fit_error") && all(te$IN_FILTERED_SET[failed]) &&
                all(is.na(te$BH_P_VALUE_FILTERED[failed])) &&
                all(!te$IN_FILTERED_SET[excluded]) && all(is.na(te$BH_P_VALUE_FILTERED[excluded])),
              paste(name, analysis, "distinguishes membership from unavailable inference"))
      }
      check(isTRUE(all.equal(without_dnam_columns, reference, tolerance = 1e-8)),
            paste(name, "preserves model outputs and full-set BH"))

      batches <- list.files(checkpoint, "^batch_[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE)
      check(length(batches) > 1L && all(vapply(batches, function(path) {
        !"IN_FILTERED_SET" %in% names(readRDS(path)$treatment_effects)
      }, logical(1))), paste(name, "keeps membership outside checkpoint batches"))
      modification_times <- file.info(batches)$mtime
      resumed <- run("DNAm", checkpoint)
      check(isTRUE(all.equal(fresh, resumed, tolerance = 1e-8)) &&
              identical(modification_times, file.info(batches)$mtime),
            paste(name, "reconstructs membership on resume without rewriting batches"))
    }
  }
})

cat("All DNAm filtered-membership tests passed.\n")
