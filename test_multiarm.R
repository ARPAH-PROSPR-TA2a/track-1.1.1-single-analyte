# Synthetic regression tests: Rscript --vanilla test_multiarm.R
# No trial data, network access, or persistent output is needed.
source("main.R")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS: ", label, "\n", sep = "")
}
near <- function(x, y, tolerance = 1e-7) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
}
expect_error <- function(expr, label, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(error, "error") &&
          (is.null(pattern) || grepl(pattern, conditionMessage(error), ignore.case = TRUE)),
        label)
}

make_synthetic <- function() {
  set.seed(11103)
  arm <- rep(c("control", "A", "B"), c(44L, 32L, 20L))
  n <- length(arm)
  subject <- rep(seq_len(n), each = 3L)
  visit <- rep(0:2, times = n)
  sex <- rep(0:1, length.out = n)
  age <- rnorm(n, 50, 8)
  pheno <- data.frame(
    SAMPLE_ID = sprintf("synthetic_%03d_FU%d", subject, visit),
    SUBJECT_ID = sprintf("synthetic_%03d", subject),
    FU = factor(visit), FEMALE = factor(sex[subject]),
    TREATMENT_GROUP = factor(match(arm[subject], c("control", "A", "B")) - 1L),
    AGE = age[subject], stringsAsFactors = FALSE
  )
  values <- t(vapply(seq_len(4L), function(j) {
    baseline <- rnorm(n, 5, 0.7)
    subject_change <- rnorm(n, 0, 0.7)
    tx <- ifelse(arm == "A", c(0.4, 0.1, 0.7, 0)[j],
                 ifelse(arm == "B", c(-0.2, 0.8, 0, 0.3)[j], 0))
    result <- 0.7 * baseline[subject] + 1.5 +
      subject_change[subject] + tx[subject] * visit +
      0.12 * visit + 0.15 * sex[subject] + 0.01 * age[subject] +
      rnorm(length(subject), 0, 0.35)
    result[visit == 0L] <- baseline
    result
  }, numeric(nrow(pheno))))
  omics <- data.frame(ANALYTE_NAME = paste0("synthetic_probe_", seq_len(4L)),
                      values, check.names = FALSE)
  names(omics)[-1L] <- pheno$SAMPLE_ID
  list(pheno = pheno, omics = omics,
       treatment_labels = c(`0` = "Usual care", `1` = "Diet", `2` = "Diet + exercise"))
}

run_synthetic <- function(data, treatment_labels = data$treatment_labels, ...) quiet(FAST_omics_WAS(
  data$pheno, data$omics, omics_type = "other", additional_covariates = "AGE",
  treatment_labels = treatment_labels, n_cores = 1L, ...
))

# Assemble reference model data directly from the input, independently of FAST's
# validation, baseline matching, model construction, and contrast helpers.
reference_effects <- function(data, response, mixed) {
  p <- data$pheno
  codes <- sort(unique(as.numeric(as.character(p$TREATMENT_GROUP))))
  baseline <- p[p$FU == 0L, , drop = FALSE]
  followup <- p[p$FU != 0L, , drop = FALSE]
  rows <- lapply(seq_len(nrow(data$omics)), function(j) {
    md <- followup
    md$baseline <- as.numeric(data$omics[j, baseline$SAMPLE_ID[
      match(md$SUBJECT_ID, baseline$SUBJECT_ID)]])
    md$response <- as.numeric(data$omics[j, md$SAMPLE_ID])
    if (response == "change") md$response <- md$response - md$baseline
    md$arm <- factor(md$TREATMENT_GROUP, levels = as.character(codes))
    md$FU <- droplevels(md$FU)
    fit <- if (mixed) quiet(lmerTest::lmer(
      response ~ arm * FU + FEMALE + baseline + AGE + (1 | SUBJECT_ID),
      data = md, REML = FALSE, control = lme4::lmerControl(calc.derivs = FALSE)
    )) else lm(response ~ arm + FEMALE + baseline + AGE, data = md)
    contrasts <- lapply(codes[codes != 0], function(code) {
      as.numeric(codes == code) - as.numeric(codes == 0)
    })
    names(contrasts) <- paste0(codes[codes != 0], "_vs_0")
    means <- if (mixed) emmeans::emmeans(fit, ~ arm | FU) else
      emmeans::emmeans(fit, ~ arm)
    result <- as.data.frame(emmeans::contrast(means, contrasts, adjust = "none"))
    if (!mixed) result$FU <- 1L
    data.frame(ANALYTE_NAME = data$omics$ANALYTE_NAME[j],
               FU = as.integer(as.character(result$FU)),
               CONTRAST = result$contrast, EFFECT_SIZE = result$estimate,
               SE = result$SE, P_VALUE = result$p.value)
  })
  do.call(rbind, rows)
}

compare_reference <- function(actual, expected, label) {
  key <- function(x) paste(x$ANALYTE_NAME, x$FU, x$CONTRAST, sep = "/")
  check(!anyDuplicated(key(actual)) && setequal(key(actual), key(expected)),
        paste(label, "has exactly the requested contrast rows"))
  expected <- expected[match(key(actual), key(expected)), , drop = FALSE]
  check(all(actual$STATUS == "ok") &&
          all(vapply(c("EFFECT_SIZE", "SE", "P_VALUE"), function(column) {
            near(actual[[column]], expected[[column]], tolerance = 1e-6)
          }, logical(1))), paste(label, "matches independent model/contrasts"))
}

check_bh <- function(table, filtered = NULL) {
  column <- if (is.null(filtered)) "BH_P_VALUE" else "BH_P_VALUE_FILTERED"
  expected <- rep(NA_real_, nrow(table))
  groups <- interaction(table$FU, table$CONTRAST, drop = TRUE)
  for (group in levels(groups)) {
    index <- which(groups == group &
                     (is.null(filtered) | table$ANALYTE_NAME %in% filtered))
    expected[index] <- p.adjust(table$P_VALUE[index], method = "BH")
  }
  near(table[[column]], expected)
}

data_multi <- make_synthetic()
original_inputs <- serialize(data_multi, NULL)
data_single <- data_multi
data_single$pheno <- subset(data_single$pheno, FU != 2L)
results <- list()
for (model in c("single", "multi")) {
  data <- if (model == "single") data_single else data_multi
  results[[model]] <- run_synthetic(data)
  for (response in c("change", "level")) {
    analysis <- results[[model]][[paste0("analysis_", response)]]
    expected <- reference_effects(data, response, model == "multi")
    compare_reference(analysis$all$treatment_effects, expected,
                      paste(model, response))
    check(all(vapply(analysis, function(x) check_bh(x$treatment_effects), logical(1))),
          paste(model, response, "BH is per follow-up and contrast"))
    filtered <- data$omics$ANALYTE_NAME[c(1L, 3L)]
    corrected <- .add_filtered_bh_correction(analysis, filtered)
    check(all(vapply(corrected, function(x) check_bh(x$treatment_effects, filtered),
                     logical(1))), paste(model, response, "filtered BH is per contrast"))
  }
}
check(!is.null(results$multi$metadata$schema_version) &&
        all(c("ARM", "TREATMENT", "REFERENCE", "CONTRAST") %in%
              names(results$multi$metadata$treatments)), "versioned treatment metadata exists")
check(identical(results$multi$metadata$treatment_labels, data_multi$treatment_labels) &&
        identical(results$multi$metadata$treatments$CONTRAST, c("1_vs_0", "2_vs_0")),
      "metadata keeps the full label map and stable contrast codes")
table <- results$multi$analysis_change$all$treatment_effects
check(all(table$REFERENCE == "Usual care") &&
        all(table$TREATMENT[table$CONTRAST == "1_vs_0"] == "Diet") &&
        all(table$TREATMENT[table$CONTRAST == "2_vs_0"] == "Diet + exercise"),
      "treatment effects display readable treatment and reference labels")
check(identical(serialize(data_multi, NULL), original_inputs),
      "analysis does not mutate phenotype, omics, or the supplied label map")

# A single observed measurement has undefined variance and no complete
# baseline/follow-up pair; it must not abort valid analytes during validation.
local({
  valid <- data_single
  valid$omics <- valid$omics[1L, , drop = FALSE]
  reference <- run_synthetic(valid)
  singleton <- valid
  sparse_row <- valid$omics
  sparse_row$ANALYTE_NAME <- "singleton_probe"
  sparse_row[, -1L] <- NA_real_
  sparse_row[1L, valid$pheno$SAMPLE_ID[valid$pheno$FU == 1L][1L]] <- 0.25
  singleton$omics <- rbind(valid$omics, sparse_row)
  output <- run_synthetic(singleton)
  for (branch in c("analysis_change", "analysis_level")) {
    for (stratum in c("all", "male", "female")) {
      effects <- output[[branch]][[stratum]]$treatment_effects
      sparse <- effects[effects$ANALYTE_NAME == "singleton_probe", , drop = FALSE]
      available <- effects[effects$ANALYTE_NAME != "singleton_probe", , drop = FALSE]
      row.names(available) <- NULL
      check(nrow(sparse) == 2L && all(sparse$STATUS == "not_estimable") &&
              all(is.na(as.matrix(sparse[c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")]))),
            paste(branch, stratum, "retains unavailable singleton-measurement contrasts"))
      check(isTRUE(all.equal(available, reference[[branch]][[stratum]]$treatment_effects,
                            tolerance = 1e-8)),
            paste(branch, stratum, "singleton measurement leaves valid analyte inference unchanged"))
    }
  }
})

reordered <- data_single
reordered$pheno$TREATMENT_GROUP <- factor(reordered$pheno$TREATMENT_GROUP,
                                        levels = c("2", "0", "1", "99"))
reordered$treatment_labels <- reordered$treatment_labels[c("2", "0", "1")]
reordered_result <- run_synthetic(reordered)
check(isTRUE(all.equal(results$single, reordered_result, tolerance = 1e-8)),
      "factor levels and label-map order cannot change effects or treatment mapping")

local({
  previous <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(previous), add = TRUE)
  sum_coded <- run_synthetic(data_single)
  for (response in c("change", "level")) {
    branch <- paste0("analysis_", response)
    columns <- c("EFFECT_SIZE", "SE", "P_VALUE")
    check(all(vapply(columns, function(column) {
      near(sum_coded[[branch]]$all$treatment_effects[[column]],
           results$single[[branch]]$all$treatment_effects[[column]])
    }, logical(1))), paste(response, "contrasts do not depend on global factor coding"))
  }
})
aliased <- data_single
aliased$omics <- aliased$omics[1L, , drop = FALSE]
aliased$pheno$AGE <- as.integer(aliased$pheno$TREATMENT_GROUP != 0L)
aliased_output <- run_synthetic(aliased)
check(all(vapply(c("analysis_change", "analysis_level"), function(branch) {
  table <- aliased_output[[branch]]$all$treatment_effects
  nrow(table) == 2L && all(table$STATUS == "not_estimable") &&
    all(is.na(table$EFFECT_SIZE))
}, logical(1))), "aliased covariate cannot produce arbitrary adjusted treatment effects")

validated <- quiet(.validate_pheno(data_single$pheno, "AGE",
                                   treatment_labels = data_single$treatment_labels))
check(nrow(validated$all) == nrow(data_single$pheno),
      "label validation retains every valid phenotype row")
factor_subjects <- data_single$pheno
factor_subjects$SUBJECT_ID <- factor(factor_subjects$SUBJECT_ID,
  levels = c(unique(factor_subjects$SUBJECT_ID), "unused_synthetic_subject"))
check(nrow(quiet(.validate_pheno(factor_subjects, "AGE",
                  treatment_labels = data_single$treatment_labels))$all) == nrow(factor_subjects),
      "unused subject factor levels cannot invalidate treatment assignment")
invalid_maps <- list(
  "missing multi-arm map" = NULL,
  "unnamed map" = unname(data_single$treatment_labels),
  "list map" = as.list(data_single$treatment_labels),
  "numeric labels" = c(`0` = 10, `1` = 20, `2` = 30),
  "missing control label" = c(`1` = "Diet", `2` = "Exercise"),
  "missing active label" = c(`0` = "Control", `1` = "Diet"),
  "unobserved extra code" = c(data_single$treatment_labels, `3` = "Another arm"),
  "duplicate code keys" = c(`0` = "Control", `1` = "Diet", `1` = "Exercise"),
  "duplicate labels" = c(`0` = "Same", `1` = "Same", `2` = "Exercise"),
  "blank label" = c(`0` = "Control", `1` = "Diet", `2` = "  "),
  "missing label" = c(`0` = "Control", `1` = "Diet", `2` = NA_character_),
  "blank code key" = setNames(c("Control", "Diet", "Exercise"), c("0", "1", "")),
  "missing code key" = setNames(c("Control", "Diet", "Exercise"), c("0", "1", NA_character_))
)
for (label in names(invalid_maps)) {
  expect_error(quiet(.validate_pheno(data_single$pheno, "AGE",
                treatment_labels = invalid_maps[[label]])), paste(label, "fails validation"))
}
for (invalid_code in c("-1", "0.5", "Inf", "Diet", NA_character_)) {
  bad <- data_single$pheno
  bad$TREATMENT_GROUP <- as.character(bad$TREATMENT_GROUP)
  bad$TREATMENT_GROUP[1L] <- invalid_code
  expect_error(quiet(.validate_pheno(bad, "AGE",
                treatment_labels = data_single$treatment_labels)),
               paste("invalid treatment code", invalid_code, "fails validation"))
}
bad <- data_single$pheno
active <- which(bad$TREATMENT_GROUP == 1L)
bad$TREATMENT_GROUP[active[1L]] <- "2"
expect_error(quiet(.validate_pheno(bad, "AGE", treatment_labels = data_single$treatment_labels)),
             "within-subject active-arm change fails validation")
bad <- data_single$pheno
bad$TREATMENT_GROUP[active[1L]] <- "0"
expect_error(quiet(.validate_pheno(bad, "AGE", treatment_labels = data_single$treatment_labels)),
             "within-subject control-status change fails validation")

binary <- data_single
binary$pheno$TREATMENT_GROUP <- factor(as.integer(binary$pheno$TREATMENT_GROUP != 0L))
binary$treatment_labels <- c(`0` = "Usual care", `1` = "Intervention")
legacy <- quiet(FAST_omics_WAS(binary$pheno, binary$omics, omics_type = "other",
                              additional_covariates = "AGE", n_cores = 1L))
named_binary <- run_synthetic(binary)
check(identical(legacy$metadata$treatment_labels, c(`0` = "Control", `1` = "Treatment")) &&
        all(legacy$analysis_change$all$treatment_effects$REFERENCE == "Control") &&
        all(legacy$analysis_change$all$treatment_effects$TREATMENT == "Treatment"),
      "binary calls without a map receive readable default labels")
for (response in c("change", "level")) {
  branch <- paste0("analysis_", response)
  table <- legacy[[branch]]$all$treatment_effects
  compare_reference(table, reference_effects(binary, response, FALSE),
                    paste("legacy binary", response))
  columns <- c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")
  check(isTRUE(all.equal(table[columns], named_binary[[branch]]$all$treatment_effects[columns],
                        tolerance = 1e-8)), paste("binary numeric compatibility", response))
}

sparse <- data_single
sparse$pheno$TREATMENT_GROUP <- factor(c(`0` = "0", `1` = "2", `2` = "7")[
  as.character(sparse$pheno$TREATMENT_GROUP)])
sparse$treatment_labels <- c(`0` = "Usual care", `2` = "Diet", `7` = "Diet + exercise")
sparse_output <- run_synthetic(sparse)
compare_reference(sparse_output$analysis_change$all$treatment_effects,
                  reference_effects(sparse, "change", FALSE),
                  "nonconsecutive treatment codes")
check(identical(sparse_output$metadata$treatments$ARM, c("2", "7")) &&
        identical(sparse_output$metadata$treatments$CONTRAST, c("2_vs_0", "7_vs_0")),
      "nonconsecutive arm codes are preserved without renumbering")
single_code <- binary
single_code$pheno$TREATMENT_GROUP <- factor(ifelse(single_code$pheno$TREATMENT_GROUP == 0L, 0L, 7L))
single_code$omics <- single_code$omics[1L, , drop = FALSE]
single_code_output <- run_synthetic(single_code, treatment_labels = NULL)
check(identical(single_code_output$metadata$treatment_labels,
                 c(`0` = "Control", `7` = "Treatment")) &&
        all(single_code_output$analysis_change$all$treatment_effects$CONTRAST == "7_vs_0"),
      "one active arm need not have code one or a supplied label map")

# Mapping is validated once against full input, before covariate filtering.
filtered_arm <- data_single
filtered_arm$pheno$AGE[filtered_arm$pheno$TREATMENT_GROUP == 2L] <- NA_real_
filtered_output <- run_synthetic(filtered_arm)
table <- filtered_output$analysis_change$all$treatment_effects
check(identical(filtered_output$metadata$treatment_labels, data_single$treatment_labels) &&
        all(table$STATUS[table$CONTRAST == "2_vs_0"] == "not_estimable") &&
        all(table$TREATMENT[table$CONTRAST == "2_vs_0"] == "Diet + exercise") &&
        all(table$STATUS[table$CONTRAST == "1_vs_0"] == "ok"),
      "an arm lost after covariate filtering keeps its label and contrast rows")

# An arm absent from one stratum or analyte must not remove the other contrast.
missing <- data_multi
is_b <- missing$pheno$TREATMENT_GROUP == 2L
missing$pheno$FEMALE[is_b] <- "1"
missing$omics[2L, missing$pheno$SAMPLE_ID[is_b]] <- NA_real_
missing$omics[3L, missing$pheno$SAMPLE_ID[missing$pheno$TREATMENT_GROUP == 0L]] <- NA_real_
missing$omics[4L, missing$pheno$SAMPLE_ID[is_b & missing$pheno$FU == 2L]] <- NA_real_
for (model in c("single", "multi")) {
  input <- missing
  if (model == "single") input$pheno <- subset(input$pheno, FU != 2L)
  output <- run_synthetic(input)
  for (response in c("change", "level")) {
    branch <- output[[paste0("analysis_", response)]]
    all <- branch$all$treatment_effects
    male <- branch$male$treatment_effects
    n_fu <- if (model == "single") 1L else 2L
    absent <- all$ANALYTE_NAME == "synthetic_probe_2" & all$CONTRAST == "2_vs_0"
    retained <- all$ANALYTE_NAME == "synthetic_probe_2" & all$CONTRAST == "1_vs_0"
    no_control <- all$ANALYTE_NAME == "synthetic_probe_3"
    check(nrow(all) == 8L * n_fu && nrow(male) == 8L * n_fu &&
            all(all$STATUS[absent] == "not_estimable") &&
            all(is.na(all$EFFECT_SIZE[absent])) && all(all$STATUS[retained] == "ok") &&
            all(all$STATUS[no_control] == "not_estimable") &&
            all(male$STATUS[male$CONTRAST == "2_vs_0"] == "not_estimable"),
          paste(model, response, "retains missing-arm/control rows and estimable contrast"))
    check(all(male$TREATMENT[male$CONTRAST == "2_vs_0"] == "Diet + exercise") &&
            all(male$REFERENCE == "Usual care"),
          paste(model, response, "sparse-stratum output keeps full-study arm labels"))
    if (model == "multi") {
      partial <- all[all$ANALYTE_NAME == "synthetic_probe_4" & all$CONTRAST == "2_vs_0", ]
      check(partial$STATUS[partial$FU == 1L] == "ok" &&
              partial$STATUS[partial$FU == 2L] == "not_estimable",
            paste(response, "missing arm at only one visit is contrast-specific"))
    }
  }
}

# The same batch save/resume behavior serves binary and multi-arm studies.
local({
  directory <- tempfile("track111_multiarm_test_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  binary_single <- binary
  binary_single$treatment_labels <- NULL
  binary_multi <- data_multi
  binary_multi$pheno$TREATMENT_GROUP <- factor(as.integer(binary_multi$pheno$TREATMENT_GROUP != 0L))
  binary_multi$treatment_labels <- NULL
  cases <- list(
    binary_lm = list(data = binary_single, reference = legacy),
    binary_lme4 = list(data = binary_multi, reference = run_synthetic(binary_multi)),
    multiarm_lm = list(data = data_single, reference = results$single),
    multiarm_lme4 = list(data = data_multi, reference = results$multi)
  )
  for (name in names(cases)) {
    data <- cases[[name]]$data
    reference <- cases[[name]]$reference
    checkpoint <- file.path(directory, name)
    cached <- run_synthetic(data, checkpoint_dir = checkpoint,
                             checkpoint_batch_size = 2L)
    check(isTRUE(all.equal(reference, cached, tolerance = 1e-8)),
          paste(name, "checkpointed run matches fresh output"))
    batches <- list.files(checkpoint, "^batch_[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE)
    check(length(batches) > 1L &&
            length(list.files(checkpoint, recursive = TRUE)) == length(batches),
          paste(name, "writes only batch files"))
    modification_times <- file.info(batches)$mtime

    resumed <- run_synthetic(data, checkpoint_dir = checkpoint,
                              checkpoint_batch_size = 2L)
    check(isTRUE(all.equal(reference, resumed, tolerance = 1e-8)) &&
            identical(modification_times, file.info(batches)$mtime),
          paste(name, "fully cached resume matches fresh output without rewriting batches"))

    relabeled <- data
    relabeled$treatment_labels <- reference$metadata$treatment_labels
    relabeled$treatment_labels["0"] <- "Comparison group"
    relabeled$treatment_labels[length(relabeled$treatment_labels)] <- "Exercise program"
    relabeled_fresh <- run_synthetic(relabeled)
    relabeled_cached <- run_synthetic(relabeled, checkpoint_dir = checkpoint,
                                       checkpoint_batch_size = 2L)
    check(isTRUE(all.equal(relabeled_fresh, relabeled_cached, tolerance = 1e-8)) &&
            identical(modification_times, file.info(batches)$mtime),
          paste(name, "applies current labels without rewriting cached batches"))

    unlink(batches[1L])
    untouched <- batches[-1L]
    modification_times <- file.info(untouched)$mtime
    partial <- run_synthetic(data, checkpoint_dir = checkpoint,
                              checkpoint_batch_size = 2L)
    check(file.exists(batches[1L]) &&
            isTRUE(all.equal(reference, partial, tolerance = 1e-8)) &&
            identical(modification_times, file.info(untouched)$mtime),
          paste(name, "partial resume recreates the missing batch and preserves completed batches"))
  }

  # Public input validation still runs before cached results are read.
  expect_error(quiet(FAST_omics_WAS(data_single$pheno, data_single$omics,
      omics_type = "other", additional_covariates = "AGE", n_cores = 1L,
      checkpoint_dir = checkpoint, checkpoint_batch_size = 2L)),
    "missing multi-arm labels cannot reach cached results", "label")

  # Resume the final (multi-arm, repeated-follow-up) case in a fresh process.
  unlink(batches[1L])
  untouched <- batches[-1L]
  modification_times <- file.info(untouched)$mtime
  input_path <- file.path(directory, "input.rds")
  output_path <- file.path(directory, "output.rds")
  script <- file.path(directory, "resume.R")
  saveRDS(list(data = data, checkpoint = checkpoint, output = output_path), input_path)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "source(args[1L])",
    "input <- readRDS(args[2L])",
    "result <- suppressWarnings(suppressMessages(FAST_omics_WAS(",
    "  input$data$pheno, input$data$omics, omics_type = 'other',",
    "  additional_covariates = 'AGE', treatment_labels = input$data$treatment_labels,",
    "  n_cores = 2L, checkpoint_dir = input$checkpoint, checkpoint_batch_size = 2L)))",
    "saveRDS(result, input$output)"
  ), script)
  log <- file.path(directory, "resume.log")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(script),
                      shQuote(file.path(.track111_root, "main.R")), shQuote(input_path)),
                    stdout = log, stderr = log)
  if (status != 0L) cat(readLines(log), sep = "\n")
  check(status == 0L && file.exists(output_path),
        "fresh R session resumes with a changed core count")
  check(isTRUE(all.equal(cached, readRDS(output_path), tolerance = 1e-8)) &&
          identical(modification_times, file.info(untouched)$mtime),
        "partial resume preserves cached batches and reproduces fresh output")
})

source("CALERIE/DNAm/EWAS_DNAm/dunedinpace_covariance.R")
table <- results$single$analysis_change$all$treatment_effects
expect_error(.dpace_extract(list(treatment_effects = table),
                            unique(table$ANALYTE_NAME), 1L),
             "binary covariance helper rejects multiple contrasts", "single treatment contrast")
single <- table[table$CONTRAST == "1_vs_0", ]
check(nrow(.dpace_extract(list(treatment_effects = single),
                          unique(single$ANALYTE_NAME), 1L)$effects) == nrow(single),
      "binary covariance helper accepts any one contrast")
expect_error(.dpace_extract(list(treatment_effects = rbind(single, single[1L, ])),
                            unique(single$ANALYTE_NAME), 1L),
             "binary covariance helper rejects duplicate keys", "unique analyte/FU")
cat("All synthetic multi-arm tests passed.\n")
