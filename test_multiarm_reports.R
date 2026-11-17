# Synthetic checks only: no participant files or study inputs are read.
source("reporting_helpers.R")
source("plotting_helpers.R")
source("proteomics_translation_helpers.R")
.log_111 <- function(...) invisible(NULL)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(message), grepl(pattern, message))
}

arms <- rep(c("0", "1", "2"), c(6L, 8L, 4L))
subject <- rep(seq_along(arms), each = 2L)
pheno <- data.frame(
  SAMPLE_ID = paste0("synthetic-", seq_along(subject)),
  SUBJECT_ID = paste0("person-", subject), FU = factor(rep(0:1, length(arms))),
  FEMALE = factor(subject %% 2L), TREATMENT_GROUP = factor(arms[subject]),
  age = 30 + subject, smoker = subject %% 3L == 0L,
  stringsAsFactors = FALSE
)
treatments <- data.frame(
  ARM = c("1", "2"), TREATMENT = c("Arm A", "Arm B"), REFERENCE = "Usual care",
  CONTRAST = c("1_vs_0", "2_vs_0"), stringsAsFactors = FALSE
)
values <- rbind(
  seq_along(subject) / 10 + c(0, 2, -3)[as.integer(pheno$TREATMENT_GROUP)],
  sin(seq_along(subject)) + c(0, 4, -2)[as.integer(pheno$TREATMENT_GROUP)]
)
# The first analyte is unavailable in arm B; later analytes must still be reported.
values[1, pheno$TREATMENT_GROUP == "2"] <- NA_real_
omics <- data.frame(ANALYTE_NAME = c("protein-1", "protein-2"), values,
                    check.names = FALSE)
names(omics)[-1L] <- pheno$SAMPLE_ID
# Verify sample alignment independently of input column order.
omics <- omics[, c("ANALYTE_NAME", rev(pheno$SAMPLE_ID)), drop = FALSE]
pheno_list <- list(all = pheno, male = pheno[pheno$FEMALE == "0", ],
                   female = pheno[pheno$FEMALE == "1", ], treatments = treatments)
omics_list <- list(all = omics, male = omics, female = omics)
reports <- suppressWarnings(.generate_reports(pheno_list, omics_list, c("age", "smoker")))

stopifnot(identical(reports$metadata$treatments, treatments),
          identical(reports$metadata$treatment_labels,
                    c(`0` = "Usual care", `1` = "Arm A", `2` = "Arm B")),
          all(reports$arm_counts$TREATMENT[reports$arm_counts$IS_CONTROL] == "Usual care"))
stopifnot(all(reports$pheno_summary$N_CONTROL == 3L),
          all(reports$pheno_summary$N_TREATMENT == 6L),
          sum(reports$arm_counts$N_SUBJECTS) == nrow(pheno),
          all(c("omics_FU0_Tx0", "omics_FU0_Tx1", "omics_FU0_Tx2") %in%
                names(reports$variable_summaries$all)))
randomization <- reports$randomization_reports$omics_report
stopifnot(setequal(randomization$CONTRAST, treatments$CONTRAST),
          nrow(randomization[randomization$CONTRAST == "2_vs_0", ]) == 1L,
          all(randomization$REFERENCE == "Usual care"))
for (i in seq_len(nrow(treatments))) {
  arm <- treatments$ARM[i]
  selected <- randomization[randomization$CONTRAST == treatments$CONTRAST[i], ]
  for (j in seq_len(nrow(selected))) {
    analyte <- match(selected$ANALYTE_NAME[j], omics$ANALYTE_NAME)
    control_samples <- pheno$SAMPLE_ID[pheno$FU == "0" & pheno$TREATMENT_GROUP == "0"]
    treatment_samples <- pheno$SAMPLE_ID[pheno$FU == "0" & pheno$TREATMENT_GROUP == arm]
    expected <- mean(as.numeric(omics[analyte, treatment_samples]), na.rm = TRUE) -
      mean(as.numeric(omics[analyte, control_samples]), na.rm = TRUE)
    stopifnot(isTRUE(all.equal(selected$MEAN_DIFFERENCE[j], expected)))
  }
  stopifnot(isTRUE(all.equal(selected$BH_P_VALUE, p.adjust(selected$P_VALUE, "BH"))))
}
covariates <- reports$randomization_reports$covariate_report
stopifnot(all(table(covariates$CONTRAST) == 3L),
          setequal(covariates$TREATMENT, treatments$TREATMENT))
binary_pheno <- pheno[pheno$TREATMENT_GROUP != "2", , drop = FALSE]
binary_pheno$TREATMENT_GROUP <- droplevels(binary_pheno$TREATMENT_GROUP)
binary_report <- .create_randomization_report(binary_pheno, omics)
binary_reference <- .create_pairwise_randomization_report(binary_pheno, omics)
stopifnot(isTRUE(all.equal(binary_report[, names(binary_reference)], binary_reference)),
          all(binary_report$CONTRAST == "1_vs_0"),
          all(binary_report$TREATMENT == "Treatment"),
          all(binary_report$REFERENCE == "Control"))
binary_list <- list(all = binary_pheno, male = NULL, female = NULL)
binary_reports <- suppressWarnings(.generate_reports(binary_list, omics_list, c("age", "smoker")))
stopifnot(identical(binary_reports$metadata$treatment_labels,
                    c(`0` = "Control", `1` = "Treatment")),
          all(binary_reports$arm_counts$TREATMENT[binary_reports$arm_counts$IS_CONTROL] == "Control"))
absent_list <- list(all = binary_pheno, male = NULL, female = NULL, treatments = treatments)
absent_reports <- suppressWarnings(.generate_reports(absent_list, omics_list, c("age", "smoker")))
stopifnot(all(absent_reports$arm_counts$N_SUBJECTS[absent_reports$arm_counts$ARM == "2"] == 0L),
          !"omics_FU0_Tx2" %in% names(absent_reports$variable_summaries$all),
          !"2_vs_0" %in% absent_reports$randomization_reports$omics_report$CONTRAST)
absent_covariates <- subset(absent_reports$randomization_reports$covariate_report,
                            CONTRAST == "2_vs_0")
stopifnot(nrow(absent_covariates) > 0L,
          all(is.na(absent_covariates$P_VALUE)), all(is.na(absent_covariates$STATISTIC)),
          length(absent_covariates$SUMMARY_CONTROL) == nrow(absent_covariates))
empty_report <- .report_pairwise(binary_pheno, treatments,
                                 function(...) data.frame(ANALYTE_NAME = character()))
stopifnot(nrow(empty_report) == 0L,
          all(c("TREATMENT", "REFERENCE", "CONTRAST") %in% names(empty_report)))

local({
  # Frozen pre-allocation implementations. Keep these independent of the live
  # helpers so optimizations must preserve historical values and R attributes.
  original_summary <- function(sample_ids, omics_df, verbose = FALSE,
                                        progress_label = "reports variable summary",
                                        progress_every = 100000L) {

    analyte_names <- omics_df$ANALYTE_NAME

    omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
    omics_numeric <- omics_numeric[, colnames(omics_numeric) %in% sample_ids, drop = FALSE]

    # Initialize results data.frame
    report <- data.frame(
      ANALYTE_NAME = analyte_names,
      N_NONMISSING = NA_integer_,
      MEAN = NA_real_,
      MEDIAN = NA_real_,
      SD = NA_real_,
      MIN = NA_real_,
      MAX = NA_real_,
      stringsAsFactors = FALSE
    )

    .log_111(
      verbose,
      paste0(progress_label, ": 0/", length(analyte_names), " analytes (serial)")
    )

    # Calculate per-analyte statistics
    for (i in seq_along(analyte_names)) {
      analyte_values <- as.numeric(omics_numeric[i, ])

      report$N_NONMISSING[i] <- sum(!is.na(analyte_values))
      report$MEAN[i] <- mean(analyte_values, na.rm = TRUE)
      report$MEDIAN[i] <- median(analyte_values, na.rm = TRUE)
      report$SD[i] <- sd(analyte_values, na.rm = TRUE)
      report$MIN[i] <- min(analyte_values, na.rm = TRUE)
      report$MAX[i] <- max(analyte_values, na.rm = TRUE)

      if (i %% progress_every == 0L || i == length(analyte_names)) {
        .log_111(
          verbose,
          paste0(progress_label, ": ", i, "/", length(analyte_names), " analytes")
        )
      }
    }

    return(report)
  }

  original_randomization <- function(pheno_df, omics_df, verbose = FALSE,
                                                    progress_every = 100000L) {

    # Filter to baseline (FU=0)
    pheno_baseline <- pheno_df[pheno_df$FU == 0, ]
    pheno_baseline <- pheno_baseline[pheno_baseline$SAMPLE_ID %in% colnames(omics_df), , drop = FALSE]

    # Filter omics to baseline samples
    baseline_sample_ids <- pheno_baseline$SAMPLE_ID
    omics_baseline <- omics_df[, baseline_sample_ids, drop = FALSE]

    # Initialize results data.frame
    results <- data.frame(
      ANALYTE_NAME = character(),
      MEAN_DIFFERENCE = numeric(),
      COHENS_D = numeric(),
      SE = numeric(),
      P_VALUE = numeric(),
      stringsAsFactors = FALSE
    )

    # T-test per analyte
    analyte_names <- omics_df$ANALYTE_NAME
    .log_111(
      verbose,
      paste0("reports randomization: 0/", length(analyte_names), " analytes (serial)")
    )

    for (i in seq_along(analyte_names)) {
      tryCatch({
        # Extract values for this analyte
        analyte_values <- as.numeric(omics_baseline[i, ])
        control_status <- pheno_baseline$TREATMENT_GROUP

        # Split by group
        values_control <- analyte_values[control_status == 0]
        values_treatment <- analyte_values[control_status == 1]

        # Remove NAs
        values_control <- na.omit(values_control)
        values_treatment <- na.omit(values_treatment)

        # Skip if insufficient samples in either group
        if (length(values_control) < 2 || length(values_treatment) < 2) {
          next
        }

        # Calculate statistics
        n_control <- length(values_control)
        n_treatment <- length(values_treatment)
        mean_control <- mean(values_control)
        mean_treatment <- mean(values_treatment)
        var_control <- var(values_control)
        var_treatment <- var(values_treatment)

        # Mean difference (Treatment - Control)
        mean_diff <- mean_treatment - mean_control

        # Standard Error (Welch's approach - doesn't assume equal variances)
        se <- sqrt((var_control / n_control) + (var_treatment / n_treatment))

        # Cohen's d (using pooled SD)
        pooled_sd <- sqrt(((n_control - 1) * var_control + (n_treatment - 1) * var_treatment) /
                            (n_control + n_treatment - 2))
        cohens_d <- mean_diff / pooled_sd

        # Welch's t-test
        tt <- t.test(values_treatment, values_control, var.equal = FALSE)

        # Add to results
        results <- rbind(results, data.frame(
          ANALYTE_NAME = analyte_names[i],
          MEAN_DIFFERENCE = mean_diff,
          COHENS_D = cohens_d,
          SE = se,
          P_VALUE = tt$p.value,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        warning("Error processing analyte '", analyte_names[i], "': ", e$message)
      })

      if (i %% progress_every == 0L || i == length(analyte_names)) {
        .log_111(
          verbose,
          paste0("reports randomization: ", i, "/", length(analyte_names), " analytes")
        )
      }
    }

    # Apply multiple testing correction (global across all analytes)
    if (nrow(results) > 0) {
      results$BH_P_VALUE <- p.adjust(results$P_VALUE, method = "BH")
      return(results)
    } else {
      return(NULL)
    }
  }

  # Capture interleaved warnings and progress, as well as errors, without changing
  # the helpers' control flow. identical() checks types, attributes and row order.
  capture_report <- function(expr) {
    events <- list()
    error <- NULL
    value <- tryCatch(withCallingHandlers(force(expr),
      warning = function(w) {
        events[[length(events) + 1L]] <<- list("warning", conditionMessage(w), class(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        events[[length(events) + 1L]] <<- list("message", conditionMessage(m), class(m))
        invokeRestart("muffleMessage")
      }), error = function(e) {
        error <<- list(conditionMessage(e), class(e))
        NULL
      })
    list(value = value, events = events, error = error)
  }
  expect_same_report <- function(current, original) {
    stopifnot(identical(current, original))
    current
  }
  event_text <- function(captured, type) {
    vapply(Filter(function(x) identical(x[[1L]], type), captured$events),
           function(x) x[[2L]], character(1L))
  }
  in_environment <- function(fun, env) {
    environment(fun) <- env
    fun
  }
  original <- new.env(parent = environment())
  current <- new.env(parent = environment())
  for (env in list(original, current)) {
    env$.log_111 <- function(verbose, text) {
      if (isTRUE(verbose)) message("[1.1.1] ", text)
      invisible(NULL)
    }
    env$.generate_reports <- in_environment(.generate_reports, env)
    env$.create_randomization_report <- in_environment(.create_randomization_report, env)
  }
  original$.create_omics_data_report <- in_environment(original_summary, original)
  original$.create_pairwise_randomization_report <- in_environment(original_randomization, original)
  current$.create_omics_data_report <- in_environment(.create_omics_data_report, current)
  current$.create_pairwise_randomization_report <- in_environment(.create_pairwise_randomization_report, current)

  # structure() preserves names on the ID column, unlike data.frame(). This is a
  # valid input whose historical row-name behavior also belongs to the contract.
  make_omics <- function(ids, columns) {
    stopifnot(length(ids) == nrow(columns))
    structure(c(list(ANALYTE_NAME = ids), as.list(columns)),
              class = "data.frame", row.names = seq_along(ids))
  }
  compare_summary <- function(omics, samples = names(omics)[-1L], progress_every = 1000L) {
    expect_same_report(
      capture_report(current$.create_omics_data_report(samples, omics, verbose = TRUE,
                     progress_label = "parity summary", progress_every = progress_every)),
      capture_report(original$.create_omics_data_report(samples, omics, verbose = TRUE,
                     progress_label = "parity summary", progress_every = progress_every)))
  }

  # Cross both 2,000-row chunk boundaries, including missing and nonfinite rows.
  n <- 4003L
  summary_columns <- data.frame(s1 = seq_len(n), s2 = sin(seq_len(n)),
                                s3 = cos(seq_len(n)), s4 = seq_len(n) %% 3L == 0L)
  summary_columns[c(1999L, 2000L, 4000L), ] <- NA_real_
  summary_columns[2001L, ] <- c(Inf, -Inf, NaN, NA_real_)
  summary_columns[4001L, ] <- c(NaN, 2, NA_real_, 4)
  summary_omics <- make_omics(rep(c("same", "other"), length.out = n), summary_columns)
  summary_result <- compare_summary(summary_omics)
  stopifnot(is.null(summary_result$error), nrow(summary_result$value) == n,
            identical(summary_result$value$N_NONMISSING[c(1999L, 2000L, 4000L, 4001L)],
                      c(0L, 0L, 0L, 2L)),
            length(event_text(summary_result, "warning")) == 6L,
            identical(event_text(summary_result, "message"), paste0(
              "[1.1.1] parity summary: ", c("0/4003 analytes (serial)",
              "1000/4003 analytes", "2000/4003 analytes", "3000/4003 analytes",
              "4000/4003 analytes", "4003/4003 analytes"), "\n")))
  small_columns <- data.frame(s1 = c(1L, 4L, NA_integer_), s2 = c(2L, 6L, 8L),
                              s3 = c(1.5, NaN, Inf), s4 = c(TRUE, FALSE, TRUE))
  identifiers <- list(
    duplicate = c("same", "same", "third"),
    factor = factor(c("b", "a", "b"), levels = c("unused", "b", "a")),
    integer = c(10L, 10L, 30L), numeric = c(1.5, 1.5, 3.5),
    named = setNames(c("same", "same", "third"), c("x", "x", "x1")))
  for (ids in identifiers) {
    small_omics <- make_omics(ids, small_columns)
    compare_summary(small_omics)
    one_sample <- compare_summary(small_omics, "s1", progress_every = 1L)
    stopifnot(is.null(one_sample$error), all(is.na(one_sample$value$SD)))
    no_samples <- compare_summary(small_omics, "absent", progress_every = 1L)
    stopifnot(is.null(no_samples$error), all(no_samples$value$N_NONMISSING == 0L),
              all(is.nan(no_samples$value$MEAN)), all(is.na(no_samples$value$SD)),
              identical(no_samples$value$MIN, rep(Inf, 3L)),
              identical(no_samples$value$MAX, rep(-Inf, 3L)),
              length(event_text(no_samples, "warning")) == 6L)
  }
  # The historical zero-row constructor errors before it emits any progress.
  empty <- compare_summary(make_omics(character(), small_columns[FALSE, ]))
  stopifnot(!is.null(empty$error), length(empty$events) == 0L)

  # Baseline phenotype order deliberately differs from the omics column order;
  # extra follow-up and absent samples must never enter the baseline comparison.
  sample_names <- paste0("baseline-", seq_len(8L))
  random_pheno <- data.frame(
    SAMPLE_ID = c(sample_names, "follow-up", "absent"),
    FU = factor(c(rep(0L, 8L), 1L, 0L)),
    TREATMENT_GROUP = factor(c(rep(0L, 4L), rep(1L, 4L), 1L, 0L)))
  random_pheno <- random_pheno[c(7L, 2L, 5L, 4L, 1L, 8L, 3L, 6L, 9L, 10L), ]
  compare_randomization <- function(omics, phenotype = random_pheno, progress_every = 1000L) {
    expect_same_report(
      capture_report(current$.create_pairwise_randomization_report(
        phenotype, omics, verbose = TRUE, progress_every = progress_every)),
      capture_report(original$.create_pairwise_randomization_report(
        phenotype, omics, verbose = TRUE, progress_every = progress_every)))
  }
  n <- 2003L
  random_values <- outer(seq_len(n) / 100, c(0.1, 0.2, 0.3, 0.4, 0.2, 0.4, 0.6, 0.8)) +
    matrix(rep(c(1, 2, 4, 5, 1, 3, 5, 9), each = n), nrow = n)
  random_values[2L, 1:3] <- NA_real_                 # insufficient control samples
  random_values[3L, ] <- 1                         # t.test: essentially constant
  random_values[4L, 1L] <- Inf                     # t.test: nonfinite input
  random_values[5L, c(1L, 5L)] <- c(NA_real_, NaN)  # still estimable
  random_values[1000L, ] <- 1                      # error still emits progress
  random_values[n, 5:7] <- NA_real_                # final skip emits no progress
  random_columns <- as.data.frame(random_values)
  names(random_columns) <- sample_names
  random_columns$`follow-up` <- -999
  random_columns$unused <- 999
  random_columns <- random_columns[, rev(names(random_columns)), drop = FALSE]
  random_omics <- make_omics(paste0("analyte-", seq_len(n)), random_columns)
  randomized <- compare_randomization(random_omics)
  stopifnot(is.null(randomized$error), nrow(randomized$value) == n - 5L,
            length(event_text(randomized, "warning")) == 3L,
            identical(event_text(randomized, "message"), paste0(
              "[1.1.1] reports randomization: ",
              c("0/2003 analytes (serial)", "1000/2003 analytes", "2000/2003 analytes"), "\n")),
            identical(randomized$value$BH_P_VALUE,
                      p.adjust(randomized$value$P_VALUE, method = "BH")))
  # Independently check the direction, Welch SE/test and pooled-SD effect size.
  for (row in c(1L, 5L, 2000L, 2002L)) {
    control <- na.omit(random_values[row, 1:4])
    treated <- na.omit(random_values[row, 5:8])
    nc <- length(control)
    nt <- length(treated)
    difference <- mean(treated) - mean(control)
    pooled <- sqrt(((nc - 1L) * var(control) + (nt - 1L) * var(treated)) / (nc + nt - 2L))
    result <- randomized$value[randomized$value$ANALYTE_NAME == paste0("analyte-", row), ]
    # Phenotype order changes summation order, so the independent calculation is
    # compared numerically; the frozen implementation comparison above is exact.
    stopifnot(isTRUE(all.equal(result$MEAN_DIFFERENCE, difference)),
              isTRUE(all.equal(result$SE, sqrt(var(control) / nc + var(treated) / nt))),
              isTRUE(all.equal(result$COHENS_D, difference / pooled)),
              isTRUE(all.equal(result$P_VALUE, t.test(treated, control, var.equal = FALSE)$p.value)))
  }
  for (ids in identifiers) {
    successful <- make_omics(ids, random_columns[c(1L, 6L, 7L), ])
    checked <- compare_randomization(successful, progress_every = 1L)
    stopifnot(is.null(checked$error), nrow(checked$value) == 3L,
              identical(class(checked$value$ANALYTE_NAME), class(ids)))
    if (!is.null(names(ids))) {
      stopifnot(identical(rownames(checked$value), c("x", "x1", "x11")))
    }
  }
  # Names attached to the ID vector can also cause per-row construction errors.
  for (id_names in list(c("", "x", "y"), c(NA_character_, "x", "y"))) {
    compare_randomization(make_omics(setNames(c("a", "b", "c"), id_names),
                                    random_columns[c(1L, 6L, 7L), ]), progress_every = 1L)
  }
  all_skipped <- compare_randomization(make_omics(c("low-n", "also-low-n"),
                                                  random_columns[c(2L, n), ]), progress_every = 1L)
  stopifnot(is.null(all_skipped$value), is.null(all_skipped$error),
            length(event_text(all_skipped, "warning")) == 0L,
            identical(event_text(all_skipped, "message"),
                      "[1.1.1] reports randomization: 0/2 analytes (serial)\n"))
  all_failed <- compare_randomization(make_omics(c("constant", "nonfinite"),
                                                 random_columns[c(3L, 4L), ]), progress_every = 1L)
  stopifnot(is.null(all_failed$value), is.null(all_failed$error),
            length(event_text(all_failed, "warning")) == 2L,
            length(event_text(all_failed, "message")) == 3L)
  no_analytes <- compare_randomization(make_omics(character(), random_columns[FALSE, ]))
  stopifnot(is.null(no_analytes$value), is.null(no_analytes$error))

  # Exercise the existing public report assembly with both arm configurations,
  # both visits, all sex strata, covariates, and missing per-arm analyte values.
  binary_strata <- list(all = binary_pheno,
                        male = binary_pheno[binary_pheno$FEMALE == "0", ],
                        female = binary_pheno[binary_pheno$FEMALE == "1", ])
  for (phenotypes in list(binary_strata, pheno_list)) {
    assembled <- expect_same_report(
      capture_report(current$.generate_reports(phenotypes, omics_list, c("age", "smoker"), verbose = TRUE)),
      capture_report(original$.generate_reports(phenotypes, omics_list, c("age", "smoker"), verbose = TRUE)))
    stopifnot(is.null(assembled$error),
              all(vapply(assembled$value$variable_summaries,
                         function(x) all(c("omics_FU0_Tx0", "omics_FU1_Tx1") %in% names(x)), logical(1L))))
  }
  cat("PASS: exact reporting allocation parity, conditions, progress and public assembly\n")
})

effects <- data.frame(
  ANALYTE_NAME = rep(c("protein-1", "protein-2"), 2L), FU = 1L,
  TREATMENT = rep(treatments$TREATMENT, each = 2L), REFERENCE = "Usual care",
  CONTRAST = rep(treatments$CONTRAST, each = 2L),
  EFFECT_SIZE = c(1, 2, -2, -1), SE = .2, P_VALUE = c(.1, .01, .2, .02),
  BH_P_VALUE = c(.1, .02, .2, .04), STATUS = "ok", stringsAsFactors = FALSE
)
strata <- list(all = list(treatment_effects = effects), male = NULL, female = NULL)
expect_error(plot_volcano(strata), "Multiple treatments")
expect_error(plot_qq(strata), "Multiple treatments")
expect_error(plot_volcano(strata, treatment = "missing"), "Unknown treatment")
plot <- plot_volcano(strata, treatment = "Arm B")
stopifnot(grepl("Arm B vs Usual care", plot$labels$title),
          all(na.omit(plot$data$CONTRAST) == "2_vs_0"))
legacy <- strata
legacy$all$treatment_effects <- effects[1:2, setdiff(names(effects),
                                                   c("TREATMENT", "REFERENCE", "CONTRAST"))]
legacy_plot <- plot_volcano(legacy, treatment = "Treatment")
stopifnot(inherits(legacy_plot, "ggplot"),
          grepl("Treatment vs Control", legacy_plot$labels$title))
missing_contrast <- strata
missing_contrast$all$treatment_effects$CONTRAST <- NULL
expect_error(plot_volcano(missing_contrast), "Missing CONTRAST")
missing_contrast$all$treatment_effects$TREATMENT <- NULL
missing_contrast$all$treatment_effects$REFERENCE <- NULL
expect_error(plot_volcano(missing_contrast), "Missing CONTRAST")
conflicting_label <- strata
conflicting_label$all$treatment_effects$TREATMENT[1] <- "other arm"
expect_error(plot_volcano(conflicting_label, treatment = "Arm A"), "exactly one TREATMENT")

local({
  # Explicit membership wins over both missing and misleading filtered inference.
  members <- effects
  members$IN_FILTERED_SET <- rep(c(TRUE, FALSE), 2L)
  members$BH_P_VALUE_FILTERED <- c(NA_real_, .02, NA_real_, .04)
  member_results <- list(all = list(treatment_effects = members), male = NULL, female = NULL)
  stopifnot(identical(.treatment_filtered_membership(members), members$IN_FILTERED_SET))
  member_plot <- plot_volcano(member_results, treatment = "Arm A", probe_set = "filtered")
  stopifnot(identical(unique(na.omit(member_plot$data$ANALYTE_NAME)), "protein-1"))

  legacy_members <- members
  legacy_members$IN_FILTERED_SET <- NULL
  legacy_members$STATUS <- NULL
  stopifnot(identical(.treatment_filtered_membership(legacy_members), rep(c(FALSE, TRUE), 2L)))
  legacy_results <- list(all = list(treatment_effects = legacy_members), male = NULL, female = NULL)
  legacy_member_plot <- plot_volcano(legacy_results, treatment = "Arm A", probe_set = "filtered")
  stopifnot(identical(unique(na.omit(legacy_member_plot$data$ANALYTE_NAME)), "protein-2"))
  legacy_members$BH_P_VALUE_FILTERED <- c(Inf, -0.1, 1.1, NA_real_)
  stopifnot(!any(.treatment_filtered_membership(legacy_members)))

  unavailable <- members
  included <- unavailable$IN_FILTERED_SET
  unavailable[included, c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")] <- NA_real_
  unavailable$STATUS[included] <- c("fit_error", "not_estimable")
  unavailable_results <- list(all = list(treatment_effects = unavailable), male = NULL, female = NULL)
  unavailable_plot <- plot_volcano(unavailable_results, treatment = "Arm A", probe_set = "filtered")
  stopifnot(identical(ggplot_build(unavailable_plot)$data[[1L]]$label, "No estimable results"),
            grepl("Arm A vs Usual care", unavailable_plot$labels$title))
  excluded_results <- unavailable_results
  excluded_results$all$treatment_effects$IN_FILTERED_SET <- FALSE
  excluded_plot <- plot_volcano(excluded_results, treatment = "Arm A", probe_set = "filtered")
  stopifnot(identical(ggplot_build(excluded_plot)$data[[1L]]$label, "No probes in filtered set"),
            grepl("Arm A vs Usual care", excluded_plot$labels$title))
  # An unavailable status must not be plotted even if a stale row contains numbers.
  stale <- members
  stale$STATUS[stale$IN_FILTERED_SET] <- "not_estimable"
  stopifnot(!any(.treatment_plot_estimable(stale)[stale$IN_FILTERED_SET]))

  directory <- tempfile("multiarm-unavailable-plots-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  original_qq <- plot_qq
  original_volcano <- plot_volcano
  on.exit(assign("plot_qq", original_qq, envir = .GlobalEnv), add = TRUE)
  on.exit(assign("plot_volcano", original_volcano, envir = .GlobalEnv), add = TRUE)
  calls <- character()
  assign("plot_qq", function(results, stratum = "all", fu = 1, probe_set = "full",
                              treatment = NULL, ...) {
    calls <<- c(calls, paste("qq", stratum, fu, treatment, probe_set, sep = "/"))
    original_qq(results, stratum = stratum, fu = fu, probe_set = probe_set,
                treatment = treatment, ...)
  }, envir = .GlobalEnv)
  assign("plot_volcano", function(results, stratum = "all", fu = 1, probe_set = "full",
                                   treatment = NULL, ...) {
    calls <<- c(calls, paste("volcano", stratum, fu, treatment, probe_set, sep = "/"))
    original_volcano(results, stratum = stratum, fu = fu, probe_set = probe_set,
                     treatment = treatment, ...)
  }, envir = .GlobalEnv)
  second_visit <- unavailable
  second_visit$FU <- 2L
  # No filtered BH value is available anywhere, but all contrast/FU pages remain.
  unavailable$BH_P_VALUE_FILTERED <- NA_real_
  second_visit$BH_P_VALUE_FILTERED <- NA_real_
  pages <- list(all = list(treatment_effects = rbind(unavailable, second_visit)),
                male = list(treatment_effects = second_visit),
                female = list(treatment_effects = NULL))
  device_before <- dev.cur()
  generate_all_plots(pages, figures_dir = directory)
  cells <- c("all/1/Arm A", "all/1/Arm B", "all/2/Arm A", "all/2/Arm B",
             "male/2/Arm A", "male/2/Arm B")
  expected_calls <- unlist(lapply(c("full", "filtered"), function(probe_set) {
    unlist(lapply(c("qq", "volcano"), function(kind) paste(kind, cells, probe_set, sep = "/")))
  }), use.names = FALSE)
  pdfs <- file.path(directory, c("qq_full.pdf", "volcano_full.pdf", "qq_filtered.pdf", "volcano_filtered.pdf"))
  stopifnot(identical(calls, expected_calls), all(file.info(pdfs)$size > 0L),
            identical(dev.cur(), device_before))

  # QQ also emits a labeled page for an excluded set; it no longer stops.
  pdf(file.path(directory, "excluded.pdf"))
  tryCatch(original_qq(excluded_results, treatment = "Arm A", probe_set = "filtered"),
           finally = dev.off())
  stopifnot(identical(dev.cur(), device_before))
  assign("plot_qq", function(...) stop("synthetic plotting failure"), envir = .GlobalEnv)
  expect_error(generate_all_plots(pages, figures_dir = file.path(directory, "failure")),
               "synthetic plotting failure")
  stopifnot(identical(dev.cur(), device_before))
})

temporary <- tempfile("multiarm-report-tests-")
dir.create(temporary)
tryCatch({
  generate_all_plots(strata, figures_dir = temporary)
  stopifnot(all(file.info(file.path(temporary, c("qq_full.pdf", "volcano_full.pdf")))$size > 0))
  result_path <- file.path(temporary, "synthetic-results.rds")
  translation_path <- file.path(temporary, "synthetic-translation.csv")
  saveRDS(list(analysis_change = strata), result_path)
  translation <- data.frame(AptName = c("protein-1", "protein-2"))
  for (column in c("SomaId", "TargetFullName", "Target", "UniProt", "EntrezGeneID",
                   "EntrezGeneSymbol", "Uniprot_Unique", "Symbol_Unique")) {
    translation[[column]] <- c("synthetic-a", "synthetic-b")
  }
  write.csv(translation, translation_path, row.names = FALSE)
  annotated <- get_annotated_treatment_effects(results_rds = result_path,
                                               translation_csv = translation_path)
  stopifnot(nrow(annotated) == 4L,
            setequal(annotated$CONTRAST, treatments$CONTRAST),
            all(annotated$REFERENCE == "Usual care"),
            all(c("TREATMENT", "REFERENCE", "CONTRAST", "STATUS") %in% names(annotated)))
  selected <- get_annotated_treatment_effects(results_rds = result_path,
                                              translation_csv = translation_path,
                                              treatment = "Arm B")
  stopifnot(nrow(selected) == 2L, all(selected$CONTRAST == "2_vs_0"))
  saveRDS(list(analysis_change = legacy), result_path)
  legacy_annotated <- get_annotated_treatment_effects(results_rds = result_path,
                                                     translation_csv = translation_path,
                                                     treatment = "Treatment")
  stopifnot(nrow(legacy_annotated) == 2L)
}, finally = unlink(temporary, recursive = TRUE))
cat("PASS: synthetic multi-arm counts, reports, plots, and annotations\n")
