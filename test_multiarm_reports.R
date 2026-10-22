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
