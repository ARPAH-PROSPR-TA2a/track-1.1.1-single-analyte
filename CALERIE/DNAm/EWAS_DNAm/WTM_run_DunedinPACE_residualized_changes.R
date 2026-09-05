Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

# Fast validation companion to the exact delete-one-subject jackknife.
#
# For the pooled CALERIE stratum and the available DunedinPACE probes, this
# script constructs two subject-aligned factors:
#
# 1. nuisance-residualized beta changes, whose cross-product estimates probe
#    change covariance; and
# 2. treatment-effect influence contributions, whose cross-product estimates
#    covariance among adjusted treatment-effect point estimates.
#
# The second factor uses separate adjusted linear models at each follow-up. It
# is a scalable approximation to the pooled mixed-model jackknife, not an exact
# replacement. Validate it against dunedinpace_raw_beta_jackknife_v1 before
# extending it genome-wide.

.residualized_component <- function(change, baseline, treatment,
                                    nuisance_design) {
  values <- list(
    change = as.numeric(change),
    baseline = as.numeric(baseline),
    treatment = as.numeric(treatment)
  )
  n <- length(values$change)
  if (n < 3L || length(values$baseline) != n ||
      length(values$treatment) != n || nrow(nuisance_design) != n) {
    stop("Residualized component inputs have incompatible dimensions")
  }
  if (any(!is.finite(unlist(values, use.names = FALSE))) ||
      any(!is.finite(nuisance_design))) {
    stop("Residualized component inputs must be finite")
  }

  nuisance <- cbind(
    nuisance_design,
    analyte_baseline = values$baseline
  )
  nuisance_qr <- qr(nuisance)
  residualized_change <- qr.resid(nuisance_qr, values$change)
  residualized_treatment <- qr.resid(nuisance_qr, values$treatment)
  denominator <- sum(residualized_treatment^2)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    stop("Treatment has no residual variation after nuisance adjustment")
  }

  effect <- sum(residualized_treatment * residualized_change) / denominator
  model_residual <- residualized_change - residualized_treatment * effect
  influence_hc0 <- residualized_treatment * model_residual / denominator
  full_rank <- qr(cbind(nuisance, treatment = values$treatment))$rank
  residual_df <- n - full_rank
  if (residual_df <= 0L) stop("Adjusted model has no residual degrees of freedom")
  hc1_scale <- sqrt(n / residual_df)

  list(
    effect = effect,
    residualized_change = residualized_change,
    residualized_treatment = residualized_treatment,
    model_residual = model_residual,
    influence_hc0 = influence_hc0,
    influence_hc1 = influence_hc0 * hc1_scale,
    n_obs = n,
    model_rank = full_rank,
    residual_df = residual_df,
    treatment_information = denominator,
    hc1_scale = hc1_scale
  )
}

.atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not save ", path)
  invisible(path)
}

run_dunedinpace_residualized_changes <- function(
    omics_raw_path,
    pheno_raw_path,
    control_pc_path,
    cell_pc_path,
    output_dir,
    log = message) {
  started <- proc.time()[["elapsed"]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  log("loading raw CALERIE DNAm inputs")
  omics_raw <- readRDS(omics_raw_path)
  pheno_raw <- readRDS(pheno_raw_path)
  control_pc <- read.csv(control_pc_path, stringsAsFactors = FALSE)
  cell_pcs <- read.csv(cell_pc_path, stringsAsFactors = FALSE)

  stopifnot(
    !anyDuplicated(rownames(omics_raw)),
    !anyDuplicated(colnames(omics_raw)),
    !anyNA(omics_raw),
    !anyDuplicated(pheno_raw$Barcode),
    !anyDuplicated(control_pc$filenames),
    !anyDuplicated(cell_pcs$SAMPLE_ID)
  )

  pheno_columns <- c(
    "Participant_ID", "Time_Point", "Barcode", "fu", "CR", "deidsite",
    "agebl", "female", "bmistrat", "snppc1.x", "snppc2.x", "snppc3.x"
  )
  stopifnot(all(pheno_columns %in% names(pheno_raw)))
  pheno <- pheno_raw[, pheno_columns, drop = FALSE]
  pheno <- pheno[
    !duplicated(paste(pheno$Participant_ID, pheno$Time_Point, sep = "\r")),
    , drop = FALSE
  ]

  pc_columns <- paste0("PC", 1:3)
  cell_pc_columns <- paste0("cell_PC", 1:4)
  stopifnot(
    all(pc_columns %in% names(control_pc)),
    all(cell_pc_columns %in% names(cell_pcs))
  )
  pc_index <- match(pheno$Barcode, control_pc$filenames)
  cell_pc_index <- match(pheno$Barcode, cell_pcs$SAMPLE_ID)
  pheno[pc_columns] <- control_pc[pc_index, pc_columns, drop = FALSE]
  pheno[cell_pc_columns] <- cell_pcs[
    cell_pc_index, cell_pc_columns, drop = FALSE
  ]

  scaled_columns <- c("agebl", "snppc1.x", "snppc2.x", "snppc3.x")
  pheno[scaled_columns] <- lapply(
    pheno[scaled_columns], function(value) as.numeric(scale(as.numeric(value)))
  )
  as_base_factor <- function(value) {
    if (is.numeric(value)) factor(as.numeric(value)) else factor(value)
  }
  pheno$fu <- factor(as.numeric(pheno$fu), levels = c(0, 1, 2))
  pheno$female <- as_base_factor(pheno$female)
  pheno$CR <- as_base_factor(pheno$CR)
  pheno$deidsite <- as_base_factor(pheno$deidsite)
  pheno$bmistrat <- as_base_factor(pheno$bmistrat)
  names(pheno)[match(
    c("Barcode", "Participant_ID", "fu", "female", "CR"), names(pheno)
  )] <- c("SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP")

  covariates <- c(
    "snppc1.x", "snppc2.x", "snppc3.x",
    "agebl", "deidsite", "bmistrat",
    paste0("cell_PC", 1:4), paste0("PC", 1:3)
  )
  required_pheno <- c(
    "SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP",
    covariates
  )
  stopifnot(all(required_pheno %in% names(pheno)))
  pheno <- pheno[
    pheno$SAMPLE_ID %in% colnames(omics_raw) &
      complete.cases(pheno[, required_pheno, drop = FALSE]),
    , drop = FALSE
  ]
  stopifnot(nrow(pheno) > 0L, !anyDuplicated(pheno$SAMPLE_ID))
  fu_numeric <- as.integer(as.character(pheno$FU))
  complete_subjects <- intersect(
    unique(pheno$SUBJECT_ID[fu_numeric == 0L]),
    unique(pheno$SUBJECT_ID[fu_numeric > 0L])
  )
  pheno <- pheno[pheno$SUBJECT_ID %in% complete_subjects, , drop = FALSE]
  if (!nrow(pheno)) {
    stop("No subjects have both a retained baseline and follow-up")
  }

  weights_dir <- file.path(output_dir, "dunedinpace_raw_beta_jackknife_v1")
  weights_path <- file.path(weights_dir, "sysdata.rda")
  if (!file.exists(weights_path)) {
    stop("Pinned DunedinPACE weights not found: ", weights_path)
  }
  reference <- new.env(parent = emptyenv())
  load(weights_path, envir = reference)
  probes_all <- reference$mPACE_Models$model_probes$DunedinPACE
  weights_all <- reference$mPACE_Models$model_weights$DunedinPACE
  if (length(probes_all) != 173L || length(weights_all) != 173L ||
      anyDuplicated(probes_all) || anyNA(weights_all)) {
    stop("Pinned DunedinPACE weights have an unexpected structure")
  }
  weight_table <- data.frame(
    ANALYTE_NAME = as.character(probes_all),
    WEIGHT = as.numeric(weights_all),
    stringsAsFactors = FALSE
  )
  matched <- weight_table[
    weight_table$ANALYTE_NAME %in% rownames(omics_raw), , drop = FALSE
  ]
  if (!nrow(matched)) stop("No DunedinPACE probes matched the raw beta matrix")
  probe_values <- as.matrix(
    omics_raw[match(matched$ANALYTE_NAME, rownames(omics_raw)), , drop = FALSE]
  )
  storage.mode(probe_values) <- "double"
  if (any(!is.finite(probe_values)) || min(probe_values) < -1e-8 ||
      max(probe_values) > 1 + 1e-8) {
    stop("Matched DunedinPACE inputs are not finite raw beta values in [0, 1]")
  }
  rownames(probe_values) <- matched$ANALYTE_NAME
  rm(omics_raw)
  invisible(gc())

  fu_numeric <- as.integer(as.character(pheno$FU))
  baseline <- pheno[fu_numeric == 0L, , drop = FALSE]
  followup <- pheno[fu_numeric > 0L, , drop = FALSE]
  if (anyDuplicated(baseline$SUBJECT_ID) ||
      anyDuplicated(paste(followup$SUBJECT_ID, followup$FU, sep = "\r"))) {
    stop("Expected one baseline and at most one observation per subject/follow-up")
  }
  subject_order <- unique(as.character(pheno$SUBJECT_ID))
  baseline_samples <- stats::setNames(
    as.character(baseline$SAMPLE_ID), as.character(baseline$SUBJECT_ID)
  )
  fu_order <- sort(unique(as.integer(as.character(followup$FU))))
  probes <- matched$ANALYTE_NAME
  component_index <- expand.grid(
    ANALYTE_NAME = probes,
    FU = fu_order,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  component_names <- paste0(
    component_index$ANALYTE_NAME, "::FU", component_index$FU
  )
  n_components <- nrow(component_index)
  n_subjects <- length(subject_order)

  residualized_changes <- matrix(
    NA_real_, n_components, n_subjects,
    dimnames = list(component_names, subject_order)
  )
  residualized_treatment <- residualized_changes
  model_residuals <- residualized_changes
  change_factor <- matrix(
    0, n_components, n_subjects,
    dimnames = list(component_names, subject_order)
  )
  influence_hc0 <- change_factor
  influence_hc1 <- change_factor
  component_summary <- data.frame(
    component_index,
    EFFECT_SIZE = NA_real_,
    N_OBS = NA_integer_,
    MODEL_RANK = NA_integer_,
    RESIDUAL_DF = NA_integer_,
    TREATMENT_INFORMATION = NA_real_,
    HC1_SCALE = NA_real_,
    stringsAsFactors = FALSE
  )

  nuisance_formula <- stats::reformulate(c("FEMALE", covariates))
  treatment_levels <- levels(droplevels(pheno$TREATMENT_GROUP))
  if (length(treatment_levels) != 2L) {
    stop("Expected exactly two treatment-group levels")
  }

  for (fu in fu_order) {
    rows <- followup[as.integer(as.character(followup$FU)) == fu, , drop = FALSE]
    rows <- rows[order(match(as.character(rows$SUBJECT_ID), subject_order)), , drop = FALSE]
    row_subjects <- as.character(rows$SUBJECT_ID)
    subject_columns <- match(row_subjects, subject_order)
    baseline_ids <- unname(baseline_samples[row_subjects])
    if (anyNA(baseline_ids)) stop("A follow-up subject has no retained baseline")

    nuisance_design <- stats::model.matrix(nuisance_formula, data = rows)
    treatment <- as.numeric(
      as.character(rows$TREATMENT_GROUP) == treatment_levels[[2L]]
    )
    followup_values <- probe_values[
      probes, as.character(rows$SAMPLE_ID), drop = FALSE
    ]
    baseline_values <- probe_values[probes, baseline_ids, drop = FALSE]
    changes <- followup_values - baseline_values
    component_rows <- which(component_index$FU == fu)

    log(paste0("FU", fu, ": residualizing ", length(probes),
               " probes across ", nrow(rows), " subjects"))
    for (probe_index in seq_along(probes)) {
      component_row <- component_rows[[probe_index]]
      fit <- .residualized_component(
        change = changes[probe_index, ],
        baseline = baseline_values[probe_index, ],
        treatment = treatment,
        nuisance_design = nuisance_design
      )
      residualized_changes[component_row, subject_columns] <-
        fit$residualized_change
      residualized_treatment[component_row, subject_columns] <-
        fit$residualized_treatment
      model_residuals[component_row, subject_columns] <- fit$model_residual
      change_factor[component_row, subject_columns] <-
        fit$residualized_change / sqrt(fit$n_obs - 1L)
      influence_hc0[component_row, subject_columns] <- fit$influence_hc0
      influence_hc1[component_row, subject_columns] <- fit$influence_hc1
      component_summary$EFFECT_SIZE[[component_row]] <- fit$effect
      component_summary$N_OBS[[component_row]] <- fit$n_obs
      component_summary$MODEL_RANK[[component_row]] <- fit$model_rank
      component_summary$RESIDUAL_DF[[component_row]] <- fit$residual_df
      component_summary$TREATMENT_INFORMATION[[component_row]] <-
        fit$treatment_information
      component_summary$HC1_SCALE[[component_row]] <- fit$hc1_scale
    }
  }

  weight_matrix <- matrix(
    0, nrow = length(fu_order), ncol = n_components,
    dimnames = list(paste0("FU", fu_order), component_names)
  )
  for (fu_index in seq_along(fu_order)) {
    cells <- which(component_index$FU == fu_order[[fu_index]])
    weight_matrix[fu_index, cells] <- matched$WEIGHT[
      match(component_index$ANALYTE_NAME[cells], matched$ANALYTE_NAME)
    ]
  }
  score_effects <- as.numeric(
    weight_matrix %*% component_summary$EFFECT_SIZE
  )
  names(score_effects) <- rownames(weight_matrix)
  score_influence_hc0 <- weight_matrix %*% influence_hc0
  score_influence_hc1 <- weight_matrix %*% influence_hc1

  input_paths <- c(
    omics_raw_path, pheno_raw_path, control_pc_path, cell_pc_path, weights_path
  )
  input_info <- file.info(input_paths)
  input_fingerprint <- data.frame(
    PATH = normalizePath(input_paths),
    SIZE = input_info$size,
    MODIFIED = input_info$mtime,
    stringsAsFactors = FALSE
  )
  artifact <- list(
    schema_version = "dunedinpace_residualized_changes_v1",
    scope = "all",
    interpretation = paste(
      "Separate-FU adjusted linear-model approximation to covariance among",
      "baseline-adjusted treatment-effect estimates; validate against the",
      "exact pooled mixed-model delete-one-subject jackknife."
    ),
    treatment_contrast = paste(
      treatment_levels[[2L]], "minus", treatment_levels[[1L]]
    ),
    normalization = "none",
    input_fingerprint = input_fingerprint,
    covariates = covariates,
    coverage = transform(
      weight_table,
      PRESENT_RAW_BETA = ANALYTE_NAME %in% probes
    ),
    subject_order = subject_order,
    component_index = component_index,
    component_summary = component_summary,
    residualized_changes = residualized_changes,
    residualized_treatment = residualized_treatment,
    model_residuals = model_residuals,
    residualized_change_factor = change_factor,
    residualized_change_covariance = tcrossprod(change_factor),
    effect_influence_factor_hc0 = influence_hc0,
    effect_influence_factor_hc1 = influence_hc1,
    effect_covariance_hc0 = tcrossprod(influence_hc0),
    effect_covariance_hc1 = tcrossprod(influence_hc1),
    score = list(
      effects = score_effects,
      weight_matrix = weight_matrix,
      influence_hc0 = score_influence_hc0,
      influence_hc1 = score_influence_hc1,
      covariance_hc0 = tcrossprod(score_influence_hc0),
      covariance_hc1 = tcrossprod(score_influence_hc1)
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )

  path <- file.path(
    output_dir, "dunedinpace_residualized_changes_v1.rds"
  )
  .atomic_save_rds(artifact, path)
  log(paste0(
    "saved ", path, "; ", length(probes), "/173 probes; ",
    n_subjects, " subjects; elapsed ",
    round(artifact$elapsed_seconds / 60, 1), " min"
  ))
  invisible(artifact)
}

if (sys.nframe() == 0L) {
  omics_raw_path <- path.expand(
    "~/FAST/Data/CALERIE/Raw/DNAm/GRSet_fully_filtered_bmiq_chunk.rds"
  )
  pheno_raw_path <- path.expand(
    "~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds"
  )
  control_pc_path <- path.expand(
    "~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_control_pcs_rgset_goodsamples.csv"
  )
  cell_pc_path <- path.expand(
    "~/FAST/Data/CALERIE/Raw/DNAm/Cell_PCs.csv"
  )
  pipeline_repo <- path.expand("~/FAST/GitHub/track-1.1.1")
  output_dir <- path.expand("~/FAST/Outputs/1.1.1")
  log_path <- file.path(
    output_dir, "dunedinpace_residualized_changes_v1.log"
  )
  log_status <- function(text) {
    line <- paste0(
      "[residualized-change] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | ", text
    )
    message(line)
    cat(line, "\n", file = log_path, append = TRUE, sep = "")
  }

  run_dunedinpace_residualized_changes(
    omics_raw_path = omics_raw_path,
    pheno_raw_path = pheno_raw_path,
    control_pc_path = control_pc_path,
    cell_pc_path = cell_pc_path,
    output_dir = output_dir,
    log = log_status
  )
}
