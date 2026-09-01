.DPACE_REVISION <- "4b569983543e51d1022aecec9a25e694bb3a336a"
.DPACE_SHA256 <-
  "356f3c311918a9c45233e43ae57aeff202f075d970bab59c863bd67aebbfd781"
.DPACE_URL <- paste0(
  "https://raw.githubusercontent.com/danbelsky/DunedinPACE/",
  .DPACE_REVISION,
  "/R/sysdata.rda"
)

.dpace_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not save ", path)
  invisible(path)
}

.dpace_weights <- function(output_dir, log) {
  path <- file.path(output_dir, "sysdata.rda")
  if (!file.exists(path)) {
    log("DunedinPACE: downloading pinned weights")
    temporary <- paste0(path, ".download-", Sys.getpid())
    utils::download.file(.DPACE_URL, temporary, mode = "wb", quiet = FALSE)
    if (digest::digest(temporary, "sha256", file = TRUE) != .DPACE_SHA256) {
      stop("Downloaded DunedinPACE weights failed the SHA-256 check")
    }
    if (!file.rename(temporary, path)) stop("Could not cache DunedinPACE weights")
  }
  if (digest::digest(path, "sha256", file = TRUE) != .DPACE_SHA256) {
    stop("Cached DunedinPACE weights failed the SHA-256 check: ", path)
  }

  reference <- new.env(parent = emptyenv())
  load(path, envir = reference)
  probes <- reference$mPACE_Models$model_probes$DunedinPACE
  weights <- reference$mPACE_Models$model_weights$DunedinPACE
  if (length(probes) != 173L || length(weights) != 173L ||
      anyDuplicated(probes) || anyNA(weights) ||
      !identical(as.character(probes), names(weights))) {
    stop("Pinned DunedinPACE weights have an unexpected structure")
  }
  data.frame(
    ANALYTE_NAME = as.character(probes),
    WEIGHT = as.numeric(weights),
    stringsAsFactors = FALSE
  )
}

.dpace_extract <- function(results, probes, fu_order) {
  table <- results$treatment_effects
  expected <- as.vector(outer(
    probes, fu_order, function(probe, fu) paste(probe, fu, sep = "\r")
  ))
  observed <- paste(table$ANALYTE_NAME, table$FU, sep = "\r")
  index <- match(expected, observed)
  if (anyNA(index)) {
    stop("FAST omitted ", sum(is.na(index)), " required probe/FU result(s)")
  }
  ordered <- table[index, , drop = FALSE]
  if (any(!is.finite(ordered$EFFECT_SIZE)) || any(!is.finite(ordered$SE))) {
    stop("FAST returned non-finite DunedinPACE effects or standard errors")
  }
  dimensions <- list(probes, paste0("FU", fu_order))
  list(
    effects = matrix(ordered$EFFECT_SIZE, nrow = length(probes),
                     dimnames = dimensions),
    standard_errors = matrix(ordered$SE, nrow = length(probes),
                             dimnames = dimensions)
  )
}

.dpace_fit <- function(pheno, omics, covariates, probes, fu_order) {
  results <- suppressMessages(.perform_analysis(
    pheno_df = pheno,
    omics_df = omics,
    omics_type = "other",
    mixed_effects = TRUE,
    additional_covariates = covariates,
    response_type = "change",
    checkpoint_dir = NULL,
    checkpoint_batch_size = nrow(omics),
    verbose = FALSE
  ))
  .dpace_extract(results, probes, fu_order)
}

.dpace_assemble <- function(checkpoints, reference, probes, weights, fu_order) {
  n_subjects <- length(checkpoints)
  n_components <- length(probes) * length(fu_order)
  leave_one_out <- vapply(
    checkpoints, function(checkpoint) as.vector(checkpoint$effects),
    numeric(n_components)
  )
  component_index <- expand.grid(
    ANALYTE_NAME = probes, FU = fu_order,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  component_names <- paste0(component_index$ANALYTE_NAME, "::FU", component_index$FU)
  rownames(leave_one_out) <- component_names
  colnames(leave_one_out) <- vapply(
    checkpoints, `[[`, character(1), "subject_id"
  )

  jackknife_factor <- sqrt((n_subjects - 1) / n_subjects) * sweep(
    leave_one_out, 1L, rowMeans(leave_one_out), "-"
  )
  covariance <- tcrossprod(jackknife_factor)
  decomposition <- svd(jackknife_factor, nv = 0L)
  tolerance <- max(dim(jackknife_factor)) * max(decomposition$d) *
    .Machine$double.eps
  rank <- sum(decomposition$d > tolerance)
  covariance_factor <- sweep(
    decomposition$u[, seq_len(rank), drop = FALSE], 2L,
    decomposition$d[seq_len(rank)], "*"
  )
  rownames(covariance_factor) <- component_names
  colnames(covariance_factor) <- paste0("factor_", seq_len(rank))

  weight_matrix <- matrix(0, nrow = length(fu_order), ncol = n_components)
  for (fu_index in seq_along(fu_order)) {
    cells <- (fu_index - 1L) * length(probes) + seq_along(probes)
    weight_matrix[fu_index, cells] <- weights
  }
  rownames(weight_matrix) <- paste0("FU", fu_order)
  colnames(weight_matrix) <- component_names
  score_effects <- as.numeric(weight_matrix %*% as.vector(reference$effects))
  names(score_effects) <- rownames(weight_matrix)
  score_covariance <- weight_matrix %*% covariance %*% t(weight_matrix)
  legacy_covariance <- weight_matrix %*%
    diag(as.vector(reference$standard_errors)^2, nrow = n_components) %*%
    t(weight_matrix)
  diagonal_jackknife_covariance <- weight_matrix %*%
    diag(diag(covariance), nrow = n_components) %*% t(weight_matrix)

  score_contrasts <- diag(length(fu_order))
  dimnames(score_contrasts) <- list(
    paste0("FU", fu_order), paste0("FU", fu_order)
  )
  if (length(fu_order) > 1L) {
    difference <- numeric(length(fu_order))
    difference[c(1L, length(fu_order))] <- c(-1, 1)
    score_contrasts <- rbind(score_contrasts, difference)
    rownames(score_contrasts)[nrow(score_contrasts)] <- paste0(
      "FU", tail(fu_order, 1L), "-FU", fu_order[1L]
    )
  }
  contrast_effects <- as.numeric(score_contrasts %*% score_effects)
  projected_se <- function(value) sqrt(pmax(
    0, diag(score_contrasts %*% value %*% t(score_contrasts))
  ))
  se_legacy <- projected_se(legacy_covariance)
  se_diagonal <- projected_se(diagonal_jackknife_covariance)
  se_full <- projected_se(score_covariance)
  t_statistic <- contrast_effects / se_full
  inference <- data.frame(
    CONTRAST = rownames(score_contrasts),
    EFFECT_SIZE = contrast_effects,
    SE_LEGACY_DIAGONAL = se_legacy,
    SE_JACKKNIFE_DIAGONAL = se_diagonal,
    SE_JACKKNIFE_FULL = se_full,
    T_STAT = t_statistic,
    DF = n_subjects - 1L,
    P_VALUE_JACKKNIFE_T = 2 * stats::pt(
      -abs(t_statistic), df = n_subjects - 1L
    ),
    P_VALUE_JACKKNIFE_Z = 2 * stats::pnorm(-abs(t_statistic)),
    stringsAsFactors = FALSE
  )

  component_index$WEIGHT <- weights[
    match(component_index$ANALYTE_NAME, probes)
  ]
  list(
    n_subjects = n_subjects,
    subject_order = colnames(leave_one_out),
    fu_order = fu_order,
    component_index = component_index,
    full_probe_effects = reference$effects,
    legacy_probe_standard_errors = reference$standard_errors,
    leave_one_out_effects = leave_one_out,
    jackknife_factor = jackknife_factor,
    covariance = covariance,
    svd_factor = list(
      factor = covariance_factor,
      singular_values = decomposition$d[seq_len(rank)],
      numerical_rank = rank,
      usage = "Use factor[, 1:r, drop = FALSE] for a rank-r approximation"
    ),
    score = list(
      effects = score_effects,
      weight_matrix = weight_matrix,
      contrast_matrix = score_contrasts,
      covariance = score_covariance,
      legacy_diagonal_covariance = legacy_covariance,
      jackknife_diagonal_covariance = diagonal_jackknife_covariance,
      inference = inference,
      p_value_note = "Approximate t calibration uses df = n_subjects - 1"
    )
  )
}

.dpace_run_stratum <- function(stratum, pheno, omics, reference, covariates,
                                probes, weights, signature, output_dir,
                                n_cores, log) {
  fu_order <- sort(unique(as.integer(as.character(pheno$FU))))
  fu_order <- fu_order[fu_order > 0L]
  subjects <- unique(as.character(pheno$SUBJECT_ID))
  checkpoint_dir <- file.path(output_dir, "checkpoints", stratum)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(checkpoint_dir, sprintf("loo_%04d.rds", seq_along(subjects)))
  checkpoints <- vector("list", length(subjects))

  for (index in which(file.exists(paths))) {
    checkpoint <- readRDS(paths[index])
    valid <- identical(checkpoint$signature, signature) &&
      identical(checkpoint$stratum, stratum) &&
      identical(checkpoint$index, index) &&
      identical(checkpoint$subject_id, subjects[index]) &&
      identical(dimnames(checkpoint$effects), dimnames(reference$effects))
    if (!isTRUE(valid) || any(!is.finite(checkpoint$effects))) {
      stop("Incompatible DunedinPACE checkpoint: ", paths[index])
    }
    checkpoints[[index]] <- checkpoint
  }

  pending <- which(vapply(checkpoints, is.null, logical(1)))
  cached <- length(subjects) - length(pending)
  log(paste0("DunedinPACE ", stratum, ": ", cached, "/", length(subjects),
             " subject jackknifes cached; ", length(pending), " pending"))
  waves <- split(pending, ceiling(seq_along(pending) / max(1L, n_cores)))
  started <- proc.time()[["elapsed"]]
  completed <- 0L

  for (wave_number in seq_along(waves)) {
    wave <- waves[[wave_number]]
    if (wave_number == 1L) {
      log(paste0("DunedinPACE ", stratum,
                 ": jackknife loop starting; dispatching ", length(wave),
                 " subject deletions across ", min(n_cores, length(wave)),
                 " fork workers"))
    }
    results <- parallel::mclapply(wave, function(index) {
      tryCatch({
        reduced_pheno <- pheno[as.character(pheno$SUBJECT_ID) != subjects[index], ]
        fit <- .dpace_fit(reduced_pheno, omics, covariates, probes, fu_order)
        list(
          signature = signature, stratum = stratum, index = index,
          subject_id = subjects[index], effects = fit$effects,
          worker_pid = Sys.getpid()
        )
      }, error = function(error) {
        list(index = index, worker_pid = Sys.getpid(), error = conditionMessage(error))
      })
    }, mc.cores = min(n_cores, length(wave)))

    if (wave_number == 1L) {
      pids <- sort(unique(vapply(results, `[[`, integer(1), "worker_pid")))
      log(paste0("DunedinPACE ", stratum, ": first wave worker PIDs: ",
                 paste(pids, collapse = ", ")))
    }
    failed <- vapply(results, function(result) !is.null(result$error), logical(1))
    for (result in results[!failed]) {
      .dpace_save(result, paths[result$index])
      checkpoints[[result$index]] <- result
    }
    if (any(failed)) {
      first <- results[[which(failed)[1L]]]
      stop("DunedinPACE subject ", subjects[first$index], " failed: ", first$error)
    }

    completed <- completed + length(results)
    elapsed <- proc.time()[["elapsed"]] - started
    remaining <- length(pending) - completed
    eta_minutes <- if (completed > 0L) elapsed / completed * remaining / 60 else NA_real_
    log(paste0("DunedinPACE ", stratum, ": ", cached + completed, "/",
               length(subjects), " complete; ETA ", round(eta_minutes, 1), " min"))
  }
  .dpace_assemble(checkpoints, reference, probes, weights, fu_order)
}

run_dunedinpace_covariance <- function(pheno, omics, additional_covariates,
                                        output_dir, n_cores, log = message) {
  required <- c("digest", "future", "furrr", "lme4", "lmerTest", "emmeans")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))
  if (!exists(".perform_analysis") || !exists(".validate_pheno") ||
      !exists(".validate_omics")) stop("Source track-1.1.1/main.R first")

  started <- proc.time()[["elapsed"]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  weight_table <- .dpace_weights(output_dir, log)
  weight_table$PRESENT_RAW_BETA <-
    weight_table$ANALYTE_NAME %in% as.character(omics$ANALYTE_NAME)
  matched <- weight_table[weight_table$PRESENT_RAW_BETA, , drop = FALSE]
  if (!nrow(matched)) stop("No DunedinPACE probes matched the raw beta data")
  target <- omics[
    match(matched$ANALYTE_NAME, as.character(omics$ANALYTE_NAME)), , drop = FALSE
  ]
  values <- as.matrix(target[, setdiff(names(target), "ANALYTE_NAME"), drop = FALSE])
  storage.mode(values) <- "double"
  value_range <- range(values)
  if (any(!is.finite(values)) || value_range[1L] < -1e-8 ||
      value_range[2L] > 1 + 1e-8) {
    stop("Matched DunedinPACE inputs are not finite raw beta values in [0, 1]")
  }

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(target, pheno_list, verbose = FALSE)
  signature <- digest::digest(list(
    schema = "dunedinpace_raw_beta_jackknife_v1",
    probes = matched$ANALYTE_NAME,
    weights = matched$WEIGHT,
    covariates = additional_covariates,
    pheno = pheno_list$all,
    omics = omics_list$all
  ), algo = "sha256")
  config_path <- file.path(output_dir, "config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path)$signature, signature)) {
      stop("Existing DunedinPACE checkpoints belong to different inputs")
    }
  } else {
    .dpace_save(list(signature = signature), config_path)
  }
  log(paste0("DunedinPACE: ", nrow(matched), "/173 raw-beta probes matched; range [",
             signif(value_range[1L], 4L), ", ", signif(value_range[2L], 4L), "]"))

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_cores > 1L) {
    future::plan(future::multicore, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }
  references <- list()
  for (stratum in c("all", "male", "female")) {
    if (is.null(pheno_list[[stratum]])) next
    fu_order <- sort(unique(as.integer(as.character(pheno_list[[stratum]]$FU))))
    fu_order <- fu_order[fu_order > 0L]
    log(paste0("DunedinPACE ", stratum, ": full-data reference starting"))
    references[[stratum]] <- .dpace_fit(
      pheno_list[[stratum]], omics_list[[stratum]], additional_covariates,
      matched$ANALYTE_NAME, fu_order
    )
    log(paste0("DunedinPACE ", stratum, ": full-data reference complete"))
  }

  future::plan(future::sequential)
  strata <- list()
  for (stratum in names(references)) {
    strata[[stratum]] <- .dpace_run_stratum(
      stratum, pheno_list[[stratum]], omics_list[[stratum]],
      references[[stratum]], additional_covariates, matched$ANALYTE_NAME,
      matched$WEIGHT, signature, output_dir, n_cores, log
    )
    inference <- strata[[stratum]]$score$inference
    log(paste0("DunedinPACE ", stratum, " inference: ", paste0(
      inference$CONTRAST, " p=", signif(inference$P_VALUE_JACKKNIFE_T, 4L),
      collapse = "; "
    )))
  }

  artifact <- list(
    schema_version = "dunedinpace_raw_beta_jackknife_v1",
    method = "delete-one-subject jackknife",
    normalization = "none",
    weight_renormalization = FALSE,
    score_interpretation = paste(
      "Weighted sum of per-probe baseline-adjusted FAST effects; comparison",
      "with the separately fitted DunedinPACE score is empirical."
    ),
    reference_revision = .DPACE_REVISION,
    reference_sha256 = .DPACE_SHA256,
    signature = signature,
    coverage = weight_table,
    strata = strata,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
  .dpace_save(artifact, file.path(output_dir, "dunedinpace_covariance.rds"))
  log(paste0("DunedinPACE covariance saved; elapsed ",
             round(artifact$elapsed_seconds / 60, 1), " min"))
  artifact
}
