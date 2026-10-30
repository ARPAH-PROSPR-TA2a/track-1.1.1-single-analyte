# Follow-up codes remain model/checkpoint identities; labels are output metadata.
.validate_fu_labels <- function(fu_labels, fu_codes) {
  if (is.null(fu_labels)) return(NULL)
  if (!is.character(fu_labels) || !is.null(dim(fu_labels)) ||
      is.null(names(fu_labels)) || anyNA(names(fu_labels)) ||
      anyDuplicated(names(fu_labels)) || !setequal(names(fu_labels), fu_codes)) {
    stop("fu_labels must be a named character vector with unique code names matching every observed FU code exactly, including baseline 0.")
  }
  if (anyNA(fu_labels) || any(!nzchar(trimws(fu_labels))) || anyDuplicated(fu_labels)) {
    stop("fu_labels must contain unique, non-missing, non-blank labels.")
  }
  fu_labels[fu_codes]
}

.attach_fu_labels <- function(table, fu_labels) {
  if (is.null(table) || is.null(fu_labels)) return(table)
  table$FU_LABEL <- unname(fu_labels[as.character(table$FU)])
  columns <- setdiff(names(table), "FU_LABEL")
  table[, append(columns, "FU_LABEL", after = match("FU", columns)), drop = FALSE]
}

.label_fu_coefficients <- function(coefficients, pheno, fu_labels,
                                   additional_covariates = NULL) {
  if (is.null(coefficients) || nrow(coefficients) == 0L || is.null(fu_labels)) {
    return(coefficients)
  }
  followups <- pheno[pheno$FU != "0", , drop = FALSE]
  visits <- levels(droplevels(followups$FU))
  if (length(visits) < 2L) return(coefficients)

  # Obtain R's exact main/interaction term names, without refitting. This also
  # preserves contrast indices under sum/Helmert coding, where labels do not
  # name individual coefficients. A small reference grid suffices.
  grid <- expand.grid(
    TREATMENT_GROUP = levels(pheno$TREATMENT_GROUP), FU = visits,
    stringsAsFactors = FALSE
  )
  grid$TREATMENT_GROUP <- factor(grid$TREATMENT_GROUP, levels = levels(pheno$TREATMENT_GROUP))
  grid$FU <- factor(grid$FU, levels = visits)
  original <- colnames(stats::model.matrix(~ TREATMENT_GROUP * FU, grid))
  levels(grid$FU) <- unname(fu_labels[visits])
  labeled <- colnames(stats::model.matrix(~ TREATMENT_GROUP * FU, grid))
  changed <- original != labeled
  mapping <- setNames(labeled[changed], original[changed])
  if (!length(mapping)) return(coefficients)

  # A covariate could itself generate a column such as FU2 or FU24mo. Do not
  # reinterpret it as a visit coefficient or produce ambiguous output names.
  if (length(additional_covariates)) {
    covariate_names <- colnames(stats::model.matrix(
      stats::reformulate(additional_covariates), followups
    ))
    collisions <- intersect(c(names(mapping), unname(mapping)), covariate_names)
    if (length(collisions)) {
      stop("fu_labels create ambiguous coefficient names with additional covariates: ",
           paste(collisions, collapse = ", "), ". Choose different labels or covariate names.")
    }
  }
  terms <- unique(coefficients$COEFFICIENT)
  renamed <- terms
  matched <- match(terms, names(mapping))
  renamed[!is.na(matched)] <- unname(mapping[matched[!is.na(matched)]])
  if (anyDuplicated(renamed)) {
    stop("fu_labels create duplicate coefficient names. Choose different labels.")
  }
  coefficients$COEFFICIENT <- renamed[match(coefficients$COEFFICIENT, terms)]
  coefficients
}
