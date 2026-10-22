# Deterministic test measurements and phenotypes; no participant data are read.
# Probe names come from the public lists already bundled with this repository.
.make_synthetic_test_data <- function(n_analytes = 16L, seed = 20261012L) {
  stopifnot(length(n_analytes) == 1L, n_analytes >= 2L,
            n_analytes == as.integer(n_analytes))
  set.seed(seed)

  full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
  filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
  n_filtered <- n_analytes %/% 2L
  analyte_names <- c(
    head(intersect(full_probes, filtered_probes), n_filtered),
    head(setdiff(full_probes, filtered_probes), n_analytes - n_filtered)
  )
  stopifnot(length(analyte_names) == n_analytes, !anyDuplicated(analyte_names))

  # Cross these categories so every sex stratum has both treatment groups and
  # estimable factor covariates. Three replicates give 72 synthetic subjects.
  subjects <- expand.grid(
    TREATMENT_GROUP = 0:1, FEMALE = 0:1,
    ethnic = c("no", "yes"), race3 = c("a", "b", "c"),
    replicate = 1:3, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_subjects <- nrow(subjects)
  subjects$SUBJECT_ID <- sprintf("synthetic_%03d", seq_len(n_subjects))
  subjects$agebl <- runif(n_subjects, 35, 65)
  subjects$mbmi <- runif(n_subjects, 20, 32)
  subjects$TREATMENT_GROUP <- factor(subjects$TREATMENT_GROUP, levels = 0:1)
  subjects$FEMALE <- factor(subjects$FEMALE, levels = 0:1)
  subjects$ethnic <- factor(subjects$ethnic)
  subjects$race3 <- factor(subjects$race3)
  subjects$replicate <- NULL

  pheno <- subjects[rep(seq_len(n_subjects), each = 3L), ]
  fu <- rep(0:2, times = n_subjects)
  pheno$FU <- factor(fu, levels = 0:2)
  pheno$SAMPLE_ID <- paste0(pheno$SUBJECT_ID, "_FU", fu)
  # Visit-specific variation avoids agevis being exactly agebl + FU.
  pheno$agevis <- pheno$agebl + fu + rnorm(nrow(pheno), sd = 0.25)
  row.names(pheno) <- NULL

  baseline <- matrix(rnorm(n_analytes * n_subjects, 0.5, 0.07),
                     nrow = n_analytes)
  subject_change <- matrix(rnorm(n_analytes * n_subjects, sd = 0.04),
                           nrow = n_analytes)
  treatment_effect <- seq(-0.025, 0.025, length.out = n_analytes)
  values <- matrix(NA_real_, nrow = n_analytes, ncol = nrow(pheno))
  for (i in seq_len(n_subjects)) {
    active <- as.integer(as.character(subjects$TREATMENT_GROUP[i]))
    for (visit in 0:2) {
      j <- (i - 1L) * 3L + visit + 1L
      values[, j] <- if (visit == 0L) baseline[, i] else {
        baseline[, i] + subject_change[, i] +
          visit * (0.005 + treatment_effect * active) +
          rnorm(n_analytes, sd = 0.035)
      }
    }
  }
  colnames(values) <- pheno$SAMPLE_ID
  omics <- data.frame(ANALYTE_NAME = analyte_names, values, check.names = FALSE)
  list(pheno = pheno, omics = omics,
       full_probes = full_probes, filtered_probes = filtered_probes)
}
