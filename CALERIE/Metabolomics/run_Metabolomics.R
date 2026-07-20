library(dplyr)
library(forcats)
library(readxl)
library(tibble)
library(tidyr)

# -----------------------------
# Paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/CALERIE/Metabolomics/data/OBT_log.csv")
pheno_raw_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_clinical_w_omics_crosswalks.csv")
genetic_pcs_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_GeneticPCs_20200727.xlsx")

pipeline_repo <- path.expand("~/CALERIE/repos/track-1.1.1-single-analyte")
out_dir <- path.expand("~/CALERIE/Metabolomics")

omics_type <- "Metabolomics"
results_path <- file.path(out_dir, "results_Metabolomics.rds")
reports_path <- file.path(out_dir, "reports_Metabolomics.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_Metabolomics.log")

cat("START run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Read inputs
# -----------------------------
omics_raw <- read.csv(
  omics_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "latin1"
)

pheno_raw <- read.csv(
  pheno_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

PCs <- read_excel(genetic_pcs_path)

# -----------------------------
# Input QC
# -----------------------------
stopifnot("CLIENT_IDENTIFIER" %in% colnames(omics_raw))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!any(duplicated(omics_raw$CLIENT_IDENTIFIER)))

metabolite_cols <- names(omics_raw)[25:ncol(omics_raw)]
stopifnot(length(metabolite_cols) > 0)
stopifnot(!any(is.na(omics_raw[, metabolite_cols])))

stopifnot("ID" %in% colnames(PCs))
stopifnot(!any(duplicated(PCs$ID)))

# -----------------------------
# Build phenotype table
# -----------------------------
pheno <- pheno_raw |>
  select(DEID, CR, deidsite, agebl, female, bmistrat, fumo, metabolomics_barcode) |>
  left_join(
    PCs |>
      select(ID, PC1, PC2, PC3),
    by = c("DEID" = "ID")
  ) |>
  distinct(metabolomics_barcode, .keep_all = TRUE) |>
  rename(
    SAMPLE_ID = metabolomics_barcode,
    SUBJECT_ID = DEID,
    FU = fumo,
    FEMALE = female,
    TREATMENT_GROUP = CR
  ) |>
  mutate(
    agebl = as.numeric(scale(agebl)),
    PC1 = as.numeric(scale(PC1)),
    PC2 = as.numeric(scale(PC2)),
    PC3 = as.numeric(scale(PC3))
  ) |>
  mutate(
    FEMALE = as_factor(FEMALE),
    TREATMENT_GROUP = as_factor(TREATMENT_GROUP),
    deidsite = as_factor(deidsite),
    bmistrat = as_factor(bmistrat)
  ) |>
  mutate(
    FU = case_when(
      FU == 0 ~ 0,
      FU == 6 ~ 1,
      FU == 12 ~ 2,
      FU == 24 ~ 3,
      TRUE ~ NA_real_
    ),
    FU = factor(FU, levels = c(0, 1, 2, 3))
  )

covariates <- c(
  c("PC1", "PC2", "PC3"),
  c("agebl", "deidsite", "bmistrat")
)

required_pheno_cols <- c(
  "SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP",
  covariates
)

stopifnot(all(required_pheno_cols %in% colnames(pheno)))

pheno <- pheno |>
  filter(SAMPLE_ID %in% omics_raw$CLIENT_IDENTIFIER) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x)))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% omics_raw$CLIENT_IDENTIFIER))

# -----------------------------
# Build omics table
# -----------------------------
omics <- omics_raw |>
  select(CLIENT_IDENTIFIER, all_of(metabolite_cols)) |>
  column_to_rownames("CLIENT_IDENTIFIER") |>
  t() |>
  as.data.frame(check.names = FALSE) |>
  rownames_to_column("ANALYTE_NAME")

stopifnot(!any(duplicated(omics$ANALYTE_NAME)))
stopifnot(all(pheno$SAMPLE_ID %in% colnames(omics)))

# -----------------------------
# Load pipeline functions
# -----------------------------
source(file.path(pipeline_repo, "main.R"), chdir = TRUE)

stopifnot(
  exists("FAST_omics_WAS"),
  exists("FAST_omics_WAS_reports")
)

# -----------------------------
# Analysis
# -----------------------------
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates
)

saveRDS(results, file = results_path)

cat("DONE analysis: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Reports
# -----------------------------
cat("START report run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates
)

saveRDS(reports, file = reports_path)

cat("DONE reports: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(pipeline_repo, "plotting_helpers.R"), chdir = TRUE)

fig_change <- file.path(out_dir, "Figures", "change")
fig_level <- file.path(out_dir, "Figures", "level")

dir.create(fig_change, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_level, recursive = TRUE, showWarnings = FALSE)

generate_all_plots(
  results,
  figures_dir = fig_change,
  analysis = "analysis_change"
)

generate_all_plots(
  results,
  figures_dir = fig_level,
  analysis = "analysis_level"
)