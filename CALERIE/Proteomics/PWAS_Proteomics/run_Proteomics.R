library(dplyr)
library(forcats)
library(readr)
library(readxl)
library(tibble)
library(tidyr)

# -----------------------------
# Paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_cleaned_log2_soma_matrix.csv")
pheno_raw_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_clinical_w_omics_crosswalks.csv")
genetic_pcs_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_GeneticPCs_20200727.xlsx")

pipeline_repo <- path.expand("~/CALERIE/repos/track-1.1.1-single-analyte")
plotting_repo <- path.expand("~/CALERIE/repos/track-1.1.1-single-analyte")
out_dir <- path.expand("~/CALERIE/calerie_Proteomics")

omics_type <- "Proteomics"
n_cores <- 30
checkpoint_batch_size <- 2000L

results_rds <- file.path(out_dir, "results_full_Proteomics.rds")
reports_rds <- file.path(out_dir, "reports_full_Proteomics.rds")
translation_csv <- file.path(out_dir, "CALERIE_cleaned_protein_translation_table.csv")
annotated_tables_path <- file.path(out_dir, "annotated_treatment_effect_tables.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_full_Proteomics.log")

cat("START run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Read inputs
# -----------------------------
omics_raw <- read.csv(
  omics_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
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
stopifnot("SampleId" %in% colnames(omics_raw))
stopifnot(!any(is.na(omics_raw)))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!any(duplicated(omics_raw$SampleId)))

stopifnot("ID" %in% colnames(PCs))
stopifnot(!any(duplicated(PCs$ID)))

# -----------------------------
# Build phenotype table
# -----------------------------
pheno <- pheno_raw |>
  select(DEID, CR, deidsite, agebl, female, bmistrat, fumo, proteomics_barcode) |>
  left_join(
    PCs |>
      select(ID, PC1, PC2, PC3),
    by = c("DEID" = "ID")
  ) |>
  distinct(proteomics_barcode, .keep_all = TRUE) |>
  rename(
    SAMPLE_ID = proteomics_barcode,
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
      FU == 3 ~ 1,
      FU == 6 ~ 2,
      FU == 12 ~ 3,
      FU == 24 ~ 4
    ),
    FU = factor(FU, levels = c(0, 1, 2, 3, 4))
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
  filter(SAMPLE_ID %in% omics_raw$SampleId) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x)))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% omics_raw$SampleId))

# -----------------------------
# Build omics table
# -----------------------------
omics <- omics_raw |>
  column_to_rownames("SampleId") |>
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
  additional_covariates = covariates,
  n_cores = n_cores,
  checkpoint_dir = file.path(out_dir, "checkpoints"),
  checkpoint_batch_size = checkpoint_batch_size
)

saveRDS(results, file = results_rds)

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

saveRDS(reports, file = reports_rds)

cat("DONE reports: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(plotting_repo, "plotting_helpers.R"), chdir = TRUE)

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

# -----------------------------
# Annotated treatment-effect tables
# -----------------------------
source(file.path(pipeline_repo, "proteomics_translation_helpers.R"), chdir = TRUE)
stopifnot(file.exists(translation_csv))

analysis_types <- c("change", "level")
groups <- c("all", "male", "female")
fu_values <- 1:4
fu_labels <- c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")

annotated_tables <- list()

for (analysis_type in analysis_types) {
  for (group in groups) {
    for (fu in fu_values) {
      table_name <- paste(analysis_type, group, fu_labels[as.character(fu)], sep = "_")
      
      annotated_tables[[table_name]] <- get_annotated_treatment_effects(
        analysis_type = analysis_type,
        group = group,
        fu = fu,
        results_rds = results_rds,
        translation_csv = translation_csv
      )
    }
  }
}

saveRDS(annotated_tables, annotated_tables_path)

cat("Saved ", length(annotated_tables),
    " annotated treatment-effect tables to: ", annotated_tables_path, "\n",
    file = log_file, append = TRUE, sep = "")
cat("DONE annotation: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

cat("Saved", length(annotated_tables),
    "annotated treatment-effect tables to:", annotated_tables_path, "\n")
cat("Table names:\n")
print(names(annotated_tables))