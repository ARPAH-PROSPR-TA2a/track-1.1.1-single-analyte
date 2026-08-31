library(dplyr)
library(forcats)
library(readr)

# -----------------------------
# Paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/GRSet_fully_filtered_bmiq_chunk.rds")
pheno_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds")
control_pc_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_control_pcs_rgset_goodsamples.csv")
cell_count_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_SALAS_cell_counts_chunk.csv")
cell_pcs_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/Cell_PCs.csv")

pipeline_repo <- path.expand("~/FAST/GitHub/track-1.1.1")
out_dir <- path.expand("~/FAST/Outputs/1.1.1")

omics_type <- "DNAm"
n_cores <- 15
checkpoint_batch_size <- 2000L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_full_DNAm.log")
options(track111.progress_log = log_file)

log_status <- function(text) {
  line <- paste0(
    "[1.1.1] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", text
  )
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}

run_started <- proc.time()[["elapsed"]]
log_status("run starting")

# -----------------------------
# Read inputs
# -----------------------------
input_started <- proc.time()[["elapsed"]]
log_status("loading raw inputs (serial)")
omics_raw <- readRDS(omics_raw_path)
pheno_raw <- readRDS(pheno_raw_path)
control_pc <- read.csv(control_pc_path, stringsAsFactors = FALSE)
cell_count <- read.csv(cell_count_path, stringsAsFactors = FALSE)
cell_PCs <- read_csv(cell_pcs_path, show_col_types = FALSE)
log_status(
  paste0(
    "raw inputs loaded: ", nrow(omics_raw), " omics rows, ",
    nrow(pheno_raw), " phenotype rows (",
    round(proc.time()[["elapsed"]] - input_started, 1), " sec)"
  )
)
preparation_started <- proc.time()[["elapsed"]]
log_status("input QC and table preparation starting (serial)")

# -----------------------------
# Input QC
# -----------------------------
stopifnot(!any(duplicated(rownames(omics_raw))))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!any(is.na(omics_raw)))

stopifnot(!any(duplicated(pheno_raw$Barcode)))
stopifnot(!any(duplicated(control_pc$filenames)))
stopifnot(!any(duplicated(cell_count$IDATid)))
stopifnot(!any(duplicated(cell_PCs$SAMPLE_ID)))

# -----------------------------
# Build phenotype table
# -----------------------------
pheno <- pheno_raw |>
  select(
    Participant_ID, Time_Point, Barcode, fu, CR, deidsite, agebl, female,
    bmistrat, snppc1.x, snppc2.x, snppc3.x
  ) |>
  left_join(
    cell_count |> select(-any_of("chunk")),
    by = c("Barcode" = "IDATid")
  ) |>
  left_join(
    control_pc |> select(filenames, paste0("PC", 1:20)),
    by = c("Barcode" = "filenames")
  ) |>
  group_by(Participant_ID, Time_Point) |>
  slice(1) |>
  ungroup() |>
  mutate(
    agebl = as.numeric(scale(agebl)),
    snppc1.x = as.numeric(scale(snppc1.x)),
    snppc2.x = as.numeric(scale(snppc2.x)),
    snppc3.x = as.numeric(scale(snppc3.x))
  ) |>
  mutate(
    fu = factor(haven::zap_labels(fu), levels = c(0, 1, 2)),
    female = as_factor(female),
    CR = as_factor(CR),
    deidsite = as_factor(deidsite),
    bmistrat = as_factor(bmistrat)
  ) |>
  rename(
    SAMPLE_ID = Barcode,
    SUBJECT_ID = Participant_ID,
    FU = fu,
    FEMALE = female,
    TREATMENT_GROUP = CR
  ) |>
  select(
    -Time_Point,
    -any_of(c(
      "Bas", "Bmem", "Bnv",
      "CD4mem", "CD4nv",
      "CD8mem", "CD8nv",
      "Eos", "Mono", "Neu", "NK", "Treg",
      paste0("PC", 4:20)
    ))
  ) |>
  left_join(
    cell_PCs |> select(SAMPLE_ID, any_of(paste0("cell_PC", 1:4))),
    by = "SAMPLE_ID"
  )

covariates <- c(
  c("snppc1.x", "snppc2.x", "snppc3.x"),
  c("agebl", "deidsite", "bmistrat"),
  paste0("cell_PC", 1:4),
  paste0("PC", 1:3)
)

required_pheno_cols <- c(
  "SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP",
  covariates
)

stopifnot(all(required_pheno_cols %in% colnames(pheno)))

pheno <- pheno |>
  filter(SAMPLE_ID %in% colnames(omics_raw)) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x)))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% colnames(omics_raw)))

# -----------------------------
# Build omics table
# -----------------------------
omics <- data.frame(
  ANALYTE_NAME = rownames(omics_raw),
  data.frame(omics_raw, check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(!any(duplicated(omics$ANALYTE_NAME)))
stopifnot(all(pheno$SAMPLE_ID %in% colnames(omics)))

# -----------------------------
# Load pipeline functions
# -----------------------------
source(file.path(pipeline_repo, "main.R"), chdir = TRUE)

stopifnot(
  exists("FAST_omics_WAS"),
  exists("FAST_omics_WAS_reports"),
  exists("generate_all_plots") || file.exists(file.path(pipeline_repo, "plotting_helpers.R"))
)

log_status(
  paste0(
    "inputs ready: ", nrow(omics), " analytes, ",
    nrow(pheno), " matched phenotype rows (",
    round(proc.time()[["elapsed"]] - preparation_started, 1), " sec)"
  )
)

# -----------------------------
# Analysis
# -----------------------------
log_status(
  paste0(
    "analysis call starting: n_cores=", n_cores,
    ", checkpoint batch size=", checkpoint_batch_size,
    "; parallel work occurs inside uncached batches"
  )
)
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  n_cores = n_cores,
  checkpoint_dir = file.path(out_dir, "checkpoints"),
  checkpoint_batch_size = checkpoint_batch_size,
  verbose = TRUE
)

saveRDS(
  results,
  file = file.path(out_dir, "DNAm_betas_1.1.1.rds")
)

log_status("analysis complete; results saved")

# -----------------------------
# Reports
# -----------------------------
log_status("reports starting (serial)")

reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  verbose = TRUE
)

saveRDS(
  reports,
  file = file.path(out_dir, "reports_full_DNAm.rds")
)

log_status("reports complete; report output saved")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(pipeline_repo, "plotting_helpers.R"), chdir = TRUE)

fig_change <- file.path(out_dir, "Figures", "change")
fig_level <- file.path(out_dir, "Figures", "level")

dir.create(fig_change, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_level, recursive = TRUE, showWarnings = FALSE)

log_status("change-result plotting starting (serial)")
generate_all_plots(
  results,
  figures_dir = fig_change,
  analysis = "analysis_change"
)
log_status("change-result plotting complete")

log_status("level-result plotting starting (serial)")
generate_all_plots(
  results,
  figures_dir = fig_level,
  analysis = "analysis_level"
)
log_status("level-result plotting complete")
log_status(
  paste0(
    "run complete (",
    round(proc.time()[["elapsed"]] - run_started, 1),
    " sec)"
  )
)
