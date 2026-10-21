library(readr)
library(dplyr)
library(forcats)
library(tibble)
library(future)
library(furrr)
library(emmeans)
library(lme4)
library(lmerTest)

# -----------------------------
# Create clock data
# -----------------------------

DunedinPace <- read_csv(file = "~/FAST/Data/CALERIE/Raw/DNAm/calerie_dunedinPACE_barcode.csv")
PCClocks <- read_csv(file = "~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_PC-clocks_agedwb_CPR_fixed.csv")
SystemAge <- read_rds("~/FAST/Data/CALERIE/Raw/DNAm/result_SystemAge")
GrimAgeV2 <- read_csv("~/FAST/Data/CALERIE/Raw/DNAm/GrimAgeV2.csv")
FraminghamPACE <- read_csv("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_FraminghamPACE_df_2.csv")

dup_check <- bind_rows(
  DunedinPace |> count(barcode) |> filter(n > 1) |> mutate(dataset = "DunedinPace"),
  PCClocks     |> count(barcode) |> filter(n > 1) |> mutate(dataset = "PCClocks"),
  SystemAge   |> count(barcode) |> filter(n > 1) |> mutate(dataset = "SystemAge"),
  GrimAgeV2   |> count(barcode) |> filter(n > 1) |> mutate(dataset = "GrimAgeV2")
)

dup_check

# Keep first row per barcode
DunedinPace_clean <- DunedinPace |>
  distinct(barcode, .keep_all = TRUE)

PCClocks_clean <- PCClocks |>
  select(barcode, 
         PCGrimAge,
         PCPhenoAge,
         PCHorvath1,
         PCHorvath2,
         PCHannum) |>
  distinct(barcode, .keep_all = TRUE)

GrimAgeV2_clean <- GrimAgeV2 |>
  select(barcode, GrimAgeV2) |>
  distinct(barcode, .keep_all = TRUE)

SystemAge_clean <- SystemAge |>
  select(barcode, SystemsAge) |>
  distinct(barcode, .keep_all = TRUE)

FraminghamPACE_clean <- FraminghamPACE |>
  rename(barcode=BARCODE) |>
  select(barcode, ElasticNet, Ridge) |>
  distinct(barcode, .keep_all = TRUE)


# Confirm barcode is unique
stopifnot(!anyDuplicated(DunedinPace_clean$barcode))
stopifnot(!anyDuplicated(PCClocks_clean$barcode))
stopifnot(!anyDuplicated(GrimAgeV2_clean$barcode))
stopifnot(!anyDuplicated(SystemAge_clean$barcode))

# Join all clocks
Epigenetic_Clock <- DunedinPace_clean |>
  full_join(
    PCClocks_clean,
    by = "barcode",
    relationship = "one-to-one"
  ) |>
  full_join(
    GrimAgeV2_clean,
    by = "barcode",
    relationship = "one-to-one"
  ) |>
  full_join(
    SystemAge_clean,
    by = "barcode",
    relationship = "one-to-one"
  ) |>
  full_join(
    FraminghamPACE_clean,
    by = "barcode",
    relationship = "one-to-one"
  )

write_csv(Epigenetic_Clock, file = "~/FAST/Data/CALERIE/Raw/DNAm/Epigenetic_Clock_updated.csv")

# -----------------------------
# Paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/Epigenetic_Clock_updated.csv")
pheno_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds")
control_pc_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_control_pcs_rgset_goodsamples.csv")
cell_pcs_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/Cell_PCs.csv")

pipeline_repo <- path.expand("~/FAST/GitHub/track-1.1.1")
out_dir <- path.expand("~/FAST/Outputs/1.1.1//DNAmClock_1.1.1_Outputs")

omics_type <- "other"
n_cores <- 3
results_path <- file.path(out_dir, "results.rds")
reports_path <- file.path(out_dir, "reports.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")

cat("START run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Read inputs
# -----------------------------
omics_raw <- read_csv(omics_raw_path, show_col_types = FALSE)
pheno_raw <- readRDS(pheno_raw_path)
control_pc <- read.csv(control_pc_path, stringsAsFactors = FALSE)
cell_PCs <- read_csv(cell_pcs_path, show_col_types = FALSE)

# -----------------------------
# Input QC
# -----------------------------
stopifnot("barcode" %in% colnames(omics_raw))
stopifnot(!any(duplicated(omics_raw$barcode)))
stopifnot(!any(duplicated(colnames(omics_raw))))

stopifnot(!any(duplicated(pheno_raw$Barcode)))
stopifnot(!any(duplicated(control_pc$filenames)))
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
    cell_PCs |> select(SAMPLE_ID, starts_with("cell_PC")),
    by = c("Barcode" = "SAMPLE_ID")
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
    female = haven::as_factor(female),
    CR = haven::as_factor(CR),
    deidsite = haven::as_factor(deidsite),
    bmistrat = haven::as_factor(bmistrat)
  ) |>
  rename(
    SAMPLE_ID = Barcode,
    SUBJECT_ID = Participant_ID,
    FU = fu,
    FEMALE = female,
    TREATMENT_GROUP = CR
  ) |>
  select(-Time_Point)

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
  filter(SAMPLE_ID %in% omics_raw$barcode) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x)))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% omics_raw$barcode))

# -----------------------------
# Build omics table
# -----------------------------
omics <- omics_raw |>
  column_to_rownames("barcode") |>
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
  n_cores = n_cores
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