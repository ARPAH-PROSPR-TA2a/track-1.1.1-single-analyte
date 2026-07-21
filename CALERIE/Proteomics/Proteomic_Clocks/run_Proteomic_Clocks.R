library(dplyr)
library(forcats)
library(readxl)
library(tidyr)

# -----------------------------
# Paths and parameters
# -----------------------------
out_dir <- path.expand("~/CALERIE/Proteomic_Clock")
data_dir <- file.path(out_dir, "data")

organ_path <- file.path(data_dir, "calerie_organage.csv")
tanaka_path <- file.path(data_dir, "CALERIE_tanaka_76_values.csv")
pheno_raw_path <- file.path(data_dir, "CALERIE_clinical_w_omics_crosswalks.csv")
genetic_pcs_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_GeneticPCs_20200727.xlsx")

pipeline_repo <- path.expand("~/CALERIE/repos/track-1.1.1-single-analyte")

omics_type <- "Proteomics"
clock_omics_csv <- file.path(data_dir, "omics_Proteomic_Clock.csv")
results_path <- file.path(out_dir, "results_Proteomic_Clock.rds")
reports_path <- file.path(out_dir, "reports_Proteomic_Clock.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_Proteomic_Clock.log")

cat("START run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Read clock inputs
# -----------------------------
organ_raw <- read.csv(
  organ_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

tanaka_raw <- read.csv(
  tanaka_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

pheno_clock_raw <- read.csv(
  pheno_raw_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

# -----------------------------
# Clock input QC
# -----------------------------
stopifnot(all(c("SlideId_SubArray", "Organ", "Predicted_Age") %in% colnames(organ_raw)))
stopifnot(all(c("row_names", "tanaka_76") %in% colnames(tanaka_raw)))
stopifnot(all(c("DEID", "fumo", "agebl", "female", "proteomics_barcode") %in% colnames(pheno_clock_raw)))

organ <- data.frame(
  proteomics_barcode = organ_raw$SlideId_SubArray,
  Organ = organ_raw$Organ,
  Predicted_Age = organ_raw$Predicted_Age,
  stringsAsFactors = FALSE
)

tanaka <- data.frame(
  proteomics_barcode = tanaka_raw$row_names,
  Organ = "Tanaka_76",
  Predicted_Age = tanaka_raw$tanaka_76,
  stringsAsFactors = FALSE
)

omics_long <- rbind(organ, tanaka)

stopifnot(!any(is.na(omics_long)))
stopifnot(!any(duplicated(omics_long[, c("proteomics_barcode", "Organ")])))

# -----------------------------
# Create adjusted clock omics
# -----------------------------
pheno_clock <- pheno_clock_raw |>
  filter(!is.na(proteomics_barcode)) |>
  select(DEID, fu, fumo, agebl, female, proteomics_barcode) |>
  mutate(
    fu = case_when(
      fumo == 0 ~ 0,
      fumo == 3 ~ 1,
      fumo == 6 ~ 2,
      fumo == 12 ~ 3,
      fumo == 24 ~ 4
    ),
    fu = factor(fu, levels = c(0, 1, 2, 3, 4))
  ) |>
  distinct(proteomics_barcode, .keep_all = TRUE)

clock_long <- inner_join(
  pheno_clock,
  omics_long,
  by = "proteomics_barcode"
)

stopifnot(nrow(clock_long) > 0)
stopifnot(!any(is.na(clock_long)))

clock_long$Baseline_Predicted_Age <- NA_real_

for (organ_name in sort(unique(clock_long$Organ))) {
  fit_rows <- clock_long$Organ == organ_name &
    clock_long$fu == 0 &
    !is.na(clock_long$Predicted_Age) &
    !is.na(clock_long$agebl)
  
  predict_rows <- clock_long$Organ == organ_name
  
  stopifnot(sum(fit_rows) >= 2)
  
  fit <- lm(Predicted_Age ~ agebl, data = clock_long[fit_rows, ])
  clock_long$Baseline_Predicted_Age[predict_rows] <-
    predict(fit, newdata = clock_long[predict_rows, ])
}

clock_long$Adjusted_Predicted_Age <-
  clock_long$Predicted_Age - clock_long$Baseline_Predicted_Age

clock_long <- clock_long[order(clock_long$Organ, clock_long$DEID, clock_long$fu), ]

write.csv(
  clock_long,
  clock_omics_csv,
  row.names = FALSE,
  na = ""
)

# -----------------------------
# Read phenotype covariates
# -----------------------------
pheno_raw <- read.csv(
  pheno_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

PCs <- read_excel(genetic_pcs_path)

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
  filter(SAMPLE_ID %in% clock_long$proteomics_barcode) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x)))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% clock_long$proteomics_barcode))

# -----------------------------
# Build omics table
# -----------------------------
omics <- clock_long |>
  select(Organ, proteomics_barcode, Adjusted_Predicted_Age) |>
  pivot_wider(
    names_from = proteomics_barcode,
    values_from = Adjusted_Predicted_Age
  ) |>
  rename(ANALYTE_NAME = Organ)

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

