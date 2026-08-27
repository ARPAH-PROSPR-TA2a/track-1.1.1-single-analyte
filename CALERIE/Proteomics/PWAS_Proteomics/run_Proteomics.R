setwd("~/calerie_pipeline/repos/track-1.1.1-single-analyte-branch")
source("main.R")

load("~/calerie_pipeline/calerie_Proteomics/input_Proteomic.Rdata")

stopifnot(
  exists("FAST_omics_WAS"),
  exists("FAST_omics_WAS_reports"),
  exists("pheno"),
  exists("omics"),
  exists("covariates")
)



out_dir <- path.expand("~/calerie_pipeline/calerie_Proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_full_Proteomics.log")

cat("START full Proteomics run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# analysis
# -----------------------------
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = covariates,
  n_cores = 30,
  checkpoint_dir = file.path(out_dir, "checkpoints"),
  checkpoint_batch_size = 2000L
)

saveRDS(
  results,
  file = file.path(out_dir, "results_full_Proteomics.rds")
)

cat("DONE analysis: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# reports
# -----------------------------
reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = covariates
)

saveRDS(
  reports,
  file = file.path(out_dir, "reports_full_Proteomics.rds")
)

cat("DONE reports: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")


## plot
setwd("~/calerie_pipeline/repos/track-1.1.1-single-analyte")
source("plotting_helpers.R")


fig_change <- "~/calerie_pipeline/calerie_Proteomics/Figures/change"
fig_level  <- "~/calerie_pipeline/calerie_Proteomics/Figures/level"

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
