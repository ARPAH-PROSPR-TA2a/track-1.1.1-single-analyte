setwd("~/calerie_pipeline/repos/track-1.1.1-single-analyte-branch")
source("main.R")

load("~/calerie_pipeline/calerie_Proteomics_Clock/input_Proteomics_Clock.Rdata")

stopifnot(
  exists("FAST_omics_WAS"),
  exists("FAST_omics_WAS_reports"),
  exists("pheno"),
  exists("omics"),
  exists("covariates")
)



out_dir <- path.expand("~/calerie_pipeline/calerie_Proteomics_Clock")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_Proteomics_Clock.log")

cat("START full Proteomics Clock run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# analysis
# -----------------------------
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = covariates,
)

saveRDS(
  results,
  file = file.path(out_dir, "results_Proteomics_Clock.rds")
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
  file = file.path(out_dir, "reports_Proteomics_Clock.rds")
)

cat("DONE reports: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")


## plot
setwd("~/calerie_pipeline/repos/track-1.1.1-single-analyte-branch")
source("plotting_helpers.R")


fig_change <- "~/calerie_pipeline/calerie_Proteomics_Clock/Figures/change"
fig_level  <- "~/calerie_pipeline/calerie_Proteomics_Clock/Figures/level"

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
