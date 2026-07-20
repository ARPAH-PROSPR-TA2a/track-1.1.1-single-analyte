old_wd <- setwd("..")
on.exit(setwd(old_wd), add = TRUE)

source("main.R")
source("plotting_helpers.R")

pheno <- readRDS("Examples/ExampleData/pheno_example.rds")
omics <- readRDS("Examples/ExampleData/proteomics_log2.rds")
dir.create("Examples/ExampleResults", showWarnings = FALSE)
figures_dir <- "Examples/ExampleFigures/Proteomics"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

additional_covariates <- c("agebl", "agevis", "mbmi")

results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)

reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

saveRDS(results, "Examples/ExampleResults/proteomics_results.rds")
saveRDS(reports, "Examples/ExampleResults/proteomics_reports.rds")

generate_all_plots(results$analysis_change, figures_dir = file.path(figures_dir, "analysis_change"))
generate_all_plots(results$analysis_level, figures_dir = file.path(figures_dir, "analysis_level"))

cat("Proteomics example complete\n")
cat("  Pheno samples: ", nrow(pheno), "\n", sep = "")
cat("  Proteomic analytes: ", nrow(omics), "\n", sep = "")
cat("  Change coefficients: ", nrow(results$analysis_change$all$coefficients), "\n", sep = "")
cat("  Level coefficients: ", nrow(results$analysis_level$all$coefficients), "\n", sep = "")
cat("  Report summary rows: ", nrow(reports$pheno_summary), "\n", sep = "")
cat("  Figures: ", figures_dir, "\n", sep = "")
