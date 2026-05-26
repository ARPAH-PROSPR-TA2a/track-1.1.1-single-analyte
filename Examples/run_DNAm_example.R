old_wd <- setwd("..")
on.exit(setwd(old_wd), add = TRUE)

source("main.R")

pheno <- readRDS("Examples/ExampleData/pheno_example.rds")
omics <- readRDS("Examples/ExampleData/dnam_mvalues.rds")
dir.create("Examples/ExampleResults", showWarnings = FALSE)

additional_covariates <- c("agebl", "agevis", "mbmi")

results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  n_cores = 1
)

reports <- FAST_omics_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

saveRDS(results, "Examples/ExampleResults/dnam_results.rds")
saveRDS(reports, "Examples/ExampleResults/dnam_reports.rds")

cat("DNAm example complete\n")
cat("  Pheno samples: ", nrow(pheno), "\n", sep = "")
cat("  DNAm analytes: ", nrow(omics), "\n", sep = "")
cat("  Change coefficients: ", nrow(results$analysis_change$all$coefficients), "\n", sep = "")
cat("  Level coefficients: ", nrow(results$analysis_level$all$coefficients), "\n", sep = "")
cat("  Report summary rows: ", nrow(reports$pheno_summary), "\n", sep = "")
