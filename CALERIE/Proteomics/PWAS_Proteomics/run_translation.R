setwd("~/calerie_pipeline/calerie_Proteomics")
source("~/calerie_pipeline/repos/track-1.1.1-single-analyte-branch/proteomics_translation_helpers.R")

results_rds <- "results_full_Proteomics.rds"
translation_csv <- "CALERIE_cleaned_protein_translation_table.csv"
output_rds <- "annotated_treatment_effect_tables.rds"

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

saveRDS(annotated_tables, output_rds)

cat("Saved", length(annotated_tables), "annotated treatment-effect tables to:", output_rds, "\n")
cat("Table names:\n")
print(names(annotated_tables))
