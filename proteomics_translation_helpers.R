
# Function: get_annotated_treatment_effects()
#
# Purpose:
#   Read treatment effect results from results_full_Proteomics.rds,
#   annotate them with protein naming information from
#   CALERIE_cleaned_protein_translation_table.csv, and return a data frame.
#
# Parameters:
#   analysis_type:
#     Which analysis result to use.
#     Allowed values: "change", "level"
#     Default: "change"
#
#   group:
#     Which subgroup to use.
#     Allowed values: "all", "male", "female"
#     Default: "all"
#
#   fu:
#     Follow-up timepoint filter.
#     Allowed values: 1, 2, 3, 4, or a vector such as c(1, 3)
#       1 = 3mo
#       2 = 6mo
#       3 = 12mo
#       4 = 24mo
#     Default: NULL, which returns all FU timepoints.
#
#   results_rds:
#     Path to the RDS file containing model results.
#     Default: "results_full_Proteomics.rds"
#
#   translation_csv:
#     Path to the protein translation table.
#     Default: "CALERIE_cleaned_protein_translation_table.csv"
#
# Output:
#   A data frame with treatment effects plus protein annotation columns.
#
# Examples:
#   df <- get_annotated_treatment_effects()
#
#   df_3mo <- get_annotated_treatment_effects(fu = 1)
#
#   df_level_female_12mo <- get_annotated_treatment_effects(
#     analysis_type = "level",
#     group = "female",
#     fu = 3
#   )

get_annotated_treatment_effects <- function(
    analysis_type = "change",
    group = "all",
    fu = NULL,
    results_rds = "results_full_Proteomics.rds",
    translation_csv = "CALERIE_cleaned_protein_translation_table.csv"
) {
  analysis_type <- match.arg(
    analysis_type,
    choices = c("change", "level")
  )
  
  group <- match.arg(
    group,
    choices = c("all", "male", "female")
  )
  
  if (!is.null(fu)) {
    fu <- as.integer(fu)
    
    if (any(is.na(fu)) || any(!fu %in% 1:4)) {
      stop("fu must be one or more of 1, 2, 3, 4. Leave fu unspecified to return all FU.")
    }
  }
  
  results <- readRDS(results_rds)
  
  translation <- read.csv(
    translation_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  analysis_name <- paste0("analysis_", analysis_type)
  
  treatment_effects <- results[[analysis_name]][[group]]$treatment_effects
  
  fu_labels <- c(
    "1" = "3mo",
    "2" = "6mo",
    "3" = "12mo",
    "4" = "24mo"
  )
  
  treatment_effects$FU_LABEL <- fu_labels[as.character(treatment_effects$FU)]
  
  if (!is.null(fu)) {
    treatment_effects <- treatment_effects[treatment_effects$FU %in% fu, ]
  }
  
  annotation_cols <- c(
    "AptName",
    "SomaId",
    "TargetFullName",
    "Target",
    "UniProt",
    "EntrezGeneID",
    "EntrezGeneSymbol",
    "Uniprot_Unique",
    "Symbol_Unique"
  )
  
  annotation <- translation[, annotation_cols]
  
  annotated_df <- merge(
    treatment_effects,
    annotation,
    by.x = "ANALYTE_NAME",
    by.y = "AptName",
    all.x = TRUE,
    sort = FALSE
  )
  
  output_cols <- c(
    "ANALYTE_NAME",
    "SomaId",
    "TargetFullName",
    "Target",
    "UniProt",
    "EntrezGeneID",
    "EntrezGeneSymbol",
    "Uniprot_Unique",
    "Symbol_Unique",
    "FU",
    "FU_LABEL",
    "EFFECT_SIZE",
    "SE",
    "P_VALUE",
    "BH_P_VALUE"
  )
  
  annotated_df <- annotated_df[, output_cols]
  
  annotated_df
}

