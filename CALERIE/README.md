# CALERIE analysis runners

Use the six runners below for new runs. The paths below are configured defaults; check each runner's paths for the target machine. Keep beta values, M-values, and clock scores as separate analyses.

## Canonical runners and outputs

Paths in the first column are relative to this directory. Listed result/report filenames are relative to the output directory.

| Analysis and runner | Output directory | Saved outputs |
| --- | --- | --- |
| Beta EWAS: [WTM_run_DNAm_betas.R](DNAm/EWAS_DNAm/WTM_run_DNAm_betas.R) | `~/FAST/Outputs/1.1.1` | `DNAm_betas_1.1.1.rds`, `reports_full_DNAm.rds`, `run_full_DNAm.log` |
| M-value EWAS: [run_DNAm_mvals.R](DNAm/EWAS_DNAm/run_DNAm_mvals.R) | `~/calerie_pipeline/calerie_DNAm` | `results_full_DNAm.rds`, `reports_full_DNAm.rds`, `run_full_DNAm.log` |
| Epigenetic clocks: [run_111_DNAmClocks.R](DNAm/Epigenetic_Clocks/run_111_DNAmClocks.R) | `~/FAST/Outputs/1.1.1/DNAmClock_1.1.1_Outputs` | `results.rds`, `reports.rds`, `run.log` |
| Protein PWAS: [combined_run.R](Proteomics/PWAS_Proteomics/combined_run.R) | `~/CALERIE/calerie_Proteomics` | `results_full_Proteomics.rds`, `reports_full_Proteomics.rds`, `annotated_treatment_effect_tables.rds`, `run_full_Proteomics.log` |
| Proteomic clocks: [run_Proteomic_Clocks.R](Proteomics/Proteomic_Clocks/run_Proteomic_Clocks.R) | `~/CALERIE/Proteomic_Clock` | `results_Proteomic_Clock.rds`, `reports_Proteomic_Clock.rds`, `run_Proteomic_Clock.log` |
| Metabolomics: [combined_run.R](Metabolomics/combined_run.R) | `~/CALERIE/Metabolomics` | `results_Metabolomics.rds`, `reports_Metabolomics.rds`, `run_Metabolomics.log` |

All six write plots under `Figures/change` and `Figures/level` in their output directory. Beta EWAS, M-value EWAS, and PWAS use `checkpoints` there. Use a fresh checkpoint directory after changing inputs, model settings, batch layout, or result schema; see [the checkpoint contract](../INPUTS_OUTPUTS.md).

## Inputs and follow-up coding

- **Beta EWAS:** `~/FAST/Data/CALERIE/Raw/DNAm/` contains `GRSet_fully_filtered_bmiq_chunk.rds`, `CALERIE_CPR_processed_pheno.rds`, `CALERIE_control_pcs_rgset_goodsamples.csv`, `CALERIE_SALAS_cell_counts_chunk.csv`, and `Cell_PCs.csv`.
- **M-value EWAS:** `~/CALERIE/DNAm/data/omics/GRSet_fully_filtered_bmiq_chunk_mvals.rds`; the phenotype file is under `data/pheno/`, and the same three covariate CSV filenames above are under `data/covariate/`.
- **Epigenetic clocks:** inputs under `~/FAST/Data/CALERIE/Raw/DNAm/`: `calerie_dunedinPACE_barcode.csv`, `CALERIE_PC-clocks_agedwb_CPR_fixed.csv`, `result_SystemAge`, `GrimAgeV2.csv`, and `CALERIE_FraminghamPACE_df_2.csv`. The script writes `Epigenetic_Clock_updated.csv` there, then reads it with the phenotype, control-PC, and cell-PC files listed above. Clock scores use `omics_type = "other"`.
- **PWAS:** `~/CALERIE/data/Proteomics/` contains `CALERIE_cleaned_log2_soma_matrix.csv`, `CALERIE_clinical_w_omics_crosswalks.csv`, and `CALERIE_GeneticPCs_20200727.xlsx`. Annotation additionally requires `CALERIE_cleaned_protein_translation_table.csv` in the PWAS output directory.
- **Proteomic clocks:** `~/CALERIE/Proteomic_Clock/data/` contains `calerie_organage.csv`, `CALERIE_tanaka_76_values.csv`, and `CALERIE_clinical_w_omics_crosswalks.csv`; genetic PCs use the PWAS path above. The script writes `data/omics_Proteomic_Clock.csv` with adjusted clock values and uses `omics_type = "other"`.
- **Metabolomics:** `~/CALERIE/Metabolomics/data/OBT_log.csv`, plus the clinical crosswalk and genetic-PC files at the PWAS paths above.

DNAm runners retain the supplied `fu` coding: baseline = `0`, 12 months = `1`, and 24 months = `2`. PWAS and proteomic clocks map months `0, 3, 6, 12, 24` to `FU = 0, 1, 2, 3, 4`. Metabolomics maps months `0, 6, 12, 24` to `FU = 0, 1, 2, 3`.

Beta EWAS and epigenetic clocks source `~/FAST/GitHub/track-1.1.1`; M-value EWAS sources `~/calerie_pipeline/repos/track-1.1.1-single-analyte`; PWAS, proteomic clocks, and metabolomics source `~/CALERIE/repos/track-1.1.1-single-analyte`. Check each script's path settings for the target machine.

## DunedinPACE covariance workflows

The beta runner invokes [dunedinpace_covariance.R](DNAm/EWAS_DNAm/dunedinpace_covariance.R) before its main EWAS. This computes the exact delete-one-subject jackknife of FAST probe effects, with pooled and available sex strata. Its production cache/output directory is:

```text
~/FAST/Outputs/1.1.1/DNAm_Bvals_1.1.1/dunedinpace_raw_beta_jackknife_v1
```

That directory contains pinned weights `sysdata.rda`, `config.rds`, leave-one-out batches under `checkpoints/<stratum>/`, and `dunedinpace_covariance.rds`. It is separate from the beta runner's primary EWAS results and `~/FAST/Outputs/1.1.1/checkpoints`.

[WTM_run_DunedinPACE_residualized_changes.R](DNAm/EWAS_DNAm/WTM_run_DunedinPACE_residualized_changes.R) is a separate pooled validation workflow. It uses the beta, phenotype, control-PC, and cell-PC inputs above and requires `sysdata.rda` from this cache. It writes `~/FAST/Outputs/1.1.1/dunedinpace_residualized_changes_v1.rds` and `.log`. Its separate-follow-up adjusted linear models produce residualized-change covariance and approximate treatment-effect covariance; they do not replace the mixed-model jackknife.

## Results and annotation

Keep raw nested model results for downstream analyses: `analysis_change` / `analysis_level` → `all` / `male` / `female` → `coefficients` / `treatment_effects`, plus metadata. Reports and annotated tables are separate artifacts.

PWAS writes `annotated_treatment_effect_tables.rds` as a list of 24 flat tables: two analyses × three strata × four follow-ups, with names such as `change_all_3mo`. These tables do not replace `results_full_Proteomics.rds`.

[run_translation.R](Proteomics/PWAS_Proteomics/run_translation.R) reruns annotation without model fitting. It still reads/writes under `~/calerie_pipeline/calerie_Proteomics` and sources the `track-1.1.1-single-analyte-branch` checkout there; this differs from the canonical PWAS paths. Point it at the intended existing results, translation CSV, and helper checkout before use. Its output filename remains `annotated_treatment_effect_tables.rds`.

## Historical entry points

Retained for provenance; use their canonical counterparts above for new runs:

- [PWAS run_Proteomics.R](Proteomics/PWAS_Proteomics/run_Proteomics.R): prepared-workspace wrapper; use PWAS `combined_run.R`.
- [Proteomic-clock combined_run.R](Proteomics/Proteomic_Clocks/combined_run.R): older duplicate using `omics_type = "Proteomics"`; use `run_Proteomic_Clocks.R`.
- [run_Proteomics_Clock.R](Proteomics/Proteomic_Clocks/run_Proteomics_Clock.R): prepared-workspace wrapper; use `run_Proteomic_Clocks.R`.
- [Metabolomics run_Metabolomics.R](Metabolomics/run_Metabolomics.R): prepared-workspace wrapper; use metabolomics `combined_run.R`.
- [OLD_run_Epigenetic_Clocks.R](DNAm/Epigenetic_Clocks/OLD_run_Epigenetic_Clocks.R): archived clock runner; use `run_111_DNAmClocks.R` with FraminghamPACE inputs.
