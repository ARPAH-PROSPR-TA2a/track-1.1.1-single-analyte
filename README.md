# Omics-Wide Association Study Pipeline

This repository, written in R, provides a reliable implementation of
omics-wide association studies for Proteomics, Metabolomics, and DNA
Methylation. It standardizes inference on `lm()` (single follow-up) and
`lmerTest::lmer()` (multiple follow-ups) and automatically stratifies
results by gender when data permits.

The pipeline exposes two functions: `FAST_omics_WAS()` for the statistical analysis (parallelized, resumable) and `FAST_omics_WAS_reports()` for QC and data summary reports.

## Installation

Clone the repository:

``` bash
git clone https://github.com/ARPAH-PROSPR-TA2a/track-1.1.1-single-analyte.git
```

In your R script or console, source the main module from your local checkout:

``` r
source("/path/to/track-1.1.1-single-analyte/main.R")
```

The pipeline automatically sources all necessary helper functions from
the same directory, without changing your working directory.

## FAST_omics_WAS()

The main function for conducting omics-wide association analyses.

### Quick Example

``` r
# Load your phenotype and omics data
pheno <- read.csv("my_phenotype.csv")
omics <- read.csv("my_omics_data.csv")

# Run the association analysis (parallelized across 8 cores, with checkpointing)
results <- FAST_omics_WAS(
  pheno  = pheno,
  omics  = omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi"),
  n_cores = 8,
  checkpoint_dir = "checkpoints",
  verbose = TRUE
)

# View analysis results
head(results$analysis_change$all$coefficients)
head(results$analysis_change$all$treatment_effects)

head(results$analysis_level$all$coefficients)
head(results$analysis_level$all$treatment_effects)

# Generate QC reports separately
reports <- FAST_omics_WAS_reports(
  pheno  = pheno,
  omics  = omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi"),
  verbose = TRUE
)

reports$pheno_summary
reports$variable_summaries$all$omics_FU0_Tx0
```

For the usual binary input, `TREATMENT_GROUP` is `0` for control and `1`
for treatment. Labels default to `"Control"` and `"Treatment"`; optionally pass
`treatment_labels = c("0" = "Usual care", "1" = "Diet")` to either function.

### Multiple intervention arms

Use distinct codes directly in `TREATMENT_GROUP`: `0` for the shared control,
`1` for one active arm, `2` for another, and so on. Supply a named label for
**every** observed code, including control:

```r
# pheno$TREATMENT_GROUP contains 0, 1, or 2, constant within each subject.
labels <- c(`0` = "Usual care", `1` = "Diet", `2` = "Diet + exercise")
results <- FAST_omics_WAS(pheno, omics, treatment_labels = labels, n_cores = 1)
reports <- FAST_omics_WAS_reports(pheno, omics, treatment_labels = labels)

# Separate treatment-effect data frames within this stratum.
te <- results$analysis_change$all$treatment_effects
effects <- split(te, te$TREATMENT)
effects[["Diet"]]
effects[["Diet + exercise"]]
```

Codes are categories and are never renumbered. Each analyte is fitted once
per analysis and stratum with all arms together, then each active arm is
compared with the shared control at each follow-up. Unequal group sizes are
supported; multiple follow-ups use arm-by-follow-up interactions. The original
input data frames are unchanged.

For a fully synthetic example, open
[`Examples/run_multiarm_example.R`](Examples/run_multiarm_example.R) in RStudio
and click **Source**, or run its sections in order starting with setup. The
script finds the pipeline from its own location; your working directory can
stay where it is. `pheno`, `omics`, `results`, `reports`, `effects_diet`, and
`effects_diet_exercise` remain in your R workspace for inspection.
Running individual sections uses `rstudioapi` to locate the open script;
if it is unavailable, use **Source**.

From an R console, you can also run:

```r
source("/path/to/Track1.1.1/Examples/run_multiarm_example.R")
```

### Parameters

-   **`pheno`** (data.frame): Phenotype data with required columns (see
    Data Format section below)
-   **`omics`** (data.frame): Omics data with required structure (see
    Data Format section below)
-   **`omics_type`** (character): Type of omics data being analyzed.
    Options: `"Proteomics"`, `"DNAm"`, `"Metabolomics"`, `"other"`.
    Use `"other"` for non-omics feature data; it produces no preprocessing
    reminder and applies no type-specific filtering. Default: `"Proteomics"`
-   **`additional_covariates`** (character vector, optional): Names of
    additional covariates to include in regression models. Must be
    column names in the `pheno` data frame. These columns must be
    numeric, factor, or logical.
-   **`treatment_labels`** (named character vector): Labels keyed by every
    observed `TREATMENT_GROUP` code, including `0`. Required with multiple
    active arms; optional for binary input, which defaults to `"Control"` and
    `"Treatment"`. Labels must be unique and nonblank. See
    [Multiple intervention arms](#multiple-intervention-arms).
-   **`n_cores`** (integer, optional): Number of cores to use for
    parallelization. Defaults to `NULL`, which auto-detects as
    `max(1, detectCores() - 1)` and falls back to `1` if detection is unavailable.
    Set to `1` to run serially. Note:
    `detectCores()` may overcount in HPC/container environments — set
    explicitly if running on a cluster with allocated core limits.
-   **`checkpoint_dir`** (character, optional): Path to a directory for
    saving per-batch checkpoints. If `NULL` (default), checkpointing is
    disabled. If provided, completed batches are saved to disk and
    automatically reused on resume, so a crashed run can pick up where
    it left off. The directory is created if it does not exist. See
    [Checkpointing](#checkpointing) below.
-   **`checkpoint_batch_size`** (integer): Number of analytes per
    checkpoint batch. Default: `2000`. Only relevant when
    `checkpoint_dir` is set.
-   **`verbose`** (logical): Emit `[1.1.1]` progress messages for serial
    input quality checks (every 100,000 analytes),
    preparation, the configured parallel backend and effective worker count,
    cached versus pending batches, computed-batch timing and observed worker
    PIDs, multiple-testing correction, and report stages. Default: `FALSE`.

### Return Value

A list containing:

-   **`$analysis_change`**: Change-score analysis (follow-up minus baseline, baseline-adjusted)
-   **`$analysis_level`**: Level analysis (absolute follow-up, baseline-adjusted)
-   **`$metadata`**: Result schema and treatment mapping. In
    `$metadata$treatments`, `ARM` is the input arm code, `TREATMENT` is
    its display label, `REFERENCE` is the control label, and `CONTRAST` identifies
    the comparison (`"1_vs_0"`, `"2_vs_0"`, ...).

Both analyses adjust for the same baseline measurement. On the same data and
model specification, their treatment contrasts are mathematically identical;
the baseline coefficient changes between the response parameterizations.

Each analysis contains results stratified by sex. Male and female results are only generated if both sexes are present in the data:

-   **`$all`**: Results from the full dataset
-   **`$male`**: Results from male subset
-   **`$female`**: Results from female subset

Each stratum contains:

-   **`$coefficients`**: Data frame of all model coefficients for each
    analyte (ANALYTE_NAME, COEFFICIENT, N_OBS, EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE)
-   **`$treatment_effects`**: Data frame of treatment effects at each
    follow-up level and active arm (ANALYTE_NAME, FU, CONTRAST, TREATMENT,
    REFERENCE, EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE, STATUS).

Use `TREATMENT` to select a named intervention. Original codes are retained in
the contrast identifiers; `$metadata$treatment_labels` records the full label map.
Rows for contrasts that cannot be estimated are retained with missing
estimates and `STATUS = "not_estimable"`; failed fits use `"fit_error"` and
successful contrasts use `"ok"`. BH correction is performed across
analytes separately for each contrast and follow-up, within each analysis and
stratum; the underlying contrast p-values have no cross-arm adjustment.

**DNAm only**: `BH_P_VALUE_FILTERED` is added to `coefficients` and
`treatment_effects` for a pre-specified filtered probe set (useful for
assessing significance under different multiple-testing burdens).
DNAm `treatment_effects` also includes logical `IN_FILTERED_SET`, equal to
`ANALYTE_NAME %in% filtered_probes`. This records membership even when inference
is missing: `BH_P_VALUE_FILTERED = NA` alone cannot distinguish a probe outside
the set from a member with no p-value. Membership is the same for a probe across
follow-ups, contrasts, analyses, and strata, regardless of `STATUS`. It is added
after batch assembly and is absent from coefficient tables and non-DNAm results.

---

## FAST_omics_WAS_reports()

Generates QC and data summary reports. Takes the same `pheno`, `omics`,
`omics_type`, `additional_covariates`, and `treatment_labels` arguments as `FAST_omics_WAS()`,
plus `verbose`, and runs the same input validation, but only produces reports
— no model fitting. With `verbose = TRUE`, long serial per-analyte report
steps emit coarse progress milestones.

Separating reports from analysis allows long-running analyses to be
parallelized and checkpointed without re-running reporting, and vice versa.

### Quick Example

``` r
reports <- FAST_omics_WAS_reports(
  pheno  = pheno,
  omics  = omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi")
)
```

### Return Value

A list containing:

**`$pheno_summary`** (study-level): Data frame with one row per (FU, FEMALE) cell
giving N_SUBJECTS, N_CONTROL, N_TREATMENT (subject-level) and N_SAMPLES (row-level).
`N_TREATMENT` counts all active arms together.

**`$arm_counts`** (study-level): One row per FU × FEMALE × arm, including
zero-count cells. Columns: FU, FEMALE, ARM, TREATMENT, IS_CONTROL, N_SUBJECTS,
N_SAMPLES. Use this table for separate intervention counts.

**`$metadata`**: Schema version `"2.0"`, the full `treatment_labels` map, and
the treatment comparisons used by the reports.

**`$variable_summaries`**: Sex-stratified (`$all`, `$male`, `$female`) omics and
covariate summaries. Within each stratum, one entry per FU × treatment group cell
(e.g. `$omics_FU0_Tx0`, `$omics_FU1_Tx1`, `$covariates_FU0_Tx0`, ...). Cells with
no samples are omitted. Each omics entry is a per-analyte data frame (N_NONMISSING,
MEAN, MEDIAN, SD, MIN, MAX); each covariates entry mirrors the same structure as
before (COVARIATE_NAME, TYPE, N_NA, SUMMARY).

**`$randomization_reports`** (study-level, not sex-stratified):

-   **`$omics_report`**: Per-analyte baseline balance check (Welch's t-test,
    treatment vs control; MEAN_DIFFERENCE, COHENS_D, SE, P_VALUE, BH_P_VALUE)
-   **`$covariate_report`**: Baseline balance check for `FEMALE` (if both sexes
    present) and any additional covariates. Numeric: Welch's t-test; logical:
    chi-squared; factor: chi-squared if min cell count ≥ 5, Fisher's exact
    otherwise (simulated for tables larger than 2×2). Columns: VARIABLE, TYPE,
    TEST, STATISTIC, MIN_CELL_COUNT, P_VALUE, SUMMARY_CONTROL, SUMMARY_TREATMENT.

For multiple active arms, each baseline balance table contains separate
active-arm-versus-control comparisons identified by `CONTRAST`, `TREATMENT`,
and `REFERENCE`. Omics balance BH correction is applied within each contrast.

---

## Checkpointing

For large analyses (e.g. 750k+ DNAm probes), runs can take many hours.
Passing `checkpoint_dir` enables fault-tolerant execution: analytes are
processed in batches, and each completed batch is saved atomically to disk.
If the run is interrupted, re-running with the same `checkpoint_dir` skips
already-completed batches and continues from where it left off.

``` r
# First run — saves progress to checkpoints/
results <- FAST_omics_WAS(
  pheno, omics,
  omics_type = "DNAm",
  n_cores = 16,
  checkpoint_dir = "checkpoints",
  checkpoint_batch_size = 2000
)

# If interrupted, re-run the same call — completed batches are loaded from disk
results <- FAST_omics_WAS(
  pheno, omics,
  omics_type = "DNAm",
  n_cores = 16,
  checkpoint_dir = "checkpoints",
  checkpoint_batch_size = 2000
)
```

Checkpoints are organized as:
```
checkpoints/
  change/all/    batch_1.rds, batch_2.rds, ...
  change/male/   batch_1.rds, ...
  change/female/ batch_1.rds, ...
  level/all/     batch_1.rds, ...
  ...
```

Resume with the same input data, model settings, and batch layout. Use a fresh
checkpoint directory when any of these change or when upgrading across this
result-schema change. You can change `n_cores` or display labels
in `treatment_labels` on resume; labels are applied after the batches are
assembled and do not require refitting.

## Data Format Requirements

### Phenotype Data

Phenotype data must be a data frame with the following **required
columns**:

| Column                  | Type                     | Description                                          | Valid Values                            |
|----------------|----------------|------------------------|------------------|
| `SAMPLE_ID`             | Character                | Unique identifier for each sample                    | Unique, no duplicates                   |
| `SUBJECT_ID`            | Character                | Subject identifier for tracking across follow-ups    | All SUBJECT_ID\*FU pairs must be unique |
| `FU`                    | Factor                   | Follow-up timepoint                                  | Levels: 0 (baseline), 1..max (follow-ups), consecutive integers |
| `FEMALE`                | Factor                   | Sex indicator                                        | Levels: 0 (male), 1 (female)            |
| `TREATMENT_GROUP`        | Factor or numeric        | Treatment assignment                                 | `0` = control, `1` = treatment          |
| *Additional Covariates* | *Factor/Numeric/Logical* | *Optional columns for additional analysis variables* | *Any valid values*                      |

If you wish to add additional intervention arms, use codes `2`, `3`, and so on
in `TREATMENT_GROUP` and supply
[`treatment_labels`](INPUTS_OUTPUTS.md#treatment_labels) for every group,
including control.

**Data requirements:**

-   Every subject must have both a baseline sample (`FU == 0`) and at
    least one follow-up sample (`FU >= 1`); subjects missing either are
    dropped during validation
-   Must contain both males (`FEMALE == 0`) and females (`FEMALE == 1`),
    otherwise gender stratification will be skipped
-   Must contain control (`TREATMENT_GROUP == 0`) and at least one active
    arm (`TREATMENT_GROUP > 0`)
-   No duplicate `SAMPLE_ID` values
-   All `SUBJECT_ID*FU` pairs must be unique
-   `FU`, `FEMALE`, and `TREATMENT_GROUP` must be factors (or numeric
    values that will be automatically converted to factors)
-   `FEMALE` and `TREATMENT_GROUP` must be non-missing for all samples
-   Treatment assignment must remain constant across visits for each subject

**Example structure:**

```         
SAMPLE_ID      SUBJECT_ID  FU  TREATMENT_GROUP  FEMALE  agebl  agevis
sample_001     subj_001    0   1               1       55     56
sample_002     subj_001    1   1               1       55     57
sample_003     subj_002    0   0               0       62     63
sample_004     subj_002    1   0               0       62     64
```

**Additional covariates (optional):**

-   Any columns not listed above can be included as additional
    covariates
-   Must be numeric, factor, or logical
-   Samples with NA values in any additional covariate are dropped
    during validation; subjects left without both a baseline and a
    follow-up after this filtering are then dropped as well
-   Specified using the `additional_covariates` parameter, a character
    vector that names the columns within `pheno` that should be included
    as additional covariates in the models

### Omics Data

Omics data must be a data frame with:

| Column            | Type      | Description                                                                 |
|-------------|-------------|---------------------------------------------|
| `ANALYTE_NAME`    | Character | Unique identifier for each feature (analyte, protein, metabolite, CpG site) |
| Sample ID columns | Numeric   | One column per sample, named exactly as in `pheno$SAMPLE_ID`                |

**Data requirements:**

-   `ANALYTE_NAME` column is required
-   Column names (except `ANALYTE_NAME`) must match `SAMPLE_ID` values
    from phenotype data, the rest will be dropped. The final number of
    samples will be reported
-   All measurement columns must be numeric
-   Analytes with NA values or near-zero variance will generate a
    warning but are included in analysis
-   Any omics samples not in phenotype data are automatically excluded
-   Any phenotype samples not in omics data are automatically excluded

**Example structure:**

```         
ANALYTE_NAME    sample_001   sample_002   sample_003   sample_004
cg00000029      0.602        0.515        0.684        0.598
cg00000103      0.456        0.468        0.412        0.401
cg00000109      0.721        0.735        0.691        0.702
```

## Analysis Methods

Within each of the two analysis types (`analysis_change` and `analysis_level`), the
pipeline selects the appropriate model based on follow-up structure:

```         
Check maximum follow-up level
├─ max FU == 1 → Use LM (linear regression via lm)
└─ max FU > 1 → Use LME4 (linear mixed effects via lmerTest::lmer)
```

## Plotting Results

The pipeline includes plotting functions for visualizing treatment effects.
Source `plotting_helpers.R` and call `generate_all_plots()` on either the
`change` or `level` result:

``` r
source("plotting_helpers.R")

results <- FAST_omics_WAS(pheno, omics, omics_type = "DNAm",
                          additional_covariates = c("age", "bmi"),
                          n_cores = 8)

generate_all_plots(results$analysis_change)                          # saves to Figures/
generate_all_plots(results$analysis_change, figures_dir = "my_dir")  # saves to my_dir/

generate_all_plots(results$analysis_level)                           # optional
```

This auto-detects all dimensions of your results (strata, FU levels,
contrasts, filtered probes) and outputs multi-page PDFs with labeled contrast
pages. For DNAm results, both
full and filtered probe set plots are generated.
DNAm plots use `IN_FILTERED_SET` for filtered membership and show labeled empty
plots when a cell has no available inference.

For a direct `plot_qq()` or `plot_volcano()` call with multiple treatments,
supply your label, for example `treatment = "Diet"`. With one treatment it is
selected automatically. `generate_all_plots()` includes all treatments.

Two plot types are produced:

-   **QQ plots**: Observed vs expected p-value distributions (via `qqman`)
-   **Volcano plots**: Effect size vs significance, with Up/Down/Not
    significant coloring based on BH-corrected p-values (via `ggplot2`)

Requires packages: `qqman`, `ggplot2`.

## Testing

Run these from the repository root. They generate synthetic measurements;
the DNAm checks use public probe names from the bundled lists. No study data
or `PracticeData` directory is needed.

```bash
Rscript --vanilla test_multiarm.R
Rscript --vanilla test_multiarm_reports.R
Rscript --vanilla test_dnam_filtered_membership.R
Rscript --vanilla test_comprehensive.R
Rscript --vanilla test_parallel_checkpoint.R
Rscript --vanilla test_resource_paths_and_validation_progress.R
Rscript --vanilla test_example_execution.R
```

The checks cover independent joint-model contrasts, unchanged binary estimates,
missing arms, input validation, BH families, labeled reports/plots/annotations,
explicit DNAm filtered-set membership despite unavailable inference,
fresh-session checkpoint resume, and example execution outside the repository
working directory. The parallel suite also checks a
5,000-analyte input under a restricted future-globals memory limit.

## Additional Documentation

-   [Inputs and outputs](INPUTS_OUTPUTS.md): required phenotype and omics
    input formats, result object structure, and report outputs
-   [Code walkthrough](CODE_WALKTHROUGH_v3.md): implementation overview
    for validation, analysis, reporting, batching, and DNAm probe filtering
-   [Examples](Examples/): small runnable DNAm, proteomics, and
    metabolomics examples with static toy data and saved example results

## Downstream use

Multi-arm support here applies to Track 1.1.1. Track 3.3 and external validation
need their own adaptations before accepting the combined multi-arm output.
For a summary-only consumer that supports one comparison, select one named
treatment and retain its treatment/reference mapping as provenance. Track 3.3
also needs its outcome model to adjust for the full arm factor: splitting these
tables alone does not complete the mediation extension. Do not pass the combined
table as though it contained one treatment effect per analyte and follow-up.
The CALERIE DunedinPACE covariance helper remains a single-comparison helper
and now rejects combined contrast tables instead of selecting the first row.

## Contact

For any bugs or issues running the program, please feel free to submit
an Issue in this repository or email Will Marella at
[wm2530\@cumc.columbia.edu](mailto:wm2530@cumc.columbia.edu)
