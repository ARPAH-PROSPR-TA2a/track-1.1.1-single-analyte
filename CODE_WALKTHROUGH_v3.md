# FAST Omics WAS Pipeline: Code Walkthrough (v3)

This walkthrough matches the current pipeline behavior in `main.R`. The pipeline
exposes two public functions: `FAST_omics_WAS()` runs the statistical analyses
(`analysis_change` and `analysis_level`) with parallelization and optional
checkpointing; `FAST_omics_WAS_reports()` generates QC reports independently.
Both standardize inference on `lm()` / `lmerTest::lmer()`.
Both use `TREATMENT_GROUP`: `0` for control and `1` for treatment. Additional
intervention arms use codes `2`, `3`, and so on, and require
[`treatment_labels`](INPUTS_OUTPUTS.md#treatment_labels).

## Table of Contents

1. [File Structure](#file-structure)
2. [Main Functions](#main-functions)
3. [Accepted Inputs](#accepted-inputs)
4. [Input Validation](#input-validation)
5. [High-Level Pipeline Flow](#high-level-pipeline-flow)
6. [Response Types: `change` vs `level`](#response-types-change-vs-level)
7. [Model Selection](#model-selection)
8. [Analysis Method 1: Linear Regression (LM)](#analysis-method-1-linear-regression-lm)
9. [Analysis Method 2: Linear Mixed Effects (LME4)](#analysis-method-2-linear-mixed-effects-lme4)
10. [Parallelization and Checkpointing](#parallelization-and-checkpointing)
11. [Multiple Testing Correction](#multiple-testing-correction)
12. [DNAm Probe Sets and Filtered BH](#dnam-probe-sets-and-filtered-bh)
13. [Reporting Pipeline](#reporting-pipeline)
14. [Results Output](#results-output)

---

## File Structure

```
main.R                        # Public API: FAST_omics_WAS(), FAST_omics_WAS_reports()
validation_helpers.R          # Input validation and data subsetting
treatment_helpers.R           # Arm mapping, contrast identity and metadata
analysis_helpers.R            # Statistical analysis (lm / lmerTest), parallelization, checkpointing
reporting_helpers.R           # QC and summary report generation
plotting_helpers.R            # QQ + volcano plotting from treatment_effects
test_comprehensive.R          # Functional test suite (2x2 scenario matrix)
test_parallel_checkpoint.R    # Parallelization and checkpointing tests
```

### Function Locations

| File | Key functions |
|------|---------------|
| `main.R` | `FAST_omics_WAS()`, `FAST_omics_WAS_reports()` |
| `validation_helpers.R` | `.validate_omics_type()`, `.validate_pheno()`, `.validate_omics()`, `.validate_dnam_probe_coverage()`, `.subset_omics_list()` |
| `treatment_helpers.R` | Treatment encoding, explicit control contrasts and metadata |
| `analysis_helpers.R` | `.prepare_analysis_data()`, `.perform_lm_analysis()`, `.perform_lme4_analysis()`, `.perform_analysis()`, `.run_stratified_analysis()`, `.apply_multiple_testing_correction()`, `.add_filtered_bh_correction()` |
| `reporting_helpers.R` | `.generate_reports()`, `.create_pheno_data_report()`, `.create_omics_data_report()`, `.create_addx_covariate_report()`, `.create_randomization_report()` |

---

## Main Functions

**File**: `main.R`

```r
FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type            = "Proteomics",
                           additional_covariates = NULL,
                           n_cores               = NULL,
                           checkpoint_dir        = NULL,
                           checkpoint_batch_size = 2000L,
                           verbose               = FALSE,
                           treatment_labels      = NULL)

FAST_omics_WAS_reports <- function(pheno,
                                   omics,
                                   omics_type            = "Proteomics",
                                   additional_covariates = NULL,
                                   verbose               = FALSE,
                                   treatment_labels      = NULL)
```

### FAST_omics_WAS — What it does

1. Resolves `n_cores` (auto-detects if `NULL`, falling back to one worker if unavailable), sets `future::multicore` for multiple workers or `future::sequential` for one, and restores the previous plan on exit
2. Validates `omics_type`, `pheno`, and `omics` (and creates sex strata)
3. For `omics_type == "DNAm"`: loads probe lists, validates coverage, and subsets omics once
4. Runs the analysis twice in parallel across analytes:
   - once with `response_type = "change"`
   - once with `response_type = "level"`
5. Returns:

```r
list(
  analysis_change = <analysis_by_stratum>,
  analysis_level  = <analysis_by_stratum>,
  metadata       = <schema_and_treatment_mapping>
)
```

### FAST_omics_WAS_reports — What it does

1. Runs the same validation as `FAST_omics_WAS` (steps 2–3 above)
2. Generates QC reports via `.generate_reports()`
3. Returns:

```r
list(
  pheno_summary         = <study_counts>,
  arm_counts            = <separate_arm_counts>,
  variable_summaries    = <summaries_by_stratum>,
  randomization_reports = <baseline_balance_checks>,
  metadata              = <schema_and_treatment_mapping>
)
```

Reports and analysis are intentionally separate: analyses can be parallelized
and checkpointed across many cores while reports are generated independently.

---

## Accepted Inputs

### Phenotype (`pheno`)

Required columns:

| Column | Type | Notes |
|:---|:---:|:---|
| SAMPLE_ID | character | Unique identifier for each measurement occasion |
| SUBJECT_ID | character | Subject/participant identifier (repeated across FU levels) |
| FU | factor (or numeric coerced to factor) | Must contain 0 (baseline) and at least one 1 (follow-up) |
| TREATMENT_GROUP | factor (or numeric coerced to factor) | `0` = control, `1` = treatment; assignment is constant within subject |
| FEMALE | factor (or numeric coerced to factor) | Must be binary 0/1 |
| *Additional covariates* | numeric/factor/logical | Provided via `additional_covariates`; samples with NA values are dropped during validation |

If you wish to add additional intervention arms, use codes `2`, `3`, and so on
in `TREATMENT_GROUP` and supply
[`treatment_labels`](INPUTS_OUTPUTS.md#treatment_labels) for every group,
including control.

### Omics (`omics`)

| Column | Type | Notes |
|:---|:---:|:---|
| ANALYTE_NAME | character | Unique identifier for each feature |
| Sample ID columns | numeric | One column per sample, named exactly as in `pheno$SAMPLE_ID` |

---

## Input Validation

**File**: `validation_helpers.R`

- `.validate_omics_type(omics_type)`:
  - Enforces `omics_type %in% c("DNAm", "Proteomics", "Metabolomics", "other")`
  - Prints reminders about expected preprocessing conventions for omics-specific types; `"other"` is silent
- `.validate_pheno(pheno, additional_covariates, treatment_labels)`:
  - Checks required columns and value constraints (FU must be 0,1,2,... with no gaps)
  - Coerces `FU`, `TREATMENT_GROUP`, `FEMALE` to factors if needed
  - Checks constant treatment assignment within each subject
  - Validates the label map against all observed codes before filtering or
    stratification: unique nonblank keys/labels, no missing or extra codes
  - Requires labels for multiple active arms; defaults to Control/Treatment
    for binary input
  - Uses the original codes as factor levels with control `0` as reference;
    no renumbering or fixed arm limit, and the caller's input is unchanged
  - Carries the study treatment mapping across all strata, even when an arm is
    absent in a particular stratum
  - Deduplicates duplicate `SUBJECT_ID` x `FU` pairs (keeps first row)
  - Drops samples with NA values in any additional covariate
  - Drops subjects without both a baseline (FU == 0) and at least one follow-up (FU > 0)
  - Builds phenotype strata, the mixed-model indicator, and the treatment mapping
- `.validate_omics(omics, pheno_list)`:
  - Validates that all omics measurement columns are numeric
  - Intersects omics sample columns with phenotype `SAMPLE_ID`s
  - Builds `omics_list = list(all, male, female)`

---

## High-Level Pipeline Flow

```
FAST_omics_WAS()
│
├── resolve n_cores; set future::plan(multicore or sequential); register on.exit to restore plan
├── .validate_omics_type()
├── .validate_pheno()  → pheno_list (all/male/female + requires_mixed_effects)
├── .validate_omics()  → omics_list (all/male/female)
│
├── [DNAm only] load probe lists; validate overlap; subset omics_list to full probe set
│
├── .run_stratified_analysis(..., "change", checkpoint_dir, batch_size)  → results$analysis_change
└── .run_stratified_analysis(..., "level",  checkpoint_dir, batch_size)  → results$analysis_level

FAST_omics_WAS_reports()
│
├── .validate_omics_type()
├── .validate_pheno()  → pheno_list
├── .validate_omics()  → omics_list
│
├── [DNAm only] load probe lists; validate overlap; subset omics_list to full probe set
│
└── .generate_reports()  → study reports and stratified variable summaries
```

`.generate_reports()` loops over strata (`all`, `male`, `female`) once and builds
pheno/omics/covariates/randomization summaries — these are independent of response
type and analyte-level model fits.

Inside `.run_stratified_analysis()`:

1. For each stratum (`all`, `male`, `female`): construct checkpoint subdir path (`{checkpoint_dir}/{response_type}/{stratum}/`), run `.perform_analysis(...)`
2. **DNAm only**: after batch assembly, add `BH_P_VALUE_FILTERED` to coefficients and treatment effects for the filtered probe set, and logical `IN_FILTERED_SET` to treatment effects only

---

## Response Types: `change` vs `level`

**File**: `analysis_helpers.R` (`.perform_lm_analysis()` and `.perform_lme4_analysis()`)

For each post-baseline row in `pheno` (FU != 0), the pipeline finds the
matching baseline sample for the same `SUBJECT_ID` (FU == 0) and pulls the
baseline analyte value.

Let:

- `FU_value` = omics measurement at a follow-up timepoint for that sample
- `baseline_value` = omics measurement at baseline for that subject

Then:

- **`change` analysis** uses: `response = FU_value - baseline_value`
- **`level` analysis** uses: `response = FU_value`

In both analyses, the model includes `analyte_baseline = baseline_value` as
an adjustment covariate.

For the same fitted rows and specification, subtracting baseline from the
response only shifts its coefficient. The baseline-adjusted treatment
contrasts are identical between `change` and `level` parameterizations.

---

## Model Selection

**File**: `analysis_helpers.R` (`.perform_analysis()`)

Within a stratum, model choice depends only on the maximum follow-up level in
the post-baseline data:

```
Check maximum follow-up level
├─ max FU == 1 → Use LM (lm)
└─ max FU > 1 → Use LME4 (lmerTest::lmer)
```

---

## Analysis Method 1: Linear Regression (LM)

**File**: `analysis_helpers.R` (`.perform_lm_analysis()`)

Used when there is a single post-baseline follow-up level (max FU == 1).

### Parallelization and Checkpointing

Both `.perform_lm_analysis()` and `.perform_lme4_analysis()` use the same
execution strategy:

1. Analytes are split into batches of `checkpoint_batch_size` (default 2000)
2. For each batch:
   - If `checkpoint_dir` is set and a batch file exists: load from disk and skip
   - Otherwise: run `furrr::future_map()` over the batch in parallel, then save the batch result atomically (write to `.tmp`, rename) if `checkpoint_dir` is set
3. After all batches: `do.call(rbind, ...)` assembles the full results

`furrr::future_map()` uses whatever `future::plan()` is active — set by
`FAST_omics_WAS()` before calling into the helpers. With `n_cores = 1` the
plan is `sequential` and there is no parallelization overhead.

BH correction runs once on the full assembled results after all batches
complete, so it always reflects the correct multiple-testing burden regardless
of batch boundaries.

Resume with the same input data, model settings, and batch layout. Use a fresh
checkpoint directory when any of these change or when upgrading across this
result-schema change. You can change `n_cores` or display labels
in `treatment_labels` on resume; labels are applied after the batches are
assembled and do not require refitting.

### Model

For each analyte:

```
response ~ TREATMENT_GROUP + FEMALE + analyte_baseline + additional_covariates
```

Notes:

- `FEMALE` is dropped if a stratum contains only one sex.
- `TREATMENT_GROUP` is the categorical arm factor using the original codes. All arms enter one joint fit;
  each active-minus-control effect is extracted separately.
- With default treatment coding and a full-rank design, the arm coefficients
  are already these contrasts. The extraction helper reuses them directly;
  other parameterizations fall back to explicit `emmeans` contrasts.
- LM is intended for a single FU level. If multiple FU levels are present, it
  warns and uses the maximum FU level for labeling outputs.

### Outputs

- `coefficients`: all coefficients from `summary(lm_fit)$coefficients`
- `treatment_effects`: one row per analyte × FU × active/control contrast

---

## Analysis Method 2: Linear Mixed Effects (LME4)

**File**: `analysis_helpers.R` (`.perform_lme4_analysis()`)

Used when there are multiple post-baseline follow-up levels (max FU > 1).

### Model

For each analyte:

```
response ~ TREATMENT_GROUP * FU + FEMALE + analyte_baseline + additional_covariates + (1|SUBJECT_ID)
```

Implementation details:

- Uses `lmerTest::lmer(..., REML = FALSE)` to obtain Satterthwaite p-values in
  `summary(fit)$coefficients[, "Pr(>|t|)"]`.
- Uses `emmeans` with explicit active-minus-control contrasts at each FU,
  with `adjust = "none"`. This includes each arm's main effect and FU
  interaction, with their covariance, without adding active-versus-active tests.

### Outputs

- `coefficients`: all fixed-effect coefficients from `summary(lmer_fit)$coefficients`
- `treatment_effects`: active/control contrasts per FU from `emmeans`, with
  `CONTRAST`, display labels `TREATMENT`/`REFERENCE`, and `STATUS`

The expected contrast grid is retained when inference fails. Successful rows
have `STATUS = "ok"`; nonestimable contrasts use `"not_estimable"`, and fit
failures use `"fit_error"`. Missing inference is reported as `NA`.

---

## Multiple Testing Correction

**File**: `analysis_helpers.R` (`.apply_multiple_testing_correction()`)

Benjamini-Hochberg correction is applied separately by:

- `COEFFICIENT` for the `coefficients` table
- `FU` and `CONTRAST` for the `treatment_effects` table

Corrections are separate within each analysis and stratum. This produces a
`BH_P_VALUE` column in both outputs. Filtered DNAm correction uses the same
grouping restricted to the pre-specified probe set. Raw contrast p-values do
not receive an additional adjustment across arms.

---

## DNAm Probe Sets and Filtered BH

**File**: `main.R` and `analysis_helpers.R` (`.run_stratified_analysis()`)

For `omics_type == "DNAm"`:

1. `main.R` reads two probe lists from `Data/`, validates overlap, and subsets omics
   to the full probe set — this happens once before reports or analysis run.
2. `.run_stratified_analysis()` receives `filtered_probes` and after batch assembly adds
   `BH_P_VALUE_FILTERED` to `coefficients` and `treatment_effects` for probes in
   the filtered set (and leaves it as `NA` for other probes).
3. DNAm `treatment_effects` also receives logical `IN_FILTERED_SET`, computed as
   `ANALYTE_NAME %in% filtered_probes`. Membership is independent of `STATUS` and
   missing inference, so it is the same for a probe across follow-ups, contrasts,
   analyses, and strata. This column is absent from coefficient tables and
   non-DNAm results.

This allows comparing significance under different multiple-testing burdens
without re-running the model fits. `BH_P_VALUE_FILTERED = NA` can mean either
nonmembership or a missing p-value; `IN_FILTERED_SET` distinguishes these cases
without changing the models or BH correction.

---

## Reporting Pipeline

**File**: `reporting_helpers.R` (`.generate_reports()`)

Reports are generated once (independent of response type), with the shared
treatment mapping in `$metadata`:

### `$pheno_summary` (study-level)

One row per (FU, FEMALE) cell with subject/sample counts across the full dataset.
Produced by `.create_pheno_data_report()`.
`N_TREATMENT` retains its pooled count of all active subjects.

### `$arm_counts` (study-level)

One row per FU × FEMALE × arm, including zero cells, with original `ARM` codes, display
`TREATMENT` labels, `IS_CONTROL`, `N_SUBJECTS`, and `N_SAMPLES`.

### `$variable_summaries`

Sex-stratified (`all`, `male`, `female`). Within each stratum, `.generate_reports()`
iterates over all FU × TREATMENT_GROUP (Tx) cells present in the data. For each
non-empty cell it produces:

- `omics_FU{n}_Tx{m}`: per-analyte summary (N_NONMISSING, MEAN, MEDIAN, SD, MIN, MAX)
  for the samples in that cell — produced by `.create_omics_data_report(sample_ids, omics_df)`
- `covariates_FU{n}_Tx{m}`: covariate summary for those samples (COVARIATE_NAME, TYPE,
  N_NA, SUMMARY list-column) — produced by `.create_addx_covariate_report()`; omitted
  if no additional covariates

Cells with no samples are skipped entirely. Keys are generated dynamically from the
actual FU levels in the data, so the structure adapts to single or multiple follow-ups.
Tx codes use the same mapping as analysis output: control `Tx0`, then active
`Tx1`, `Tx2`, ... .

### `$randomization_reports` (study-level, not sex-stratified)

- `omics_report`: per-analyte baseline balance check — `.create_randomization_report()`
  (Welch's t-test, MEAN_DIFFERENCE, COHENS_D, SE, P_VALUE, BH_P_VALUE)
- `covariate_report`: baseline balance for `FEMALE` (if both sexes present) and
  additional covariates — `.create_pheno_randomization_report()`. Test selection:
  - **numeric**: Welch's t-test
  - **logical**: chi-squared
  - **factor**: chi-squared if min cell count ≥ 5, Fisher's exact otherwise
    (simulation enabled for tables larger than 2×2)
  
  Columns: `VARIABLE`, `TYPE`, `TEST`, `STATISTIC`, `MIN_CELL_COUNT`, `P_VALUE`,
  `SUMMARY_CONTROL`, `SUMMARY_TREATMENT`

Both tables repeat comparisons for each active arm against control and include
`CONTRAST`, `TREATMENT`, and `REFERENCE`. Omics BH correction is applied within
each contrast.

---

## Results Output

### `FAST_omics_WAS()` return value

```r
results$analysis_change$all
results$analysis_change$male
results$analysis_change$female

results$analysis_level$all
results$analysis_level$male
results$analysis_level$female
```

Within each stratum (example: `results$analysis_change$all`):

```r
list(
  coefficients      = data.frame(...),  # All fixed-effect coefficients
  treatment_effects = data.frame(...)   # Treatment-control effects per FU
)
```

Key columns:

- `coefficients`: `ANALYTE_NAME`, `COEFFICIENT`, `N_OBS`, `EFFECT_SIZE`, `SE`, `P_VALUE`, `BH_P_VALUE` (+ `BH_P_VALUE_FILTERED` for DNAm)
- `treatment_effects`: `ANALYTE_NAME`, `FU`, `CONTRAST`, `TREATMENT`, `REFERENCE`,
  `EFFECT_SIZE`, `SE`, `P_VALUE`, `BH_P_VALUE`, `STATUS`
  (+ `BH_P_VALUE_FILTERED` and logical `IN_FILTERED_SET` for DNAm)

`results$metadata` records schema version `"2.0"`, the normalized full
`treatment_labels` map, the `treatments` comparison mapping (`ARM`, `TREATMENT`,
`REFERENCE`, `CONTRAST`), and the BH policy. `REFERENCE` uses the control display
label; `CONTRAST` retains the input codes, e.g. `2_vs_0`. Separate data frames are available without refitting:

```r
te <- results$analysis_change$all$treatment_effects
effects_by_treatment <- split(te, te$TREATMENT)
```

Consumers that assume one row per analyte/FU must select a contrast explicitly.
Track 3.3 and external validation require separate adaptations before accepting
combined multi-arm results; preserve the selected mapping as provenance.
Track 3.3 also needs full-arm adjustment in its outcome model before the agreed
shared-beta mediation analysis is available.

### `FAST_omics_WAS_reports()` return value

```r
reports$pheno_summary          # study-level: data.frame, one row per (FU, FEMALE) cell
reports$arm_counts             # study-level: separate counts for each arm
reports$variable_summaries     # sex-stratified omics + covariate summaries
reports$randomization_reports  # study-level balance checks
reports$metadata               # schema version and treatment mapping
```

`reports$variable_summaries` (example with FU=0,1 and both treatment groups):

```r
reports$variable_summaries$all$omics_FU0_Tx0       # data.frame, per-analyte stats
reports$variable_summaries$all$omics_FU0_Tx1
reports$variable_summaries$all$omics_FU1_Tx0
reports$variable_summaries$all$omics_FU1_Tx1
reports$variable_summaries$all$covariates_FU0_Tx0  # data.frame or NULL
# ... etc.; male/female follow the same structure
```

`reports$randomization_reports` (study-level, not sex-stratified):

```r
list(
  omics_report     = data.frame(...),      # Per-analyte balance checks (t-test)
  covariate_report = data.frame(...)|NULL  # FEMALE + covariate balance checks
)
```
