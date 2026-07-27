# FAST Omics WAS Pipeline: Code Walkthrough (v3)

This walkthrough matches the current pipeline behavior in `main.R`. The pipeline
exposes two public functions: `FAST_omics_WAS()` runs the statistical analyses
(`analysis_change` and `analysis_level`) with parallelization and optional
checkpointing; `FAST_omics_WAS_reports()` generates QC reports independently.
Both standardize inference on `lm()` / `lmerTest::lmer()`.

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
analysis_helpers.R            # Statistical analysis (lm / lmerTest), parallelization, checkpointing
reporting_helpers.R           # QC and summary report generation
plotting_helpers.R            # QQ + volcano plotting from treatment_effects
test_omics_type.R             # omics_type allowlist and silent "other" regression test
test_comprehensive.R          # Functional test suite (2x2 scenario matrix)
test_parallel_checkpoint.R    # Parallelization and checkpointing tests
```

### Function Locations

| File | Key functions |
|------|---------------|
| `main.R` | `FAST_omics_WAS()`, `FAST_omics_WAS_reports()` |
| `validation_helpers.R` | `.validate_omics_type()`, `.validate_pheno()`, `.validate_omics()`, `.validate_dnam_probe_coverage()`, `.subset_omics_list()` |
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
                           checkpoint_batch_size = 2000L)

FAST_omics_WAS_reports <- function(pheno,
                                   omics,
                                   omics_type            = "Proteomics",
                                   additional_covariates = NULL)
```

### FAST_omics_WAS — What it does

1. Resolves `n_cores` (auto-detects if `NULL`), sets `future::plan(multisession)`, restores previous plan on exit
2. Validates `omics_type`, `pheno`, and `omics` (and creates sex strata)
3. For `omics_type == "DNAm"`: loads probe lists, validates coverage, and subsets omics once
4. Runs the analysis twice in parallel across analytes:
   - once with `response_type = "change"`
   - once with `response_type = "level"`
5. Returns:

```r
list(
  analysis_change = <analysis_by_stratum>,
  analysis_level  = <analysis_by_stratum>
)
```

### FAST_omics_WAS_reports — What it does

1. Runs the same validation as `FAST_omics_WAS` (steps 2–3 above)
2. Generates QC reports via `.generate_reports()`
3. Returns:

```r
list(
  all    = <report>,
  male   = <report>,
  female = <report>
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
| TREATMENT_GROUP | factor (or numeric coerced to factor) | Must be binary 0/1 and include both groups |
| FEMALE | factor (or numeric coerced to factor) | Must be binary 0/1 |
| *Additional covariates* | numeric/factor/logical | Provided via `additional_covariates`; samples with NA values are dropped during validation |

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
- `.validate_pheno(pheno, additional_covariates)`:
  - Checks required columns and value constraints (FU must be 0,1,2,... with no gaps)
  - Coerces `FU`, `TREATMENT_GROUP`, `FEMALE` to factors if needed
  - Deduplicates duplicate `SUBJECT_ID` x `FU` pairs (keeps first row)
  - Drops samples with NA values in any additional covariate
  - Drops subjects without both a baseline (FU == 0) and at least one follow-up (FU > 0)
  - Builds `pheno_list = list(all, male, female, requires_mixed_effects)`
- `.validate_omics(omics, pheno_list)`:
  - Validates that all omics measurement columns are numeric
  - Intersects omics sample columns with phenotype `SAMPLE_ID`s
  - Builds `omics_list = list(all, male, female)`

---

## High-Level Pipeline Flow

```
FAST_omics_WAS()
│
├── resolve n_cores; set future::plan(multisession); register on.exit to restore plan
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
└── .generate_reports()  → reports (all/male/female)
```

`.generate_reports()` loops over strata (`all`, `male`, `female`) once and builds
pheno/omics/covariates/randomization summaries — these are independent of response
type and analyte-level model fits.

Inside `.run_stratified_analysis()`:

1. For each stratum (`all`, `male`, `female`): construct checkpoint subdir path (`{checkpoint_dir}/{response_type}/{stratum}/`), run `.perform_analysis(...)`
2. **DNAm only**: add `BH_P_VALUE_FILTERED` to coefficients and treatment effects for the filtered probe set

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
   - If `checkpoint_dir` is set and a batch file already exists: load from disk and skip
   - Otherwise: run `furrr::future_map()` over the batch in parallel, then save the batch result atomically (write to `.tmp`, rename) if `checkpoint_dir` is set
3. After all batches: `do.call(rbind, ...)` assembles the full results

`furrr::future_map()` uses whatever `future::plan()` is active — set by
`FAST_omics_WAS()` before calling into the helpers. With `n_cores = 1` the
plan is `sequential` and there is no parallelization overhead.

BH correction runs once on the full assembled results after all batches
complete, so it always reflects the correct multiple-testing burden regardless
of batch boundaries.

### Model

For each analyte:

```
response ~ TREATMENT_GROUP + FEMALE + analyte_baseline + additional_covariates
```

Notes:

- `FEMALE` is dropped if a stratum contains only one sex.
- Treatment effect is taken from the `TREATMENT_GROUP` coefficient.
- LM is intended for a single FU level. If multiple FU levels are present, it
  warns and uses the maximum FU level for labeling outputs.

### Outputs

- `coefficients`: all coefficients from `summary(lm_fit)$coefficients`
- `treatment_effects`: the `TREATMENT_GROUP` coefficient, labeled with `FU`

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
- Uses `emmeans(fit, ~ TREATMENT_GROUP | FU)` plus `pairs(..., reverse = TRUE)`
  to compute treatment-minus-control contrasts at each FU.

### Outputs

- `coefficients`: all fixed-effect coefficients from `summary(lmer_fit)$coefficients`
- `treatment_effects`: treatment-control contrasts per FU from `emmeans`

---

## Multiple Testing Correction

**File**: `analysis_helpers.R` (`.apply_multiple_testing_correction()`)

Benjamini-Hochberg correction is applied separately by:

- `COEFFICIENT` for the `coefficients` table
- `FU` for the `treatment_effects` table

This produces a `BH_P_VALUE` column in both outputs.

---

## DNAm Probe Sets and Filtered BH

**File**: `main.R` and `analysis_helpers.R` (`.run_stratified_analysis()`)

For `omics_type == "DNAm"`:

1. `main.R` reads two probe lists from `Data/`, validates overlap, and subsets omics
   to the full probe set — this happens once before reports or analysis run.
2. `.run_stratified_analysis()` receives `filtered_probes` and after fitting adds
   `BH_P_VALUE_FILTERED` to `coefficients` and `treatment_effects` for probes in
   the filtered set (and leaves it as `NA` for other probes).

This allows comparing significance under different multiple-testing burdens
without re-running the model fits.

---

## Reporting Pipeline

**File**: `reporting_helpers.R` (`.generate_reports()`)

Reports are generated once (independent of response type) and returned as a
three-element list:

### `$pheno_summary` (study-level)

One row per (FU, FEMALE) cell with subject/sample counts across the full dataset.
Produced by `.create_pheno_data_report()`.

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
- `treatment_effects`: `ANALYTE_NAME`, `FU`, `EFFECT_SIZE`, `SE`, `P_VALUE`, `BH_P_VALUE` (+ `BH_P_VALUE_FILTERED` for DNAm)

### `FAST_omics_WAS_reports()` return value

```r
reports$pheno_summary          # study-level: data.frame, one row per (FU, FEMALE) cell
reports$variable_summaries     # sex-stratified omics + covariate summaries
reports$randomization_reports  # study-level balance checks
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
