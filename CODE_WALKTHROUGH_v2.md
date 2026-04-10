# FAST Omics WAS Pipeline: Code Walkthrough (v2)

## Table of Contents

1. [File Structure](#file-structure)
2. [Main Function](#main-function)
3. [Accepted Inputs](#accepted-inputs)
4. [Input Validation](#input-validation)
5. [Analysis Pipeline Flow](#analysis-pipeline-flow)
6. [Data Preparation](#data-preparation)
7. [Analysis Method Selection](#analysis-method-selection)
8. [Analysis Method 1: Linear Regression (LM)](#analysis-method-1-linear-regression-lm)
9. [Analysis Method 2: Linear Mixed Effects (LME4)](#analysis-method-2-linear-mixed-effects-lme4)
10. [Analysis Method 3: Empirical Bayes Moderation (LIMMA)](#analysis-method-3-empirical-bayes-moderation-limma)
11. [Multiple Testing Correction](#multiple-testing-correction)
12. [Reporting Pipeline](#reporting-pipeline)
13. [Results Output](#results-output)

---

## File Structure

```
main.R                  # Public API: FAST_omics_WAS()
validation_helpers.R    # Input validation and data subsetting
analysis_helpers.R      # Statistical analysis functions
reporting_helpers.R     # QC and summary report generation
test_comprehensive.R    # Test suite covering 2x2 matrix of scenarios
```

### Function Locations

| File | Functions |
|------|-----------|
| `main.R` | `FAST_omics_WAS()` |
| `validation_helpers.R` | `.is_near_zero_variance()`, `.validate_omics_type()`, `.validate_pheno()`, `.validate_omics()`, `.validate_dnam_probe_coverage()`, `.subset_omics_list()` |
| `analysis_helpers.R` | `.apply_multiple_testing_correction()`, `.prepare_analysis_data()`, `.perform_lm_analysis()`, `.perform_lme4_analysis()`, `.perform_limma_analysis()`, `.perform_analysis()`, `.add_filtered_bh_column()`, `.add_filtered_bh_correction()`, `.run_stratified_analysis()` |
| `reporting_helpers.R` | `.create_pheno_data_report()`, `.create_omics_data_report()`, `.create_addx_covariate_report()`, `.create_randomization_report()` |

---

## Main Function

**File**: `main.R`

```r
FAST_omics_WAS <- function(pheno,
                           omics,
                           omics_type = "Proteomics",
                           additional_covariates = NULL)
```

### Function Parameters

| Parameter | Type | Required? | Description |
|:---|:---:|:---:|:---|
| `pheno` | data.frame | YES | Phenotype data with subject/sample info, treatment status, covariates |
| `omics` | data.frame | YES | Omics measurements (analytes × samples) |
| `omics_type` | character | NO | Type of omics data: "Proteomics" (default), "Metabolomics", or "DNAm" |
| `additional_covariates` | character vector | NO | Names of columns in pheno to use as additional covariates |

### Return Value

A list with three elements (all, male, female), each containing:

```r
list(
  coefficients = data.frame(...),         # All model coefficients
  treatment_effects = data.frame(...),    # Treatment effects at each FU
  omics_summary = omics_report,           # Per-analyte summary at baseline (FU=0)
  pheno_summary = pheno_report,           # Subject/sample counts per (FU, FEMALE) cell
  covariates_summary = covariates_report, # Covariate distributions at baseline (FU=0)
  randomization_summary = randomization_report  # Baseline balance check
)
```

**Stratification**: The pipeline runs three times — once on all participants, once on males only (FEMALE == 0), and once on females only (FEMALE == 1). Male and female subsets are only created if both genders are present in the input data.

---

## Accepted Inputs

### Phenotype Data Format

The `pheno` data.frame must contain these columns:

| Column | Type | Notes |
|:---|:---:|:---|
| SAMPLE_ID | character | Unique identifier for each measurement occasion |
| SUBJECT_ID | character | Subject/participant identifier (repeated across FU levels) |
| FU | factor | Follow-up timepoint: levels 0 (baseline), 1, 2, or 3 |
| CONTROL_STATUS | factor | Treatment group: levels 0 (control), 1 (treatment) |
| FEMALE | factor | Sex: levels 0 (male), 1 (female) |
| *Additional covariates* | numeric/factor/logical | Additional variables to adjust for |

### Omics Data Format

The `omics` data.frame must contain:

| Column | Type | Notes |
|:---|:---:|:---|
| ANALYTE_NAME | character | Unique identifier for each feature |
| *Sample IDs* | numeric | One column per sample, named exactly as in pheno$SAMPLE_ID |

---

## Input Validation

**File**: `validation_helpers.R`

### `.validate_omics_type(omics_type)`

Checks that `omics_type` is one of: "Proteomics", "Metabolomics", "DNAm"

### `.validate_pheno(pheno, additional_covariates)`

Checks:

- Required columns exist: SAMPLE_ID, SUBJECT_ID, FU, CONTROL_STATUS, FEMALE
- Data types are correct (converts to factor if needed)
- FU contains 0 (baseline) and at least one FU=1 sample
- CONTROL_STATUS and FEMALE are binary (0/1)
- SAMPLE_ID values are unique
- Duplicate SUBJECT_ID/FU pairs are warned and deduplicated (keeps first occurrence)
- Additional covariates exist and have valid types

**Returns**: List with `$all`, `$male`, `$female` pheno subsets, plus `$requires_mixed_effects` (TRUE if max FU > 1)

### `.validate_omics(omics, pheno_list)`

Checks:

- ANALYTE_NAME column exists
- All data columns are numeric
- Sample columns match phenotype SAMPLE_IDs (inner join)
- Warns for analytes with NA values or near-zero variance

**Returns**: List with `$all`, `$male`, `$female` omics subsets (filtered to shared samples)

### `.validate_dnam_probe_coverage(full_probes, filtered_probes, available_probes)`

For DNAm only: validates that probe lists have overlap with the data. Stops if no overlap; warns if partial overlap.

### `.subset_omics_list(omics_list, analyte_subset)`

Filters omics data to a specific set of analytes (used for DNAm probe subsetting).

---

## Analysis Pipeline Flow

```
FAST_omics_WAS()
│
├── .validate_omics_type()
├── .validate_pheno() → pheno_list (all/male/female + requires_mixed_effects)
├── .validate_omics() → omics_list (all/male/female)
│
└── .run_stratified_analysis()
    │
    ├── [DNAm only] Load probe lists, validate coverage, subset to full probes
    │
    ├── FOR EACH stratum (all, male, female):
    │   ├── Generate reports (pheno, omics, covariates, randomization)
    │   └── .perform_analysis()
    │       ├── Extract baseline data (FU=0)
    │       ├── Filter to post-baseline (FU>0)
    │       └── Dispatch to appropriate method:
    │           ├── DNAm → .perform_limma_analysis()
    │           ├── Non-DNAm, single FU → .perform_lm_analysis()
    │           └── Non-DNAm, multi FU → .perform_lme4_analysis()
    │
    └── [DNAm only] .add_filtered_bh_correction() → adds BH_P_VALUE_FILTERED column
```

---

## Data Preparation

**File**: `analysis_helpers.R`

### `.prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)`

Shared helper used by all three analysis methods. Ensures identical data preparation:

1. **Sample matching**: Find intersection of pheno SAMPLE_IDs and omics column names
2. **Pheno filtering**: Filter pheno to shared samples → `pheno_merged`
3. **Baseline mapping**: For each FU sample, find its subject's baseline column index

**Returns**:
```r
list(
  pheno_merged = ...,         # Pheno filtered to shared samples
  shared_samples = ...,       # Vector of shared sample IDs
  baseline_col_idx = ...,     # Mapping: FU sample → baseline column position
  omics_baseline_matrix = ... # Baseline omics as matrix
)
```

### Change Score Computation

All methods compute change scores: `FU_value - baseline_value`

- **LM/LME4**: Computed per-analyte in the loop
- **LIMMA**: Computed as full matrix subtraction (vectorized)

---

## Analysis Method Selection

**File**: `analysis_helpers.R` → `.perform_analysis()`

```
Is omics_type == "DNAm"?
├─ YES → Use LIMMA (vectorized, empirical Bayes)
└─ NO → Check maximum follow-up level
       ├─ max FU == 1 → Use LM (linear regression)
       └─ max FU > 1 → Use LME4 (linear mixed effects)
```

---

## Analysis Method 1: Linear Regression (LM)

**When Used**: Proteomics/Metabolomics with single follow-up (FU=0 and FU=1 only)

### Model Specification

```
CHANGE ~ CONTROL_STATUS + FEMALE + analyte_baseline + additional_covariates

where:
  CHANGE = FU1_value - baseline_value
  FEMALE included only when both sexes present
```

### Implementation

Per-analyte loop using `lm()`:

1. Call `.prepare_analysis_data()` to get sample matching and baseline mapping
2. For each analyte:
   - Extract FU values and baseline values
   - Compute change score
   - Fit linear model
   - Extract all coefficients (including intercept)
   - Extract treatment effect (CONTROL_STATUS coefficient)

---

## Analysis Method 2: Linear Mixed Effects (LME4)

**When Used**: Proteomics/Metabolomics with multiple follow-ups (max FU > 1)

### Model Specification

```
CHANGE ~ CONTROL_STATUS * FU + FEMALE + analyte_baseline + additional_covariates + (1|SUBJECT_ID)

where:
  CHANGE = FU_value - baseline_value
  FU = categorical follow-up level
  CONTROL_STATUS * FU = main effect + interaction
  (1|SUBJECT_ID) = random intercept per subject
```

### Implementation

Per-analyte loop using `lmer()` with `REML=FALSE`:

1. Call `.prepare_analysis_data()`
2. For each analyte:
   - Compute change scores
   - Fit mixed model
   - Extract coefficients with approximate p-values
   - Extract treatment effects at each FU using `emmeans()`

### P-value Computation

Uses `lmerTest::lmer()` to obtain Satterthwaite degrees of freedom and p-values directly from `summary()$coefficients[, "Pr(>|t|)"]`. `emmeans` automatically inherits the Satterthwaite df from the `lmerTest` model object.

---

## Analysis Method 3: Empirical Bayes Moderation (LIMMA)

**When Used**: DNAm (DNA methylation) — any follow-up structure

### Model Specification

**Single FU**:
```
CHANGE ~ CONTROL_STATUS + FEMALE + additional_covariates
```

**Multiple FU**:
```
CHANGE ~ CONTROL_STATUS * FU_factor + FEMALE + additional_covariates
```

With repeated measures handled via `duplicateCorrelation()`:
```r
cor <- duplicateCorrelation(analyte_change, design, block = SUBJECT_ID)
fit <- lmFit(analyte_change, design, block = SUBJECT_ID, correlation = cor$consensus.correlation)
```

### Implementation

Vectorized (all analytes simultaneously):

1. Call `.prepare_analysis_data()`
2. Compute full change matrix: `omics_values - omics_baseline_merged`
3. Build design matrix based on FU structure
4. Fit using `lmFit()` + `eBayes()`
5. Extract coefficients for all analytes
6. Compute treatment effects using `contrasts.fit()`

### P-value Computation

P-values are extracted directly from `eBayes()` output (`fit$p.value`). These are moderated t-test p-values computed using empirical Bayes shrunken variance estimates.

### Note on Baseline Adjustment

LIMMA does not include `analyte_baseline` as a covariate. With massive DNAm datasets, LIMMA's vectorized approach requires a single design matrix for all analytes, which precludes per-analyte baseline adjustments.

---

## Multiple Testing Correction

**File**: `analysis_helpers.R`

### Standard BH Correction

Applied separately for each coefficient/FU level:

```r
.apply_multiple_testing_correction <- function(results_df, group_col) {
  for (grp in unique(results_df[[group_col]])) {
    idx <- which(results_df[[group_col]] == grp)
    results_df$BH_P_VALUE[idx] <- p.adjust(results_df$P_VALUE[idx], method = "BH")
  }
}
```

### DNAm Filtered BH Correction

For DNAm, an additional column `BH_P_VALUE_FILTERED` is added:

- **Full probes**: All probes get `BH_P_VALUE` (BH across all probes)
- **Filtered probes**: Probes in the filtered list also get `BH_P_VALUE_FILTERED` (BH across filtered probes only)
- **Non-filtered probes**: `BH_P_VALUE_FILTERED = NA`

This allows assessment of significance under different multiple testing burdens while running LIMMA only once (better variance estimation from full probe set).

---

## Reporting Pipeline

**File**: `reporting_helpers.R`

Four reports generated for each stratum:

### `.create_pheno_data_report(pheno_df)`

Returns a data frame with one row per (FU, FEMALE) cell. Columns: FU, FEMALE, N_SUBJECTS, N_CONTROL, N_TREATMENT, N_SAMPLES. Subject-level counts dedupe by SUBJECT_ID; N_SAMPLES is row-level so technical replicates remain visible.

### `.create_omics_data_report(pheno_df, omics_df)`

Per-analyte at baseline (FU=0), as a pre-treatment reference distribution: N_NONMISSING, MEAN, MEDIAN, SD, MIN, MAX

### `.create_addx_covariate_report(pheno_df, covariate_names)`

Per-covariate at baseline (FU=0), as a pre-treatment reference distribution: TYPE, N_NA, SUMMARY (type-specific statistics)

### `.create_randomization_report(pheno_df, omics_df)`

At baseline (FU=0), per-analyte Welch's t-test comparing treatment groups:

- MEAN_DIFFERENCE, COHENS_D, SE, P_VALUE, BH_P_VALUE

---

## Results Output

### Coefficients Table

All model coefficients for all analytes:

```
ANALYTE_NAME    COEFFICIENT              EFFECT_SIZE    SE        P_VALUE     BH_P_VALUE
cg00000029      (Intercept)              0.0234         0.0089    0.0087      0.156
cg00000029      CONTROL_STATUS1          -0.0129        0.0084    0.127       0.941
cg00000029      FEMALE1                  0.0021         0.0076    0.782       0.982
...
```

### Treatment Effects Table

Treatment effects at each follow-up level:

```
ANALYTE_NAME    FU    EFFECT_SIZE    SE        P_VALUE     BH_P_VALUE
cg00000029      1     -0.0129        0.0084    0.127       0.941
cg00000029      2     -0.0071        0.0091    0.435       0.967
...
```

### DNAm-Specific Columns

For DNAm, both tables include an additional column:

```
ANALYTE_NAME    FU    EFFECT_SIZE    SE    P_VALUE    BH_P_VALUE    BH_P_VALUE_FILTERED
cg00000029      1     -0.0129        ...   0.127      0.941         0.523
cg00000103      1     0.0087         ...   0.298      0.967         NA
```

- `BH_P_VALUE`: BH correction across all probes
- `BH_P_VALUE_FILTERED`: BH correction across filtered probes only (NA for non-filtered)

### Sorting Order

Both tables sorted by ANALYTE_NAME (primary), then COEFFICIENT/FU (secondary).
