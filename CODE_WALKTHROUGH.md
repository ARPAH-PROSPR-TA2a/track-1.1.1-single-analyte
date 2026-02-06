# FAST Omics WAS Pipeline: Code Walkthrough

## Table of Contents

1. [Entry Point & Function Signature](#entry-point--function-signature)
2. [Input Acceptance: What Data is Required](#input-acceptance-what-data-is-required)
3. [Input Validation: What is Enforced](#input-validation-what-is-enforced)
4. [Data Preparation](#data-preparation)
5. [Reporting Pipeline: QC and Summary Reports](#reporting-pipeline-qc-and-summary-reports)
6. [Analysis Method Selection](#analysis-method-selection)
7. [Analysis Method 1: Linear Regression (LM)](#analysis-method-1-linear-regression-lm)
8. [Analysis Method 2: Linear Mixed Effects (LME4)](#analysis-method-2-linear-mixed-effects-lme4)
9. [Analysis Method 3: Empirical Bayes Moderation (LIMMA)](#analysis-method-3-empirical-bayes-moderation-limma)
10. [Multiple Testing Correction](#multiple-testing-correction)
11. [Results Output](#results-output)
12. [Stratified Analysis: All, Male, Female](#stratified-analysis-all-male-female)

---

## Entry Point & Function Signature

**File**: `main.R`

```r
FAST_omics_WAS <- function(pheno, 
                           omics, 
                           omics_type = "Proteomics",
                           additional_covariates = NULL)
```

### Function Parameters

| Parameter | Type | Required? | Description |
|-----------|------|-----------|-------------|
| `pheno` | data.frame | YES | Phenotype data with subject/sample info, treatment status, covariates |
| `omics` | data.frame | YES | Omics measurements (analytes × samples) |
| `omics_type` | character | NO | Type of omics data: "Proteomics" (default), "Metabolomics", or "DNAm" |
| `additional_covariates` | character vector | NO | Names of columns in pheno to use as additional covariates |

### Return Value

A list with three elements (all, male, female), each containing:

```r
list(
  results = analysis_results,              # Main results table
  omics_summary = omics_report,            # Summary of analytes
  pheno_summary = pheno_report,            # Sample size, demographics
  covariates_summary = covariates_report,  # Covariate distributions
  randomization_summary = randomization_report  # Baseline balance check
)
```

---

## Input Acceptance: What Data is Required

### Phenotype Data Format

The `pheno` data.frame must contain these columns:

| Column | Type | Notes |
|--------|------|-------|
| SAMPLE_ID | character | Unique identifier for each measurement occasion |
| SUBJECT_ID | character | Subject/participant identifier (repeated across FU levels) |
| FU | numeric | Follow-up timepoint: 0 (baseline), 1, 2, ... |
| CONTROL_STATUS | numeric | Treatment group: 0 (control), 1 (treatment) |
| FEMALE | numeric | Sex: 0 (male), 1 (female) |
| *Any additional covariates specified* | numeric/factor/logical | Additional variables to adjust for |

**Minimal Example**:
```
SAMPLE_ID      SUBJECT_ID  FU  CONTROL_STATUS  FEMALE  agebl  agevis
sample_001     subj_001    0   1               1       55     56
sample_002     subj_001    1   1               1       55     57
sample_003     subj_002    0   0               0       62     63
sample_004     subj_002    1   0               0       62     64
```

### Omics Data Format

The `omics` data.frame must contain:

| Column | Type | Notes |
|--------|------|-------|
| ANALYTE_NAME | character | Unique identifier for each feature (gene, protein, metabolite, CpG site) |
| *Sample IDs* | numeric | One column per sample, named exactly as in pheno$SAMPLE_ID |

**Minimal Example**:
```
ANALYTE_NAME    sample_001   sample_002   sample_003   sample_004
cg00000029      0.602        0.515        0.684        0.598
cg00000103      0.456        0.468        0.412        0.401
cg00000109      0.721        0.735        0.691        0.702
```

---

## Input Validation: What is Enforced

**File**: `validation_helpers.R`

**Omics Type Validation**

Checks that `omics_type` is one of: "Proteomics", "Metabolomics", "DNAm"

**Phenotype Data Validation**

Checks:

- Required columns exist: SAMPLE_ID, SUBJECT_ID, FU, CONTROL_STATUS, FEMALE, and any `additional_covariates`
- Data types are correct
- FU >= 0, CONTROL_STATUS in {0, 1}, FEMALE in {0, 1}
- SAMPLE_ID values are unique
- No duplicate SUBJECT_ID/FU pairs
- Sex stratification: creates all/male/female groups (both genders must be present)

**Returns**: List with `pheno$all`, `pheno$male`, `pheno$female`, plus `pheno$requires_mixed_effects` (TRUE if max FU > 1)

**Omics Data Validation**

Checks:

- ANALYTE_NAME column exists
- All columns (except ANALYTE_NAME) are numeric
- Sample column names match phenotype SAMPLE_IDs (inner join to shared samples)
- Issues warnings for analytes with NA values
- Issues warnings for analytes with near-zero variance

**Returns**: List with `omics$all`, `omics$male`, `omics$female` (each with only shared samples)

---

## Data Preparation

After validation, data is prepared for analysis in `main.R`.

**Baseline Separation**

Baseline (FU=0) phenotype samples are separated and their omics values are extracted for use as adjustment covariates. The analysis dataset uses only follow-up measurements (FU > 0).

---

## Reporting Pipeline: QC and Summary Reports

**File**: `reporting_helpers.R`

Four reports are generated for each dataset (all, male, female):

**Phenotype Summary Report**

Reports: Total samples, samples per FU level, number of subjects, treatment group distribution, sex distribution

**Omics Summary Report**

For each analyte: number of non-missing values, mean, median, standard deviation, and value range (min/max)"

**Covariate Summary Report**

For each additional covariate: data type, N missing, summary statistics (mean/SD/range for numeric; unique values for character)

**Randomization/Baseline Balance Report**

For each analyte at baseline (FU=0 only):

- Performs Welch's t-test comparing CONTROL_STATUS groups
- Records: mean difference, Cohen's d, SE, p-value
- Applies Benjamini-Hochberg FDR correction

---

## Analysis Method Selection

**File**: `analysis_helpers.R`

Automatic selection based on omics type and follow-up structure:

```
Is omics_type == "DNAm"?
├─ YES → Use LIMMA (vectorized analysis with empirical Bayes)
└─ NO → Check maximum follow-up level
       ├─ max FU == 1 → Use LM (linear regression)
       └─ max FU > 1 → Use LME4 (linear mixed effects)
```

---

## Analysis Method 1: Linear Regression (LM)

**File**: `analysis_helpers.R`

**When Used**:

- Omics type: Proteomics or Metabolomics
- Follow-up structure: Only FU=0 and FU=1 (single follow-up)

### Model Specification

```
CHANGE ~ CONTROL_STATUS + FEMALE + analyte_baseline + additional_covariates

where CHANGE = FU_1_value - baseline_value
```

### Implementation

Per-analyte loop:

1. Get baseline values for analyte
2. Get FU values
3. Compute change score
4. Fit: `lm(analyte ~ CONTROL_STATUS + analyte_baseline + covariates, data = model_data)`
5. Extract all fixed effect coefficients (skip intercept)

---

## Analysis Method 2: Linear Mixed Effects (LME4)

**File**: `analysis_helpers.R`

**When Used**:

- Omics type: Proteomics or Metabolomics
- Follow-up structure: Multiple FU (max FU > 1)

### Model Specification

```
CHANGE ~ CONTROL_STATUS * factor(FU) + FEMALE + analyte_baseline + additional_covariates + (1|SUBJECT_ID)

where:
  CHANGE = FU_value - baseline_value for that subject
  factor(FU) = categorical follow-up level
  CONTROL_STATUS * factor(FU) = treatment main effect AND treatment × FU interaction
  (1|SUBJECT_ID) = random intercept per subject
```

Fitted using: `lmer()` with `REML=FALSE` (maximum likelihood estimation)

### Implementation

Per-analyte loop:

1. Pre-allocate data template with columns: analyte, analyte_baseline, SUBJECT_ID, FU, CONTROL_STATUS, covariates
2. Get baseline and FU values
3. Compute change score
4. Fill data template
5. Fit: `lmer(analyte ~ CONTROL_STATUS * factor(FU) + analyte_baseline + covariates + (1|SUBJECT_ID))`
6. Extract coefficients using summary()

### P-value Computation

```r
df_approx <- nrow(model_data) - n_fixed_effects
t_stat <- coefficient / SE
p_value <- 2 * pt(-abs(t_stat), df = df_approx)
```

P-values computed from t-statistics using approximate degrees of freedom (lmer with REML=FALSE doesn't automatically return p-values).

---

## Analysis Method 3: Empirical Bayes Moderation (LIMMA)

**File**: `analysis_helpers.R`

**When Used**:

- Omics type: DNAm (DNA methylation)
- Follow-up structure: Multiple FU (max FU > 1)

### Model Specification

```
CHANGE ~ CONTROL_STATUS * factor(FU) + additional_covariates
```

With repeated measures handled via:
```r
cor <- duplicateCorrelation(analyte_change, design, block = SUBJECT_ID)
fit <- lmFit(analyte_change, design, 
             block = SUBJECT_ID, 
             correlation = cor$consensus.correlation)
fit <- eBayes(fit)
```

### Implementation

**Vectorized change score computation**:

```r
# Pre-compute baseline indices
baseline_idx <- match(pheno_merged$SUBJECT_ID, baseline_subject_ids)
baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], 
                          colnames(omics_baseline_matrix))

# Get baseline values as matrix (n_analytes × n_samples)
omics_baseline_merged <- omics_baseline_matrix[, baseline_col_idx, drop = FALSE]

# Compute all change scores at once (matrix subtraction)
analyte_change <- omics_values - omics_baseline_merged
```

**Vectorized model fitting**:

```r
# Estimate within-subject correlation (one value across all analytes)
cor <- duplicateCorrelation(analyte_change, design, block = pheno_merged$SUBJECT_ID)

# Fit all analytes at once
fit <- lmFit(analyte_change, design, 
             block = pheno_merged$SUBJECT_ID, 
             correlation = cor$consensus.correlation)

# Apply empirical Bayes moderation
fit <- eBayes(fit)
```

**Coefficient extraction**:

```r
for (coef_name in colnames(design)) {
  if (coef_name == "(Intercept)") next
  
  coef_idx <- which(colnames(design) == coef_name)
  
  # Extract across all analytes at once
  effect_sizes <- fit$coefficients[, coef_idx]
  ses <- fit$stdev.unscaled[, coef_idx] * fit$sigma
  p_values <- fit$p.value[, coef_idx]
  
  # Append to results
  for (j in seq_len(n_analytes)) {
    results$ANALYTE_NAME[row_idx] <- omics_df$ANALYTE_NAME[j]
    results$COEFFICIENT[row_idx] <- coef_name
    results$EFFECT_SIZE[row_idx] <- effect_sizes[j]
    results$SE[row_idx] <- ses[j]
    results$P_VALUE[row_idx] <- p_values[j]
    row_idx <- row_idx + 1
  }
}
```

---

## Multiple Testing Correction

**File**: `analysis_helpers.R`

After all analyses complete, global Benjamini-Hochberg FDR correction applied:

```r
bh_p <- p.adjust(p_values, method = "BH")
```

All p-values corrected together (all analytes × all coefficients).

---

## Results Output

After analysis completes and correction is applied, results are returned in this format:

### Results Table Structure

```
ANALYTE_NAME    COEFFICIENT              EFFECT_SIZE       SE          P_VALUE        BH_P_VALUE
cg00000029      analyte_baseline        -0.5644           0.0465      5.8e-29       3.5e-28
cg00000029      agebl                    0.0001            0.0001      0.156         0.941
cg00000029      CONTROL_STATUS          -0.0129            0.0084      0.127         0.941
cg00000029      CONTROL_STATUS:factor(FU)2  0.0058        0.0087      0.508         0.941
cg00000103      analyte_baseline        -0.7285            0.0535      1.2e-34       2.2e-33
...
```

### Sorting Order

Sorted by ANALYTE_NAME (primary), then COEFFICIENT (secondary).

### Column Definitions

| Column | Definition |
|--------|-----------|
| ANALYTE_NAME | Feature name |
| COEFFICIENT | Regression coefficient name |
| EFFECT_SIZE | Estimated regression coefficient |
| SE | Standard error |
| P_VALUE | Two-tailed p-value (unadjusted) |
| BH_P_VALUE | Benjamini-Hochberg FDR-adjusted p-value |

### Coefficient Types Reported

**Treatment effects**:

- `CONTROL_STATUS`: Treatment effect
- `CONTROL_STATUS:factor(FU)N`: Interaction (treatment effect modification at timepoint N)

**Time effects** (LME4 and LIMMA only):

- `factor(FU)N`: Main effect of timepoint N

**Baseline adjustment**:

- `analyte_baseline`: Effect of baseline analyte level on change

**Covariate effects**:

- Any column names specified in `additional_covariates`

---

## Stratified Analysis: All, Male, Female

**File**: `main.R`

The entire pipeline runs three times:
1. All participants
2. Males only (if N >= 20)
3. Females only (if N >= 20)

```r
for (dataset in c("all", "male", "female")) {
  if (is.null(pheno_list[[dataset]])) {
    next  # Skip if N < 20 or no data
  }
  
  # Run validation, reporting, and analysis for this dataset
  outputs[[dataset]] <- list(
    results = analysis_results,
    omics_summary = omics_report,
    pheno_summary = pheno_report,
    covariates_summary = covariates_report,
    randomization_summary = randomization_report
  )
}
```

Each stratified analysis is completely independent with its own results, reports, and quality control checks.
