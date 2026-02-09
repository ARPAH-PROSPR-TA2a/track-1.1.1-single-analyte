# FAST Omics WAS Pipeline: Code Walkthrough

## Table of Contents

1. [Main Function](#main-function)
2. [Accepted Inputs](#accepted-inputs)
3. [Input Validation](#input-validation)
4. [Reporting Pipeline: QC and Summary Reports](#reporting-pipeline-qc-and-summary-reports)
5. [Data Preparation](#data-preparation)
6. [Analysis Method Selection](#analysis-method-selection)
7. [Analysis Method 1: Linear Regression (LM)](#analysis-method-1-linear-regression-lm)
8. [Analysis Method 2: Linear Mixed Effects (LME4)](#analysis-method-2-linear-mixed-effects-lme4)
9. [Analysis Method 3: Empirical Bayes Moderation (LIMMA)](#analysis-method-3-empirical-bayes-moderation-limma)
10. [Multiple Testing Correction](#multiple-testing-correction)
11. [Results Output](#results-output)

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
  results = analysis_results,              # Main results table
  omics_summary = omics_report,            # Summary of analytes
  pheno_summary = pheno_report,            # Sample size, demographics
  covariates_summary = covariates_report,  # Covariate distributions
  randomization_summary = randomization_report  # Baseline balance check
)
```

**Stratification**: The pipeline runs three times — once on all participants, once on males only (FEMALE == 0), and once on females only (FEMALE == 1). Male and female subsets are only created if both genders are present in the input data; if only one gender exists, the corresponding elements (male or female) will be NULL.

---

## Accepted Inputs

### Phenotype Data Format

The `pheno` data.frame must contain these columns:

| Column | Type | Notes |
|:---|:---:|:---|
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
|:---|:---:|:---|
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

## Input Validation

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

Sex stratification: creates all/male/female groups (not performed if only one sex is present)

**Returns**: List with `pheno$all`, `pheno$male`, `pheno$female`, plus `pheno$requires_mixed_effects` (TRUE if max FU > 1)

**Omics Data Validation**

Checks:

- ANALYTE_NAME column exists
- All columns (except ANALYTE_NAME) are numeric
- Sample column names match phenotype SAMPLE_IDs (inner join to shared samples)
- Issues warnings for analytes with NA values
- Issues warnings for analytes with near-zero variance

Sex stratification: creates all/male/female groups (not performed if only one sex is present)

**Returns**: List with `omics$all`, `omics$male`, `omics$female` (each with only shared samples)

---

## Reporting Pipeline: QC and Summary Reports

**File**: `reporting_helpers.R`

Four reports are generated for each dataset (all, male, female):

**Phenotype Summary Report**

Reports: Total samples, samples per FU level, number of subjects, treatment group distribution, sex distribution

**Omics Summary Report**

For each analyte: number of non-missing values, mean, median, standard deviation, and value range (min/max)

**Covariate Summary Report**

For each additional covariate: data type, N missing, summary statistics (mean/SD/range for numeric; unique values for factor/logical)

**Randomization/Baseline Balance Report**

For each analyte at baseline (FU=0 only):

- Performs Welch's t-test comparing CONTROL_STATUS groups
- Records: mean difference, Cohen's d, SE, p-value
- Applies Benjamini-Hochberg FDR correction

---

## Data Preparation

After validation, data is prepared for analysis in `main.R`.

**Baseline Separation**

Baseline (FU=0) phenotype samples are separated and their omics values are extracted for use as adjustment covariates. The analysis dataset uses only follow-up measurements (FU > 0).

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
4. Fit model
5. Extract all fixed effect coefficients (including intercept)

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
5. Fit model
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
- Follow-up structure: Any (single or multiple FU)

### Model Specification

LIMMA adapts its model based on follow-up structure:

**Multiple FU (max FU > 1)**:
```
CHANGE ~ CONTROL_STATUS * factor(FU) + FEMALE + additional_covariates
```

With repeated measures handled via:
```r
cor <- duplicateCorrelation(analyte_change, design, block = SUBJECT_ID)
fit <- lmFit(analyte_change, design, 
             block = SUBJECT_ID, 
             correlation = cor$consensus.correlation)
fit <- eBayes(fit)
```

**Single FU (max FU = 1)**:
```
CHANGE ~ CONTROL_STATUS + FEMALE + additional_covariates
```

Without repeated measures (no duplicateCorrelation):
```r
fit <- lmFit(analyte_change, design)
fit <- eBayes(fit)
```

This adaptive approach ensures LIMMA works correctly whether data has repeated measurements or not.

**Note on `analyte_baseline`**: LIMMA does not include `analyte_baseline` as a covariate. With massive DNAm datasets, LIMMA's vectorized approach requires a single design matrix for all analytes, which precludes per-analyte baseline adjustments like those in LM/LME4.

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

After analysis completes and correction is applied, results are returned in this format. All available fixed effects were recorded.

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

- **ANALYTE_NAME**: Feature name
- **COEFFICIENT**: Regression coefficient name
- **EFFECT_SIZE**: Estimated regression coefficient
- **SE**: Standard error
- **P_VALUE**: Two-tailed p-value (unadjusted)
- **BH_P_VALUE**: Benjamini-Hochberg FDR-adjusted p-value


