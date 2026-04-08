# Omics-Wide Association Study Pipeline

This repository, written in R, provides a reliable implementation of
omics-wide association studies for Proteomics, Metabolomics, and DNA
Methylation. It supports multiple analysis methods (linear models,
linear mixed-effects models, and limma) and automatically stratifies
results by gender when data permits. It also provides detailed data
documentation.

## Installation

Clone the repository:

``` bash
git clone https://github.com/ARPAH-PROSPR-TA2a/track-1.1.1-single-analyte.git
```

In your R script, set your working directory to the package directory
and source the main module:

``` r
setwd("path/to/track-1.1.1-single-analyte")
source("main.R")
```

The pipeline automatically sources all necessary helper functions from
the same directory.

## FAST_omics_WAS()

The main function for conducting omics-wide association analyses.

### Quick Example

``` r
# Load your phenotype and omics data
pheno <- read.csv("my_phenotype.csv")
omics <- read.csv("my_omics_data.csv")

# Run the association analysis
results <- FAST_omics_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi")
)

# View results
head(results$all$coefficients)
head(results$all$treatment_effects)

# Access data quality summaries
results$all$pheno_summary
results$all$omics_summary
```

### Parameters

-   **`pheno`** (data.frame): Phenotype data with required columns (see
    Data Format section below)
-   **`omics`** (data.frame): Omics data with required structure (see
    Data Format section below)
-   **`omics_type`** (character): Type of omics data being analyzed.
    Options: `"Proteomics"`, `"DNAm"`, `"Metabolomics"`. Default:
    `"Proteomics"`
-   **`additional_covariates`** (character vector, optional): Names of
    additional covariates to include in regression models. Must be
    column names in the `pheno` data frame. These columns must be
    numeric, factor, or logical.

### Return Value

A list containing results stratified by gender:

-   **`$all`**: Results from the full dataset
-   **`$male`**: Results from male subset (NULL if no males in data)
-   **`$female`**: Results from female subset (NULL if no females in
    data)

Each stratum contains:

-   **`$coefficients`**: Data frame of all model coefficients for each
    analyte (ANALYTE_NAME, COEFFICIENT, EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE)
-   **`$treatment_effects`**: Data frame of treatment effects at each
    follow-up level (ANALYTE_NAME, FU, EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE)
-   **`$pheno_summary`**: Sample size, sex, treatment, and timepoint breakdowns
-   **`$omics_summary`**: Per-analyte summary statistics
-   **`$covariates_summary`**: Summary of additional covariates (NULL if
    no additional covariates provided)
-   **`$randomization_summary`**: Baseline balance check per analyte
    (Welch's t-test comparing treatment groups)

## Data Format Requirements

### Phenotype Data

Phenotype data must be a data frame with the following **required
columns**:

| Column                  | Type                     | Description                                          | Valid Values                            |
|----------------|----------------|------------------------|------------------|
| `SAMPLE_ID`             | Character                | Unique identifier for each sample                    | Unique, no duplicates                   |
| `SUBJECT_ID`            | Character                | Subject identifier for tracking across follow-ups    | All SUBJECT_ID\*FU pairs must be unique |
| `FU`                    | Factor                   | Follow-up timepoint                                  | Levels: 0 (baseline), 1-3 (follow-ups)  |
| `FEMALE`                | Factor                   | Sex indicator                                        | Levels: 0 (male), 1 (female)            |
| `CONTROL_STATUS`        | Factor                   | Treatment assignment                                 | Levels: 0 (control), 1 (treatment)      |
| *Additional Covariates* | *Factor/Numeric/Logical* | *Optional columns for additional analysis variables* | *Any valid values*                      |

**Data requirements:**

-   Must contain at least one baseline sample (`FU == 0`) and one
    follow-up sample (`FU >= 1`)
-   Must contain both males (`FEMALE == 0`) and females (`FEMALE == 1`),
    otherwise gender stratification will be skipped
-   Must contain both control (`CONTROL_STATUS == 0`) and treatment
    (`CONTROL_STATUS == 1`) groups
-   No duplicate `SAMPLE_ID` values
-   All `SUBJECT_ID*FU` pairs must be unique
-   `FU`, `FEMALE`, and `CONTROL_STATUS` must be factors (or numeric
    values that will be automatically converted to factors)

**Example structure:**

```         
SAMPLE_ID      SUBJECT_ID  FU  CONTROL_STATUS  FEMALE  agebl  agevis
sample_001     subj_001    0   1               1       55     56
sample_002     subj_001    1   1               1       55     57
sample_003     subj_002    0   0               0       62     63
sample_004     subj_002    1   0               0       62     64
```

**Additional covariates (optional):**

-   Any columns not listed above can be included as additional
    covariates
-   Must be numeric, factor, or logical
-   NA values are allowed but will reduce analysis sample size with a
    warning
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

The pipeline automatically selects the appropriate analysis method based
on omics type and follow-up structure:

```         
Is omics_type == "DNAm"?
├─ YES → Use LIMMA (vectorized analysis with empirical Bayes)
└─ NO → Check maximum follow-up level
       ├─ max FU == 1 → Use LM (linear regression)
       └─ max FU > 1 → Use LME4 (linear mixed effects)
```

## More Information

For comprehensive examples and detailed documentation of the analysis
pipeline, see [CODE_WALKTHROUGH_v2.md](CODE_WALKTHROUGH_v2.md).

## Plotting Results

The pipeline includes plotting functions for visualizing treatment effects.
Source `plotting_helpers.R` and call `generate_all_plots()` on your results:

``` r
source("plotting_helpers.R")

results <- FAST_omics_WAS(pheno, omics, omics_type = "DNAm",
                          additional_covariates = c("age", "bmi"))

generate_all_plots(results)                        # saves to Figures/
generate_all_plots(results, figures_dir = "my_dir") # saves to my_dir/
```

This auto-detects all dimensions of your results (strata, FU levels,
filtered probes) and outputs multi-page PDFs. For DNAm results, both
full and filtered probe set plots are generated.

Two plot types are produced:

-   **QQ plots**: Observed vs expected p-value distributions (via `qqman`)
-   **Volcano plots**: Effect size vs significance, with Up/Down/Not
    significant coloring based on BH-corrected p-values (via `ggplot2`)

Requires packages: `qqman`, `ggplot2`.

## Contact

For any bugs or issues running the program, please feel free to submit
an Issue in this repository or email Will Marella at
[wm2530\@cumc.columbia.edu](mailto:wm2530@cumc.columbia.edu)
