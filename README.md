---
editor_options: 
  markdown: 
    wrap: 72
---

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

# View association results
head(results$all$results)

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

-   **`$results`**: Data frame of association test results with columns
    for each coefficient, p-values, effect sizes, and confidence
    intervals
-   **`$pheno_summary`**: Data quality report for phenotype data
-   **`$omics_summary`**: Data quality report for omics data
-   **`$covariates_summary`**: Summary of additional covariates (NULL if
    no additional covariates provided)
-   **`$randomization_summary`**: Report on randomization/balance of key
    variables

## Data Format Requirements

### Phenotype Data

Phenotype data must be a data frame with the following **required
columns**:

| Column                  | Type                     | Description                                          | Valid Values                           |
|--------------|--------------|----------------------------|-----------------|
| `SAMPLE_ID`             | Character                | Unique identifier for each sample                    | Unique, no duplicates                  |
| `SUBJECT_ID`            | Character                | Subject identifier for tracking across follow-ups    | Any unique or repeated value           |
| `FU`                    | Factor                   | Follow-up timepoint                                  | Levels: 0 (baseline), 1-3 (follow-ups) |
| `FEMALE`                | Factor                   | Sex indicator                                        | Levels: 0 (male), 1 (female)           |
| `CONTROL_STATUS`        | Factor                   | Treatment assignment                                 | Levels: 0 (control), 1 (treatment)     |
| *Additional Covariates* | *Factor/Numeric/Logical* | *Optional columns for additional analysis variables* | *Any valid values*                     |

**Data requirements:**

-   Must contain at least one baseline sample (`FU == 0`) and one
    follow-up sample (`FU >= 1`)
-   Must contain both males (`FEMALE == 0`) and females (`FEMALE == 1`),
    or gender stratification will be skipped
-   Must contain both control (`CONTROL_STATUS == 0`) and treatment
    (`CONTROL_STATUS == 1`) groups
-   No duplicate `SAMPLE_ID` values
-   `FU`, `FEMALE`, and `CONTROL_STATUS` must be factors (or numeric
    values that will be automatically converted to factors)

**Additional covariates (optional):**

-   Any columns not listed above can be included as additional
    covariates
-   Must be numeric, factor, or logical
-   NA values are allowed but will reduce analysis sample size with a
    warning
-   Specified using the `additional_covariates` parameter

### Omics Data

Omics data must be a data frame with:

| Element            | Description                                                                |
|-----------------|-------------------------------------------------------|
| `ANALYTE_NAME`     | Column containing identifiers for each analyte/protein/site (first column) |
| Additional columns | Named by `SAMPLE_ID` (must match `SAMPLE_ID` values from phenotype data)   |
| Data values        | All analyte measurements must be numeric                                   |

**Data requirements:** - Column names (except `ANALYTE_NAME`) must match
`SAMPLE_ID` values from phenotype data - All data columns must be
numeric - Analytes with NA values or near-zero variance will generate a
warning but are included in analysis - Any omics samples not in
phenotype data are automatically excluded - Any phenotype samples not in
omics data are automatically excluded

**Example structure:**

```         
ANALYTE_NAME  SAMPLE_001  SAMPLE_002  SAMPLE_003
Protein_A     0.45        0.52        0.48
Protein_B     1.23        1.31        1.19
...
```

## Analysis Methods

The pipeline automatically selects the appropriate analysis method based
on your data structure:

-   **Linear Model (LM)**: Used when all follow-up values are 0-1
    (baseline and single follow-up only)
-   **Linear Mixed-Effects Model (LME4)**: Used when multiple follow-ups
    are present (FU values of 2-3)
-   **Limma**: Used for DNA methylation data with appropriate precision
    weighting

## More Information

For comprehensive examples and detailed documentation of the analysis
pipeline, see [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md).

For any bugs or issues running the program, please feel free to submit
an Issue in this repository or email Will Marella at
[wm2530\@cumc.columbia.edu](mailto:wm2530@cumc.columbia.edu){.email}
