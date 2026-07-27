# FAST Omics WAS: Inputs and Outputs

---

## Inputs

### `pheno` — Phenotype Data

A data frame with one row per sample (one subject can appear multiple times, once per follow-up timepoint).

**Required columns:**

| Column | Type | Constraints |
|:---|:---:|:---|
| `SAMPLE_ID` | character | Unique across the entire data frame — no duplicates |
| `SUBJECT_ID` | character | Repeated across timepoints for the same person; every `SUBJECT_ID × FU` pair must be unique |
| `FU` | factor (or numeric, auto-coerced) | Must contain `0` (baseline) and at least one follow-up (`1`, `2`, …); values must be consecutive integers with no gaps |
| `TREATMENT_GROUP` | factor (or numeric, auto-coerced) | Binary: `0` = control, `1` = treatment; both groups must be present |
| `FEMALE` | factor (or numeric, auto-coerced) | Binary: `0` = male, `1` = female |
| *Additional covariates* | numeric / factor / logical | Any extra columns named in `additional_covariates`; NAs are allowed but reduce sample size |

**Example:**

```
SAMPLE_ID    SUBJECT_ID   FU   TREATMENT_GROUP   FEMALE   age   bmi
s001         subj_01      0    1                1        55    27.1
s002         subj_01      1    1                1        55    27.1
s003         subj_02      0    0                0        62    24.8
s004         subj_02      1    0                0        62    24.8
s005         subj_03      0    1                0        48    31.2
s006         subj_03      1    1                0        48    31.2
s007         subj_03      2    1                0        48    31.2
```

---

### `omics` — Omics Data

A data frame with one row per analyte (protein, metabolite, CpG site, etc.) and one column per sample.

| Column | Type | Notes |
|:---|:---:|:---|
| `ANALYTE_NAME` | character | Unique identifier for each feature; must be the first column |
| Sample ID columns | numeric | One column per sample, named exactly as values in `pheno$SAMPLE_ID` |

Samples present in omics but not in pheno are dropped silently. Samples present in pheno but not in omics are dropped with a warning.

**Example:**

```
ANALYTE_NAME   s001    s002    s003    s004    s005    s006    s007
protein_A      1.204   1.318   0.987   1.052   1.401   1.289   1.335
protein_B      0.523   0.489   0.601   0.578   0.512   0.531   0.498
protein_C      2.017   2.103   1.889   1.944   2.201   2.088   2.156
```

---

### `omics_type`

Character string indicating the data type. Controls validation messaging and DNAm-specific probe filtering.

- `"Proteomics"` (default)
- `"Metabolomics"`
- `"DNAm"` — triggers loading of probe lists from `Data/` and an additional filtered BH correction column in results
- `"other"` — accepts non-omics feature data without a preprocessing reminder or type-specific filtering

---

### `additional_covariates`

Optional character vector of column names in `pheno` to include as additional covariates in every model. Must be numeric, factor, or logical.

```r
additional_covariates = c("age", "bmi", "race")
```

---

### `n_cores` *(FAST_omics_WAS only)*

Integer. Number of cores to use for parallelizing the per-analyte model fits. Defaults to `NULL`, which auto-detects as `max(1, detectCores() - 1)`. Set to `1` to run serially.

Note: `detectCores()` may overcount in HPC/container environments — set explicitly if running on a cluster with allocated core limits.

---

### `checkpoint_dir` *(FAST_omics_WAS only)*

Optional character string. Path to a directory for saving per-batch checkpoint files. If `NULL` (default), checkpointing is disabled and behavior is unchanged. If provided, completed batches are written atomically to disk and reloaded on resume, so an interrupted run can continue from where it left off.

The directory structure created under `checkpoint_dir`:

```
checkpoint_dir/
  change/all/    batch_1.rds, batch_2.rds, ...
  change/male/   batch_1.rds, ...
  change/female/ batch_1.rds, ...
  level/all/     batch_1.rds, ...
  level/male/    batch_1.rds, ...
  level/female/  batch_1.rds, ...
```

**Important**: checkpoint files are tied to analyte ordering and batch size. Do not change the `omics` data or `checkpoint_batch_size` between a run and its resume.

---

### `checkpoint_batch_size` *(FAST_omics_WAS only)*

Integer. Number of analytes per checkpoint batch. Default: `2000`. Only relevant when `checkpoint_dir` is set. Larger batches mean fewer files and less I/O overhead; smaller batches mean less lost work if a crash occurs mid-batch.

---

## Outputs

### `FAST_omics_WAS()`

Returns a list with two top-level elements:

```r
list(
  analysis_change = ...,   # Change-score models (FU value − baseline value)
  analysis_level  = ...    # Level models (absolute FU value)
)
```

Each element is stratified the same way:

```r
$all     # Full dataset
$male    # Male subjects only  (NULL if dataset is single-sex)
$female  # Female subjects only (NULL if dataset is single-sex)
```

### `FAST_omics_WAS_reports()`

Returns a list of QC and data summary reports, stratified the same way:

```r
$all     # Full dataset
$male    # Male subjects only  (NULL if dataset is single-sex)
$female  # Female subjects only (NULL if dataset is single-sex)
```

---

## `FAST_omics_WAS()` — Analysis Outputs

### `$analysis_change` and `$analysis_level`

Each stratum contains two tables: `$coefficients` and `$treatment_effects`.

The `change` response is `FU_value − baseline_value`; the `level` response is
`FU_value`. Both include `analyte_baseline` as a covariate, so both are
baseline-adjusted. Interpretation of `TREATMENT_GROUP` effects differs:
`change` gives the treatment effect on the *change from baseline*, while
`level` gives the treatment effect on the *absolute post-treatment level*.

---

### `$coefficients`

All fixed-effect coefficients from the fitted model, for every analyte.
For a single follow-up (LM), the model is:

```
response ~ TREATMENT_GROUP + FEMALE + analyte_baseline + [additional covariates]
```

For multiple follow-ups (LME4):

```
response ~ TREATMENT_GROUP * FU + FEMALE + analyte_baseline + [additional covariates] + (1 | SUBJECT_ID)
```

One row per analyte × coefficient combination.

**Example (single follow-up, Proteomics):**

```
ANALYTE_NAME   COEFFICIENT          EFFECT_SIZE    SE       P_VALUE   BH_P_VALUE
protein_A      (Intercept)          0.231          0.089    0.010     0.051
protein_A      TREATMENT_GROUP1      0.312          0.067    0.000     0.002
protein_A      FEMALE1             -0.118          0.064    0.067     0.201
protein_A      analyte_baseline     0.794          0.041    0.000     0.000
protein_A      age                  0.004          0.003    0.189     0.389
protein_A      bmi                  0.007          0.005    0.211     0.389
protein_B      (Intercept)         -0.041          0.072    0.569     0.712
protein_B      TREATMENT_GROUP1     -0.028          0.061    0.647     0.712
protein_B      FEMALE1              0.053          0.058    0.362     0.589
protein_B      analyte_baseline     0.831          0.038    0.000     0.000
protein_B      age                 -0.002          0.003    0.501     0.712
protein_B      bmi                  0.003          0.004    0.488     0.712
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `ANALYTE_NAME` | Analyte identifier|
| `COEFFICIENT` | Name of the model term, e.g. `(Intercept)`, `TREATMENT_GROUP1`, `analyte_baseline`. |
| `EFFECT_SIZE` | Coefficient estimate (β) from the model. |
| `SE` | Standard error of the coefficient estimate |
| `P_VALUE` | Two-sided p-value. From `summary(lm_fit)$coefficients` for LM; from `summary(lmer_fit)$coefficients` with Satterthwaite degrees of freedom for LME4 |
| `BH_P_VALUE` | Benjamini-Hochberg corrected p-value, applied separately within each `COEFFICIENT` group across all analytes (e.g., all `TREATMENT_GROUP1` p-values are corrected together) |
| `BH_P_VALUE_FILTERED` | *(DNAm only)* BH correction applied only within the filtered probe set for probes in that set; `NA` for all other probes. See [DNAm Filtered BH](#dnam-filtered-bh) |

---

### `$treatment_effects`

The treatment-vs-control contrast at each follow-up level, one row per analyte × FU.
For LM (single FU), this is identical to the `TREATMENT_GROUP1` row in `$coefficients`.
For LME4 (multiple FUs), these are estimated marginal means contrasts from `emmeans`,
which properly account for the covariance between the `TREATMENT_GROUP` main effect
and the `TREATMENT_GROUP:FU` interaction.

**Example (two follow-up timepoints, LME4):**

```
ANALYTE_NAME   FU   EFFECT_SIZE    SE       P_VALUE   BH_P_VALUE
protein_A      1    0.289          0.071    0.000     0.002
protein_A      2    0.341          0.083    0.000     0.001
protein_B      1   -0.031          0.064    0.629     0.812
protein_B      2   -0.019          0.079    0.811     0.900
protein_C      1    0.105          0.068    0.123     0.287
protein_C      2    0.198          0.077    0.011     0.044
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `ANALYTE_NAME` | Identifier for the analyte |
| `FU` | Follow-up timepoint at which the contrast is evaluated |
| `EFFECT_SIZE` | Estimated treatment − control difference at that FU. For `change` models: difference in change from baseline. For `level` models: difference in absolute level |
| `SE` | Standard error of the contrast. For LME4, derived from `emmeans` and accounts for parameter covariance |
| `P_VALUE` | Two-sided p-value for the treatment vs. control contrast |
| `BH_P_VALUE` | Benjamini-Hochberg corrected p-value, applied separately within each `FU` level across all analytes |
| `BH_P_VALUE_FILTERED` | *(DNAm only)* As above; see [DNAm Filtered BH](#dnam-filtered-bh) |

---

## `FAST_omics_WAS_reports()` — Report Outputs

Each stratum (`$all`, `$male`, `$female`) contains four tables.

---

### `$pheno_summary`

One row per `(FU, FEMALE)` cell. Subject-level counts (`N_SUBJECTS`, `N_CONTROL`, `N_TREATMENT`) deduplicate by `SUBJECT_ID` so repeated samples for the same person at the same timepoint don't inflate counts. `N_SAMPLES` is the raw row count.

**Example:**

```
FU   FEMALE   N_SUBJECTS   N_CONTROL   N_TREATMENT   N_SAMPLES
0    0        85           42          43            85
0    1        93           46          47            93
1    0        83           41          42            83
1    1        91           45          46            91
2    0        81           40          41            81
2    1        89           44          45            89
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `FU` | Follow-up timepoint (0 = baseline) |
| `FEMALE` | Sex indicator: 0 = male, 1 = female |
| `N_SUBJECTS` | Number of unique subjects in this cell |
| `N_CONTROL` | Number of unique control subjects (`TREATMENT_GROUP == 0`) |
| `N_TREATMENT` | Number of unique treatment subjects (`TREATMENT_GROUP == 1`) |
| `N_SAMPLES` | Raw number of rows (samples) in this cell; equals `N_SUBJECTS` when each person has one sample per timepoint |

---

### `$omics_summary`

Per-analyte baseline distribution, computed at baseline, `FU == 0`, only. Describes the pre-treatment reference distribution for each feature.

**Example:**

```
ANALYTE_NAME   N_NONMISSING   MEAN    MEDIAN   SD      MIN     MAX
protein_A      178            1.204   1.189    0.342   0.445   2.187
protein_B      178            0.523   0.511    0.187   0.102   0.981
protein_C      175            2.017   2.003    0.401   1.102   3.214
protein_D      178            0.781   0.769    0.224   0.213   1.445
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `ANALYTE_NAME` | Identifier for the analyte |
| `N_NONMISSING` | Number of non-NA baseline samples for this analyte |
| `MEAN` | Mean value across baseline samples |
| `MEDIAN` | Median value across baseline samples |
| `SD` | Standard deviation across baseline samples |
| `MIN` | Minimum value across baseline samples |
| `MAX` | Maximum value across baseline samples |

---

### `$covariates_summary`

Only present when `additional_covariates` is specified; `NULL` otherwise. One row per covariate. The `SUMMARY` column is a **list column** — each cell contains a named list whose contents depend on the covariate type.

**Example:**

```
COVARIATE_NAME   TYPE      N_NA   SUMMARY
age              numeric   0      list(mean=55.2, median=55.0, sd=8.3, min=35.0, max=78.0)
bmi              numeric   3      list(mean=27.1, median=26.8, sd=4.2, min=18.5, max=41.2)
race             factor    0      list(n_levels=3, level_names=c("1","2","3"), counts=c(89,54,32))
smoker           logical   1      list(n_true=74, n_false=101)
```

To access a summary for a specific covariate:

```r
reports$all$covariates_summary$SUMMARY[[1]]  # by position
# or
idx <- which(reports$all$covariates_summary$COVARIATE_NAME == "age")
reports$all$covariates_summary$SUMMARY[[idx]]
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `COVARIATE_NAME` | Name of the covariate as provided in `additional_covariates` |
| `TYPE` | Detected R type: `"numeric"`, `"factor"`, or `"logical"` |
| `N_NA` | Number of missing values at baseline (`FU == 0`) |
| `SUMMARY` | List column. For **numeric**: `mean`, `median`, `sd`, `min`, `max`. For **factor**: `n_levels`, `level_names`, `counts` (one count per level). For **logical**: `n_true`, `n_false` |

---

### `$randomization_summary`

Baseline balance check: a Welch's t-test comparing treatment vs. control at `FU == 0` for each analyte. This is a diagnostic — in a well-randomized trial, no analyte should show a significant difference at baseline.

**Example:**

```
ANALYTE_NAME   MEAN_DIFFERENCE   COHENS_D   SE      P_VALUE   BH_P_VALUE
protein_A      0.012             0.035      0.041   0.770     0.912
protein_B     -0.023             0.123      0.044   0.601     0.912
protein_C      0.041             0.102      0.051   0.421     0.851
protein_D     -0.007             0.031      0.039   0.858     0.912
```

**Column descriptions:**

| Column | Description |
|:---|:---|
| `ANALYTE_NAME` | Identifier for the analyte |
| `MEAN_DIFFERENCE` | Mean(treatment) − Mean(control) at baseline |
| `COHENS_D` | Effect size: mean difference divided by pooled standard deviation |
| `SE` | Welch standard error: `sqrt(var_control/n_control + var_treatment/n_treatment)` |
| `P_VALUE` | Two-sided p-value from Welch's t-test (does not assume equal variances) |
| `BH_P_VALUE` | Benjamini-Hochberg corrected p-value across all analytes |

---

## DNAm Filtered BH

When `omics_type == "DNAm"`, two additional columns appear in `$coefficients` and `$treatment_effects` in both `analysis_change` and `analysis_level`:

| Column | Description |
|:---|:---|
| `BH_P_VALUE_FILTERED` | BH-corrected p-value computed only within the filtered probe set. `NA` for probes not in the filtered set |

This allows comparing significance under two different multiple-testing burdens — the full analyzed probe set (`BH_P_VALUE`) vs. a more restricted, pre-specified probe set (`BH_P_VALUE_FILTERED`) — without re-running the models. The filtered probe set is loaded from `Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds`.
