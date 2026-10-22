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
| `TREATMENT_GROUP` | factor (or numeric, auto-coerced) | `0` = control, `1` = treatment. Both groups must be present for a standard two-group study; assignment is constant within subject |
| `FEMALE` | factor (or numeric, auto-coerced) | Binary: `0` = male, `1` = female |
| *Additional covariates* | numeric / factor / logical | Any extra columns named in `additional_covariates`; NAs are allowed but reduce sample size |

If you wish to add additional intervention arms, use codes `2`, `3`, and so on
in `TREATMENT_GROUP` and supply [`treatment_labels`](#treatment_labels) for
every group, including control.

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

Model fitting uses matching sample IDs and reports the number of unmatched
samples. Missing baseline or follow-up measurements can further reduce the
observations usable for an individual analyte. Phenotype counts in the reports
describe retained phenotype rows, not each model's complete-case sample size.

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

### `treatment_labels`

Optional for the standard `0` = control, `1` = treatment study. Leaving it
as `NULL` uses the labels `"Control"` and `"Treatment"`. To supply custom names,
pass a named character vector to either public function:

```r
treatment_labels = c(`0` = "Usual care", `1` = "Diet")
```

If you wish to add additional intervention arms, give each its own code in
`TREATMENT_GROUP` (`2`, `3`, and so on). `treatment_labels` is then required
and must name every group, including control:

```r
treatment_labels = c(`0` = "Usual care", `1` = "Diet", `2` = "Diet + exercise")
```

When a map is supplied, its keys must match **exactly** the codes observed in the full
phenotype input, including control. Keys and labels must be unique and
nonblank, with no missing values. Validation happens before filtering or sex
stratification, so retain labels for arms later lost from a subset.

Codes are categorical: `0` is always control, positive integer codes identify
active arms, and the pipeline does not renumber them. Codes need not be
consecutive, and there is no fixed maximum number of arms. Neither the
caller’s phenotype nor omics data frame is modified. The normalized full map
is returned in `$metadata$treatment_labels`.

---

### `n_cores` *(FAST_omics_WAS only)*

Integer. Number of cores to use for parallelizing the per-analyte model fits.
Defaults to `NULL`, which auto-detects as `max(1, detectCores() - 1)` and falls
back to one worker if detection is unavailable. Set to `1` to run serially.

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

Resume with the same input data, model settings, and batch layout. Use a fresh
checkpoint directory when any of these change or when upgrading across this
result-schema change. You can change `n_cores` or display labels
in `treatment_labels` on resume; labels are applied after the batches are
assembled and do not require refitting.

---

### `checkpoint_batch_size` *(FAST_omics_WAS only)*

Integer. Number of analytes per checkpoint batch. Default: `2000`. Only relevant when `checkpoint_dir` is set. Larger batches mean fewer files and less I/O overhead; smaller batches mean less lost work if a crash occurs mid-batch.

---

### `verbose`

Logical. Default: `FALSE`. When `TRUE`, analysis and reporting emit `[1.1.1]` progress messages. Analysis messages include serial input-quality milestones every 100,000 analytes, distinguish serial preparation and correction from computed parallel batches, report configured and observed workers, and summarize cached versus pending checkpoints. Reporting messages identify serial stages and emit coarse per-analyte milestones.

---

## Outputs

### `FAST_omics_WAS()`

Returns:

```r
list(
  analysis_change = ...,   # Change-score models (FU value − baseline value)
  analysis_level  = ...,   # Level models (absolute FU value)
  metadata       = ...    # Schema version and treatment mapping
)
```

Each analysis element is stratified the same way:

```r
$all     # Full dataset
$male    # Male subjects only  (NULL if dataset is single-sex)
$female  # Female subjects only (NULL if dataset is single-sex)
```

### `FAST_omics_WAS_reports()`

Returns study-level counts and balance checks, plus stratified variable summaries:

```r
$pheno_summary
$arm_counts
$variable_summaries    # Contains $all, $male, $female
$randomization_reports
$metadata
```

---

## `FAST_omics_WAS()` — Analysis Outputs

### `$analysis_change` and `$analysis_level`

Each stratum contains two tables: `$coefficients` and `$treatment_effects`.

The `change` response is `FU_value − baseline_value`; the `level` response is
`FU_value`. Both include `analyte_baseline` as a covariate, so both are
baseline-adjusted. On the same data and model specification, these two
response parameterizations give the same treatment contrasts; subtracting
baseline changes the baseline coefficient rather than the treatment estimand.

Each model is a joint fit of all available arms in the stratum. Separate
arm-versus-control tables are views of those fitted contrasts, not separately
fitted two-arm models. Unequal allocation ratios are supported.

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

`TREATMENT_GROUP` in these formulas is a categorical arm factor with the
original input codes and control as reference. There is one coefficient for each active arm against
control, and arm-by-FU interactions in the repeated-follow-up model.

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

The active-arm-vs-control contrasts at each follow-up level, one row per
analyte × FU × contrast. For LM (single FU), each estimable contrast matches
its active-arm coefficient in `$coefficients` under the default treatment
coding. Other factor coding or aliased models use explicit contrasts instead.
For LME4 (multiple FUs), these are estimated marginal means contrasts from `emmeans`,
which properly account for the covariance between the `TREATMENT_GROUP` main effect
and the `TREATMENT_GROUP:FU` interaction.

Only active-minus-control contrasts are returned; A-minus-B is not added.
The returned effects share controls and are correlated. Their cross-contrast
covariance is not exported here; a formal A-minus-B test requires that
covariance, rather than treating the two standard errors as independent.
The raw `P_VALUE` uses no adjustment across contrasts (`adjust = "none"`).
Each contrast receives its own BH correction across analytes at each FU,
within each analysis and stratum. Nonestimable contrasts retain their expected
rows with `NA` estimates, uncertainty, and p-values.

```r
te <- results$analysis_change$all$treatment_effects
effects_by_treatment <- split(te, te$TREATMENT)
effects_diet <- effects_by_treatment[["Diet"]]
effects_diet_exercise <- effects_by_treatment[["Diet + exercise"]]
```

**Example (binary analysis, two follow-ups; identity/status columns omitted):**

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
| `CONTRAST` | Stable comparison identifier built from input codes, e.g. `1_vs_0` or `2_vs_0` |
| `TREATMENT` | Active-arm display label; `"Treatment"` for the default binary call |
| `REFERENCE` | Control display label; `"Control"` for the default binary call |
| `EFFECT_SIZE` | Estimated active-arm − control difference at that FU, adjusted for baseline and model covariates |
| `SE` | Standard error of the contrast. For LME4, derived from `emmeans` and accounts for parameter covariance |
| `P_VALUE` | Two-sided p-value for the treatment vs. control contrast |
| `BH_P_VALUE` | Benjamini-Hochberg corrected p-value across analytes within each `FU` and `CONTRAST`, separately for each analysis and stratum |
| `BH_P_VALUE_FILTERED` | *(DNAm only)* As above; see [DNAm Filtered BH](#dnam-filtered-bh) |
| `IN_FILTERED_SET` | *(DNAm treatment effects only)* Logical `TRUE`/`FALSE`: `ANALYTE_NAME %in% filtered_probes`, regardless of `STATUS` or missing inference |
| `STATUS` | `ok`, `not_estimable`, or `fit_error`; the latter two retain rows with missing inference |

### `$metadata`

Both analysis and report objects include result schema version `"2.0"`, the
normalized full `treatment_labels` map, and a `treatments` data frame containing
`ARM` (original active code), `TREATMENT` (active label), `REFERENCE` (control
label), and `CONTRAST` (code-based identifier). Metadata also records the BH
correction policy. Preserve the mapping alongside any exported table.

---

## `FAST_omics_WAS_reports()` — Report Outputs

### `$pheno_summary`

Study-level table with one row per `(FU, FEMALE)` cell. `N_SUBJECTS`,
`N_CONTROL`, and `N_TREATMENT` count unique subjects; `N_SAMPLES` counts rows.
`N_TREATMENT` includes all active arms. Use `$arm_counts` for individual arms.

### `$arm_counts`

One row per FU × FEMALE × arm, including zero-count cells.

| Column | Description |
|:---|:---|
| `FU`, `FEMALE` | Visit and sex cell |
| `ARM` | Original input code: `0` control, positive codes for active arms |
| `TREATMENT` | Display label from `treatment_labels`, including the control label |
| `IS_CONTROL` | Whether this is the shared control arm |
| `N_SUBJECTS` | Unique subjects in this cell |
| `N_SAMPLES` | Sample rows in this cell |

The control count is reported once per cell, not repeated as a new set of
participants for each comparison.

### `$variable_summaries`

Contains `$all`, `$male`, and `$female` strata where applicable. Within each
stratum, summaries are separated by FU and original arm code:

```r
reports$variable_summaries$all$omics_FU0_Tx0       # Control, baseline
reports$variable_summaries$all$omics_FU1_Tx1       # First active arm, FU1
reports$variable_summaries$all$omics_FU1_Tx2       # Second active arm, FU1
reports$variable_summaries$all$covariates_FU0_Tx0
reports$metadata$treatment_labels                 # Labels for every Tx code
```

Empty cells are omitted. Omics entries have one row per analyte with
`ANALYTE_NAME`, `N_NONMISSING`, `MEAN`, `MEDIAN`, `SD`, `MIN`, and `MAX`.
Covariate entries are present when additional covariates were requested and
have `COVARIATE_NAME`, `TYPE`, `N_NA`, and a `SUMMARY` list column. Numeric
summaries contain descriptive statistics; factor summaries contain levels
and counts; logical summaries contain true/false counts.

### `$randomization_reports`

Study-level baseline balance checks, separately for each active arm versus
shared control. Each table includes `CONTRAST`, `TREATMENT`, and `REFERENCE`.

- `$omics_report`: one row per analyte and contrast. Welch's t-test with
  `MEAN_DIFFERENCE`, `COHENS_D`, `SE`, `P_VALUE`, and `BH_P_VALUE`; BH correction
  is across analytes within each contrast.
- `$covariate_report`: `FEMALE` when both sexes are present, plus requested
  additional covariates. Numeric variables use Welch's t-test; logical
  variables and factors use chi-squared or Fisher's exact test
  according to cell counts. Columns include `VARIABLE`, `TYPE`, `TEST`,
  `STATISTIC`, `MIN_CELL_COUNT`, `P_VALUE`, `SUMMARY_CONTROL`, and
  `SUMMARY_TREATMENT`.

These are baseline diagnostics, not treatment-effect estimates. Omics rows
with fewer than two available values in either arm are omitted. If a baseline
arm is absent, covariate summaries remain available but inferential values
are `NA`.

---

## DNAm Filtered BH

When `omics_type == "DNAm"`, the following columns appear in both
`analysis_change` and `analysis_level`:

| Column | Tables | Description |
|:---|:---|:---|
| `BH_P_VALUE_FILTERED` | `$coefficients` and `$treatment_effects` | BH-corrected p-value computed only within the filtered probe set. `NA` for probes outside the set or with missing p-values |
| `IN_FILTERED_SET` | `$treatment_effects` only | Logical `TRUE`/`FALSE`, equal to `ANALYTE_NAME %in% filtered_probes` |

The grouping matches the full-set correction: coefficient-wise for
`coefficients`, and FU-by-contrast for `treatment_effects`, separately within
each analysis and stratum.

`IN_FILTERED_SET` explicitly records membership, which cannot be inferred from
`BH_P_VALUE_FILTERED = NA` when inference is missing. It is independent of
`STATUS` and is the same for a probe across all follow-ups, contrasts, analyses,
and strata. Membership is attached after batch assembly without changing model
fits or BH correction. It is absent from coefficient tables and non-DNAm results.
DNAm plots use `IN_FILTERED_SET` for filtered membership and show labeled empty
plots when a cell has no available inference.

This allows comparing significance under two different multiple-testing burdens — the full analyzed probe set (`BH_P_VALUE`) vs. a more restricted, pre-specified probe set (`BH_P_VALUE_FILTERED`) — without re-running the models. The filtered probe set is loaded from `Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds`.


## Existing downstream consumers

This result schema adds multi-arm support to Track 1.1.1 only. Summary consumers
that support one comparison must receive one selected `TREATMENT`/`CONTRAST`,
with its mapping retained as provenance. Track 3.3 additionally needs its outcome
model to adjust for the full arm factor; splitting these tables alone does not
complete the mediation extension. Neither consumer has been upgraded here to
accept combined multi-arm tables.
