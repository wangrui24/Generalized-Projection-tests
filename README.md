# Generalized projection tests

R code accompanying **Generalized projection tests for function-valued parameters with applications to testing structural causal assumptions**.

The two main simulation scripts, [`cluster_simu.R`](cluster_simu.R) and [`cluster_simu_high.R`](cluster_simu_high.R), implement the simulation studies for testing mean exchangeability in data fusion and the instrumental-variable (IV) compatibility condition. Both use the shared data-generating processes, nuisance estimation routines, and test statistics in [`functions.R`](functions.R).

## Repository structure

```text
Generalized-Projection-tests/
├── README.md
├── functions.R             # Shared simulation and testing functions
├── cluster_simu.R          # Simulations with two continuous covariates
├── cluster_simu_high.R     # Simulations with ten covariates
└── old code/               # Earlier implementations retained for reference
    ├── README.md
    ├── gp_test.R
    ├── example_1_new.R
    ├── example_2_new.R
    └── real_data.R
```

The root-level scripts are the entry points for the simulation tables described below. The `old code/` directory contains earlier simulation and real-data code; its README describes that earlier implementation. The current real-data analysis scripts and datasets are not included in this folder.

## Correspondence with the paper

Table numbers below refer to the manuscript revision supplied with this repository update (revision 33).

| Paper table | Simulation script | Covariates | GP `basis_type` | GP `test_type` |
| --- | --- | --- | --- | --- |
| Table 1 | `cluster_simu.R` | Two continuous | `fourier` | `standardized` |
| Table S2 | `cluster_simu.R` | Two continuous | `legendre` | `standardized` |
| Table S3 | `cluster_simu.R` | Two continuous | `fourier` | `unstandardized` |
| Table S4 | `cluster_simu.R` | Two continuous | `legendre` | `unstandardized` |
| Table S5 | `cluster_simu_high.R` | Five continuous and five binary | `fourier` | `standardized` |
| Table S6 | `cluster_simu_high.R` | Five continuous and five binary | `legendre` | `standardized` |
| Table S7 | `cluster_simu_high.R` | Five continuous and five binary | `fourier` | `unstandardized` |
| Table S8 | `cluster_simu_high.R` | Five continuous and five binary | `legendre` | `unstandardized` |
| Table S9 | Both scripts | Usual and higher-dimensional settings | Method timing: `method_elapsed_seconds` | See timing notes below |
| Table S10 | Both scripts | Usual and higher-dimensional settings | Nuisance timing: `nuisance_elapsed_seconds` | See timing notes below |

For Tables 1 and S2–S8:

- **Panel A** is Example 1: testing mean exchangeability for the control potential outcome across two data sources (Section 5.2; higher-dimensional extension in Supplement S8.2). Select `example == "example_1"`.
- **Panel B** is Example 2: testing the IV compatibility condition by comparing conditional Wald estimands from two instruments (Section 5.3; higher-dimensional extension in Supplement S8.2). Select `example == "example_2"`.
- **LP test** rows have `method == "projection"` and `test_type == "wald"`.
- **GP test** rows have `method == "gp"`; use the basis and test-type filters in the table above. A fixed-dimension result has `combined == FALSE`; the **GP test (Combined)** result has `combined == TRUE` and uses a Bonferroni-adjusted minimum p-value across the candidate dimensions, separately for each basis and test type.
- **Racine** rows have `method == "racine"` and are produced only for Example 1. The LP and Racine results are shared across the corresponding basis/test-type tables and have no GP basis label.

Each script generates replicate-level CSV files. Rejection proportions must be aggregated across replicates to form the paper's tables; the scripts do not directly produce formatted manuscript tables. This mapping identifies the simulation designs and output selections, rather than asserting that a new run has been verified to reproduce every published numerical entry.

### Scenarios and basis dimensions

The five scenarios have the same ordering in both scripts. `parameter_1` and `parameter_2` contain the two components of `alpha` for Example 1 and `beta` for Example 2.

| `scenario_index` | Paper scenario | Interpretation | Example 1: `(alpha_1, alpha_2)` | Example 2: `(beta_1, beta_2)` |
| --- | --- | --- | --- | --- |
| 1 | I | Null hypothesis | `(0, 0)` | `(0, 0)` |
| 2 | II | Weak nonlinear misalignment | `(0.2, 0)` | `(0.3, 0)` |
| 3 | III | Weak linear misalignment | `(0, 0.2)` | `(0, 0.3)` |
| 4 | IV | Strong, predominantly nonlinear misalignment | `(0.4, 0.2)` | `(0.6, 0.3)` |
| 5 | V | Strong, predominantly linear misalignment | `(0.2, 0.4)` | `(0.3, 0.6)` |

The active sample-size and dimension grids in both scripts are:

| Example | Sample size | Candidate `J` values (paper's J*) |
| --- | --- | --- |
| 1 | 250, 500 | 4, 8, 16 |
| 1 | 1000, 1500 | 5, 10, 20 |
| 2 | 1000, 2000 | 5, 10, 20 |
| 2 | 3000, 5000 | 5, 10, 20, 40 |

`J` is the number of nonconstant basis functions **per continuous covariate**, not the total number of columns. The basis is additive across covariates, with one shared intercept and no interaction terms. The usual setting has `1 + 2 * J` columns. The higher-dimensional setting expands `x1`–`x5` and includes `x6`–`x10` once each as binary terms, giving `1 + 5 * J + 5` columns.

## Running the simulations

Install R and the packages used by the drivers and their SuperLearner libraries:

```r
install.packages(c("SuperLearner", "caret", "np", "randomForest", "gam"))
```

Each invocation runs **one replicate across both examples, all sample sizes, and all five scenarios**. The replicate identifier comes from `SLURM_ARRAY_TASK_ID`, and the random seed is `2000 + SLURM_ARRAY_TASK_ID`. Set this variable even when running locally.

From the repository directory, run one full replicate of each design:

```bash
SLURM_ARRAY_TASK_ID=1 Rscript cluster_simu.R
SLURM_ARRAY_TASK_ID=1 Rscript cluster_simu_high.R
```

For a Slurm cluster, submit separate arrays after configuring R and any cluster-specific resource options:

```bash
sbatch --array=1-1000 --wrap="Rscript cluster_simu.R"
sbatch --array=1-1000 --wrap="Rscript cluster_simu_high.R"
```

The paper uses 1,000 simulated datasets per configuration in the usual setting. The timing-table notes report 999 retained higher-dimensional replicates after one SuperLearner fitting failure. Track completed replicates and fitting failures when aggregating a new run. Exact results can depend on package versions and the random-number environment, and elapsed times depend on hardware.

The scripts locate `functions.R` relative to their own location and create these output folders there:

| Script | Example 1 output | Example 2 output |
| --- | --- | --- |
| `cluster_simu.R` | `output1/sim_<id>.csv` | `output2/sim_<id>.csv` |
| `cluster_simu_high.R` | `output1_high/sim_<id>.csv` | `output2_high/sim_<id>.csv` |

Reusing an identifier overwrites its CSV files. Both examples are computed before the final CSV-writing calls, so an uncaught fitting failure can leave the invocation without output files.

### Reduced test run

The higher-dimensional driver provides a reduced mode:

```bash
CLUSTER_SIMU_TEST_MODE=1 SLURM_ARRAY_TASK_ID=1 Rscript cluster_simu_high.R
```

This uses only the null scenario, sample sizes 120 and 160, `J = 2`, GLM-only nuisance libraries, and 100 weighted chi-square draws. It writes to `output1_high_test/` and `output2_high_test/`. This mode checks the workflow; it does not reproduce the paper's settings. The usual driver has no corresponding test-mode switch.

## Reading and summarizing the outputs

Each row records the replicate, example, sample size, scenario, parameter values, method, basis, test type, dimension, combination indicator, statistic, p-value, timing, and any caught Racine-test error. Rejection indicators are provided at levels 0, 0.05, ..., 1; for example, `reject_alpha_0_05` records rejection at the 5% level.

For example, after running the usual-setting simulations, this R code summarizes the Table 1 Panel A GP results at 5%:

```r
files <- list.files("output1", pattern = "^sim_[0-9]+\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0L)
results <- do.call(rbind, lapply(files, read.csv))

selected <- subset(results,
  method == "gp" & basis_type == "fourier" & test_type == "standardized"
)
stopifnot(all(is.finite(selected$p_value)))
table1_panel_a_gp <- aggregate(
  reject_alpha_0_05 ~ sample_size + scenario_index + J + combined,
  data = selected,
  FUN = mean
)
table1_panel_a_gp
```

Select `method == "projection"` or `method == "racine"` separately for the comparison methods. Check missing p-values and `error_message` before averaging: the driver's rejection indicator is zero when a p-value is missing, so blindly averaging failed rows would count them as non-rejections. Also check that each configuration has the intended number of distinct replicates.

For **Table S9**, count each `method_elapsed_seconds` value once per replicate, example, sample size, scenario, method, and basis. A GP block covers basis construction, all candidate dimensions, both test types, 10,000 weighted chi-square draws per dimension, and the combined tests. Its time is repeated across the resulting GP rows and must not be summed across them. For **Table S10**, count nuisance time once per replicate, example, sample size, and scenario. The paper pools the five scenarios within each example/sample-size/design group and reports mean, median, and the 25th and 75th percentiles.

## Using `gp_test()` with your own data

`gp_test()` is the main reusable testing function. It takes an already constructed pseudo-outcome and user-supplied basis evaluations. The target conditional moment is **E[g(O) | X] = 0**; the GP statistics measure projections of this moment onto the supplied basis. Power depends on the alternatives represented by that basis. An additive basis, for example, may miss alternatives that require interaction terms.

The function computes statistics only. It does not fit nuisance models, cross-fit observations, construct a basis automatically, or return p-values. For applications with estimated nuisance functions, construct the appropriate orthogonalized, cross-fitted pseudo-outcome first, following the assumptions in the paper. An arbitrary response vector does not automatically give a valid test of a structural causal assumption.

### Arguments and accepted inputs

```r
gp_test(
  psedo_outcome,
  transformed_covariates,
  wald_covariates = NULL,
  wald_include_intercept = TRUE,
  basis_degrees = NULL
)
```

| Argument | How to supply it |
| --- | --- |
| `psedo_outcome` | A finite numeric vector of length `n`, with one pseudo-outcome per observation. The spelling `psedo_outcome` is the actual argument name; `pseudo_outcome` is not an accepted named argument. |
| `transformed_covariates` | An `n`-row numeric matrix/data frame for one GP test; a nonempty list of such matrices for multiple tests; or a structured list with `series_covariates` and optional `wald_covariates` and `basis_degrees` fields. Each matrix may have a different number of columns. |
| `wald_covariates` | Optional `n`-row numeric covariates for the Wald comparison. An explicitly supplied value overrides the structured input's `wald_covariates`. If neither is supplied, the first GP basis matrix is used directly as the Wald design. |
| `wald_include_intercept` | Adds a column of ones to the supplied Wald covariates when `TRUE`. Use `FALSE` if those covariates already contain an intercept. This argument does not change the GP bases, and has no effect when the Wald design falls back to the first GP basis. |
| `basis_degrees` | Optional labels, one per candidate GP basis, returned with the results. Explicit labels override those in a structured input. Labels do not construct, truncate, or otherwise modify a basis matrix. |

All rows must refer to the same observations in exactly the same order. Supply complete, finite numeric inputs; encode categorical variables before calling the function. `gp_test()` does not impute missing values or create factor contrasts.

Include any desired constant column in each GP basis yourself. The Legendre helper below includes one shared constant by default. `gp_test()` does not orthogonalize, center, or rescale supplied bases: their normalization affects the GP statistics, so use a basis and scaling consistent with the intended test in the paper.

### Example 1: a complete Legendre-based call

This self-contained example uses only base R and `functions.R`. It illustrates a conditional-mean test with a directly observed signal; it does not require fitting a causal model or installing the simulation packages.

```r
source("functions.R")
set.seed(2026)

n <- 600
x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
g <- 0.4 * cos(pi * x$x1) + rnorm(n)
# For a null example, replace the preceding line with: g <- rnorm(n)

J_values <- c(2L, 4L, 8L)
bases <- legendre_gp_transformed_covariates(
  covariates = x,
  k_vec = J_values,
  include_default = FALSE
)
fit <- gp_test(psedo_outcome = g, transformed_covariates = bases)

# One upper-tail normal p-value per standardized GP statistic.
p_standardized <- pnorm(fit$stat_series, lower.tail = FALSE)
gp_results <- data.frame(
  J = fit$basis_degrees,
  statistic = fit$stat_series,
  p_value = p_standardized,
  reject_at_0_05 = p_standardized <= 0.05
)
gp_results

# The helper supplies x as the Wald covariates; an intercept is added.
p_wald <- pchisq(as.numeric(fit$stat_Wald),
                 df = fit$wald_df, lower.tail = FALSE)
p_wald

# Bonferroni combination across the three candidate dimensions.
p_combined_standardized <- min(1, length(p_standardized) * min(p_standardized))
p_combined_standardized
```

`include_default = FALSE` ensures that exactly the requested dimensions are tested. Otherwise, the helper prepends `ceiling(n^(1/3) / log10(n))`, even when it duplicates an entry in `k_vec`. The continuous covariates here already lie in `[-1, 1]`, the domain used for the Legendre basis; for other domains, choose and document an appropriate transformation before constructing the basis.

Large positive `stat_series` values provide evidence against the projected null. These normal p-values use the paper's standardized asymptotic calibration; they are not exact finite-sample probabilities. The Wald comparison tests the coefficients for the entire supplied Wald design, including its intercept.

### Example 2: supply one basis or a list of bases

Continuing with `x`, `g`, and `J_values` above, the following calls illustrate the other accepted input formats:

```r
# A single matrix: one GP statistic.
B <- legendre_transformed_covariates(x, k = 4)
fit_one <- gp_test(
  psedo_outcome = g,
  transformed_covariates = B,
  wald_covariates = x,
  basis_degrees = 4L
)

# A list of matrices: one GP statistic for each element, in list order.
B_list <- lapply(J_values, function(J) legendre_transformed_covariates(x, k = J))
fit_list <- gp_test(
  psedo_outcome = g,
  transformed_covariates = B_list,
  wald_covariates = x,
  basis_degrees = J_values
)

# The equivalent structured input needs no special class attribute.
fit_structured <- gp_test(
  psedo_outcome = g,
  transformed_covariates = list(
    series_covariates = B_list,
    wald_covariates = x,
    basis_degrees = J_values
  )
)
```

You can replace `B` or the elements of `B_list` with your own numeric basis evaluations, including appropriately constructed Fourier or interaction bases. Source only `functions.R` when using the reusable function: sourcing a cluster driver also starts its simulation run.

Always specify `wald_covariates = x` when you want the LP comparison on the original covariates. If you omit it with a plain matrix/list input, the Wald statistic tests coefficients on the **first supplied series basis**, which is a different comparison. The function always computes this Wald statistic before the GP statistics; its design and sandwich covariance must be nonsingular even if you only want the GP output. Avoid duplicate constant columns, redundant predictors, and designs with too many columns relative to the sample size. A zero covariance matrix also makes GP standardization undefined.

### Example 3: unstandardized p-values

`fit$S_vec[j]` is calibrated against a weighted sum of independent chi-square variables with one degree of freedom; the weights are the eigenvalues of `fit$M_list[[j]]`. Continuing from Example 1:

```r
set.seed(2027)
n_draws <- 10000L
p_unstandardized <- vapply(seq_along(fit$S_vec), function(j) {
  weights <- eigen(fit$M_list[[j]], symmetric = TRUE, only.values = TRUE)$values
  # M is positive semidefinite; clip numerical negative eigenvalues to zero.
  weights <- pmax(weights, 0)
  draws <- matrix(rchisq(n_draws * length(weights), df = 1), nrow = n_draws)
  null_statistics <- as.numeric(draws %*% weights)
  mean(null_statistics >= fit$S_vec[j])
}, numeric(1))

data.frame(J = fit$basis_degrees, statistic = fit$S_vec,
           p_value = p_unstandardized)
p_combined_unstandardized <- min(1,
  length(p_unstandardized) * min(p_unstandardized)
)
p_combined_unstandardized
```

These are Monte Carlo estimates of tail probabilities under the estimated reference distribution, with resolution `1 / n_draws`; a reported zero means no simulated reference statistic exceeded the observed statistic. The example matches the drivers' use of 10,000 draws. Combine dimensions separately for each chosen basis family and calibration, as in the simulation scripts.

### Connecting a fitted pseudo-outcome to the test

For the repository's causal examples, replace the toy construction of `g` with the relevant estimator, then use its returned covariates:

```r
# Requires the nuisance-fitting packages and an appropriately formatted dataset.
# Example 1 expects a, s, y, x1, x2 by default.
# pseudo <- estimate_mean_exchangeability_pseudo_outcome(dataset, V = 2)
# Example 2 instead expects d, z1, z2, y, x1, x2 by default.
# pseudo <- estimate_iv_compatibility_pseudo_outcome(dataset, V = 2)

# bases <- legendre_gp_transformed_covariates(
#   pseudo$x_vec, k_vec = c(2, 4, 8), include_default = FALSE
# )
# fit <- gp_test(psedo_outcome = pseudo$psedo_outcome,
#                transformed_covariates = bases)
```

These estimators concatenate held-out folds, so their returned `x_vec` must be paired with their returned pseudo-outcome. Using the original dataset's covariates without restoring matching row order can silently test the wrong pairings. Basis choice, nuisance estimation, and calibration should follow the assumptions and design appropriate to your application; the toy dimensions above are illustrative.

## Structure of `functions.R`

Sourcing `functions.R` defines reusable functions without starting simulations. Its organization follows the analysis workflow:

```text
Generate a dataset
  -> cross-fit nuisance models and construct a pseudo-outcome
  -> build candidate basis matrices
  -> compute LP and GP statistics
  -> calibrate p-values, combine dimensions, and save rows in the driver
```

### 1. Utilities, basis construction, and model fitting

- `require_package()`, `expit()`, and `logit()` provide package checking and link-function utilities.
- `legendre_orthonormal_matrix(x, k)` constructs degrees 0 through `k` of the orthonormal Legendre basis on `[-1, 1]`.
- `legendre_transformed_covariates()` builds an additive basis across covariates, retaining only one intercept.
- `legendre_gp_transformed_covariates()` packages multiple basis matrices, the original covariates for the Wald test, and dimension labels in a `gp_transformed_covariates` object.
- `fit_super_learner()`, `predict_super_learner()`, and `make_cv_folds()` support nuisance fitting, prediction, and sample splitting. The simulation drivers use two outer cross-fitting folds; SuperLearner uses five internal cross-validation folds by default.

Fourier basis construction is defined in the simulation drivers. The higher-dimensional driver also defines `mixed_transformed_covariates()` to handle continuous and binary covariates together.

### 2. Example 1: mean exchangeability

- `dgp_mean_exchangeability()` and `dgp_mean_exchangeability_high_dim()` generate the usual and higher-dimensional data-fusion settings.
- `estimate_mean_exchangeability_pseudo_outcome()` fits outcome, source-membership, and treatment-assignment models and constructs the cross-fitted, augmented inverse-probability-weighted difference between sources for a specified treatment level (default: control, `treatment = 0`).
- `run_mean_exchangeability_gp_test()` is a convenience wrapper for generating one usual-setting dataset and applying Legendre-based tests.
- `mean_exchangeability_alpha_list()` defines the five scenarios. `mean_exchangeability_sample_size_list()` is a separate convenience grid; the replication drivers use their own `example_1_j_map`.

The default outcome library contains `SL.glm` and `SL.randomForest`; source and treatment probabilities use `SL.glm`.

### 3. Example 2: IV compatibility

- `dgp_iv_compatibility()` and `dgp_iv_compatibility_high_dim()` generate two instruments, latent principal strata, treatment receipt, and potential outcomes under the five compatibility scenarios.
- `estimate_iv_compatibility_pseudo_outcome()` fits instrument probabilities and instrument-specific treatment and outcome regressions. It constructs an orthogonalized pseudo-outcome for each conditional Wald estimand and returns their difference.
- `run_iv_compatibility_gp_test()` provides a usual-setting, Legendre-based convenience wrapper.
- `iv_compatibility_beta_list()` defines the five scenarios. `iv_compatibility_sample_size_list()` is a separate convenience grid; the replication drivers use `example_2_j_map`.

Instrument and treatment models use `SL.glm` and `SL.randomForest`; outcome models additionally use `SL.gam`.

Both pseudo-outcome estimators return `pseudo_outcome`, the backward-compatible alias `psedo_outcome`, and the aligned covariates `x_vec`. The returned rows follow held-out fold order, so use `x_vec` from the same returned object when testing. The convenience wrappers have their own dimension defaults; use the cluster drivers for the paper's grids.

### 4. Test statistics and calibration helpers

`as_numeric_matrix()` checks matrix dimensions. `normalize_gp_transformed_covariates()` accepts a single basis matrix, a list of matrices, or a structured `gp_transformed_covariates` object. It standardizes the input format, not the numerical scale of the covariates.

`gp_test(psedo_outcome, transformed_covariates, ...)` computes the linear-projection Wald statistic with a residual-based sandwich covariance and the GP statistics for each supplied basis. For pseudo-outcome `g` and basis matrix `B`, it computes:

```text
projection = B' g / n
S          = n * sum(projection^2)
M          = B' diag(g^2) B / n
T          = (S - trace(M)) / (sqrt(2) * ||M||_F)
```

| Returned component | Meaning |
| --- | --- |
| `stat_Wald` | Linear-projection Wald statistic |
| `stat_series` | Standardized GP statistics `T`, one per candidate basis |
| `S_vec` | Unstandardized quadratic statistics `S` |
| `M_list` | Matrices whose eigenvalues determine the weighted chi-square reference distributions |
| `wald_df` | Number of coefficients in the Wald test |
| `basis_degrees` | Candidate dimension labels |

The drivers turn standardized statistics into upper-tail normal p-values and simulate weighted chi-square reference distributions for the unstandardized statistics. They also calculate Bonferroni combinations and rejection indicators.

`test_summary_wald()`, `test_summary_series()`, and `test_summary_series2()` are legacy helpers returning rejection indicators over a grid of reference-distribution quantiles. Their grid represents quantile probabilities; the driver output `reject_alpha_*` instead uses significance levels directly.

Finally, `test_conditional_independence_np_example_1()` implements the Racine comparison among control participants using `np`, and `extract_np_significance_component()` extracts the source-membership test by variable name. The higher-dimensional Racine helper is defined in `cluster_simu_high.R`. Aliases `dat_gen_1`, `test_a0_1`, `dat_gen`, and `nptest` preserve compatibility with earlier scripts.
