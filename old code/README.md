# Overview

## Repository Structure

| File | Purpose |
|---|---|
| `README.md` | Provides an overview of the repository and explains how to reproduce the results in the paper. |
| `gp_test.R` | Contains the main functions for the generalized projection test. |
| `example_1_new.R` | Contains the simulation code for Example 1. |
| `example_2_new.R` | Contains the simulation code for Example 2. |
| `real_data.R` | Contains the code used to replicate the real data analysis results. The original dataset is not included if it is not publicly available. |


# Generalized Projection Test functions in R

This file `gp_test.R` contains R functions for constructing orthonormal Legendre basis functions and implementing a generalized projection test. The code also includes helper functions for summarizing rejection behavior using different reference approximations.

## Overview

The main functions in this code are:

- `legendre_orthonormal_matrix()`: constructs orthonormal Legendre basis functions
- `gp_test()`: computes the Wald statistic and generalized projection test statistics
- `test_summary_wald()`: summarizes rejection based on the chi-square approximation for the Wald test
- `test_summary_series()`: summarizes rejection based on the normal approximation for the generalized projection test
- `test_summary_series2()`: summarizes rejection based on a simulation-based chi-square approximation

This code is useful for simulation studies and empirical applications where the user wants to test whether a pseudo-outcome is associated with a set of covariates using either a standard linear projection approach or a more flexible series-based approach.

---

## Functions

### 1. `legendre_orthonormal_matrix(x, k)`

Constructs the orthonormal Legendre basis evaluated at a numeric vector `x` up to degree `k`.

#### Arguments

- `x`: a numeric vector
- `k`: a nonnegative integer giving the maximum degree of the basis

#### Returns

A matrix with `length(x)` rows and `k + 1` columns.  
The columns correspond to the orthonormal Legendre basis functions:
- `phi0`
- `phi1`
- ...
- `phik`

#### Notes

- The function first builds the standard Legendre polynomials using a recurrence relation.
- It then rescales them so that the columns form an orthonormal basis on `[-1,1]`.
- In practice, `x` should usually be scaled to lie in `[-1,1]` before using this function.

---

### 2. `gp_test(psedo_outcome, covariates, k_vec = c(2,3,5))`

Computes the generalized projection test statistics. Note that here the default setting is to create basis function for each covariate without considering their interactions.

#### Arguments

- `psedo_outcome`: a numeric vector containing the pseudo-outcome
- `covariates`: a matrix or data frame of covariates
- `k_vec`: a vector of user-specified basis dimensions for the series test

#### Returns

A list with the following elements:

- `stat_Wald`: the Wald test statistic from linear regression of the pseudo-outcome on the covariates
- `stat_series`: a vector of standardized generalized projection test statistics
- `S_vec`: a vector of unstandardized quadratic-form statistics
- `M_list`: a list of covariance-related matrices used for the chi-square approximation

