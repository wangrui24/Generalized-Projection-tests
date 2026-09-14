touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(touse+2000)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
setwd(script_dir)

source("functions.R")

suppressPackageStartupMessages({
  library(SuperLearner)
  library(caret)
  library(np)
})

dir.create("output1", showWarnings = FALSE, recursive = TRUE)
dir.create("output2", showWarnings = FALSE, recursive = TRUE)

alpha_grid <- seq(0, 1, 0.05)
chi_square_simu_time <- 10000

time_expression <- function(expr) {
  start_time <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(
    value = value,
    elapsed_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  )
}

# example_1_j_map <- list(
#   `250` = c(3, 6, 12),
#   `500` = c(3, 6, 12),
#   `1000` = c(3, 6, 12, 24),
#   `1500` = c(4, 8, 16)
# )
# 
# example_2_j_map <- list(
#   `1000` = c(3, 6, 12, 24),
#   `2000` = c(4, 8 ,16, 32),
#   `3000` = c(4, 8 ,16, 32),
#   `5000` = c(4, 8 ,16, 32)
# )
# 

example_1_j_map <- list(
  `250` = c(4, 8, 16),
  `500` = c(4, 8, 16),
  `1000` = c(5, 10, 20),
  `1500` = c(5, 10, 20)
)

example_2_j_map <- list(
  `1000` = c(5, 10, 20),
  `2000` = c(5, 10, 20),
  `3000` = c(5, 10, 20,40),
  `5000` = c(5, 10, 20,40)
)

fourier_orthonormal_matrix <- function(x, J) {
  stopifnot(is.numeric(x), J >= 0, J == as.integer(J))

  J <- as.integer(J)
  out <- matrix(0, nrow = length(x), ncol = J)
  if (J == 0) {
    return(out)
  }

  column_names <- character(J)
  for (m in seq_len(J)) {
    frequency <- (m + 1L) %/% 2L
    if (m %% 2L == 1L) {
      out[, m] <- cos(frequency * pi * x)
      column_names[m] <- paste0("cos", frequency)
    } else {
      out[, m] <- sin(frequency * pi * x)
      column_names[m] <- paste0("sin", frequency)
    }
  }

  colnames(out) <- column_names
  out
}

fourier_transformed_covariates <- function(
    covariates,
    J,
    include_intercept = TRUE,
    intercept_value = 1 / sqrt(2)) {
  x_vec <- as.matrix(covariates)
  fourier_basis <- NULL

  for (j in seq_len(ncol(x_vec))) {
    fj_basis <- fourier_orthonormal_matrix(x_vec[, j], J)
    colnames(fj_basis) <- paste0("x", j, "_fourier_", colnames(fj_basis))
    fourier_basis <- cbind(fourier_basis, fj_basis)
  }

  if (include_intercept) {
    fourier_basis <- cbind(rep(intercept_value, nrow(x_vec)), fourier_basis)
    colnames(fourier_basis)[1] <- ""
  }

  fourier_basis
}

make_transformed_covariates <- function(covariates, J_values, basis_type) {
  if (basis_type == "legendre") {
    basis_list <- lapply(
      J_values,
      function(J) legendre_transformed_covariates(covariates, k = J)
    )
  } else if (basis_type == "fourier") {
    basis_list <- lapply(
      J_values,
      function(J) fourier_transformed_covariates(covariates, J = J)
    )
  } else {
    stop("Unknown basis_type: ", basis_type, call. = FALSE)
  }

  names(basis_list) <- paste0("J_", J_values)
  structure(
    list(
      wald_covariates = as.matrix(covariates),
      series_covariates = basis_list,
      basis_degrees = J_values
    ),
    class = "gp_transformed_covariates"
  )
}

projection_test <- function(psedo_outcome, covariates) {
  psedo_outcome <- as.numeric(psedo_outcome)
  x_vec <- as.matrix(covariates)
  N <- length(psedo_outcome)
  b_vec <- cbind(1, x_vec)

  m0 <- stats::lm(psedo_outcome ~ b_vec - 1)
  coef_0 <- stats::coef(m0)
  residual_0 <- as.numeric(stats::residuals(m0))
  SXX <- crossprod(b_vec) / N
  weighted_b_vec <- sweep(b_vec, 1, residual_0, FUN = "*")
  meat <- crossprod(weighted_b_vec) / N
  df <- ncol(b_vec)

  statistic <- N * t(coef_0) %*%
    solve(solve(SXX) %*% meat %*% solve(SXX)) %*%
    coef_0

  list(
    statistic = as.numeric(statistic),
    p_value = stats::pchisq(as.numeric(statistic), df = df, lower.tail = FALSE),
    df = df
  )
}

weighted_chisq_pvalue <- function(Sn, M, simu_time = chi_square_simu_time) {
  e_values <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  e_values <- pmax(as.numeric(Re(e_values)), 0)

  null_sample <- matrix(
    stats::rchisq(simu_time * length(e_values), df = 1),
    nrow = simu_time,
    ncol = length(e_values)
  )
  null_sample <- rowSums(sweep(null_sample, 2, e_values, FUN = "*"))
  mean(null_sample >= as.numeric(Sn))
}

bonferroni_pvalue <- function(p_values) {
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0) {
    return(NA_real_)
  }
  min(1, length(p_values) * min(p_values))
}

reject_columns <- function(p_value) {
  reject <- as.integer(!is.na(p_value) & p_value <= alpha_grid)
  names(reject) <- paste0(
    "reject_alpha_",
    gsub("\\.", "_", sprintf("%.2f", alpha_grid))
  )
  as.data.frame(as.list(reject), check.names = FALSE)
}

make_test_row <- function(
    example,
    sample_size,
    scenario_index,
    parameter_1,
    parameter_2,
    method,
    basis_type = NA_character_,
    test_type = NA_character_,
    J = NA_character_,
    combined = FALSE,
    statistic = NA_real_,
    p_value = NA_real_,
    df = NA_real_,
    nuisance_elapsed_seconds = NA_real_,
    method_elapsed_seconds = NA_real_,
    error_message = NA_character_) {
  row <- data.frame(
    replicate = touse,
    example = example,
    sample_size = sample_size,
    scenario_index = scenario_index,
    parameter_1 = parameter_1,
    parameter_2 = parameter_2,
    method = method,
    basis_type = basis_type,
    test_type = test_type,
    J = as.character(J),
    combined = combined,
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    df = as.numeric(df),
    nuisance_elapsed_seconds = as.numeric(nuisance_elapsed_seconds),
    method_elapsed_seconds = as.numeric(method_elapsed_seconds),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  cbind(row, reject_columns(p_value))
}

gp_test_rows <- function(
    example,
    sample_size,
    scenario_index,
    parameter_1,
    parameter_2,
    psedo_outcome,
    covariates,
    J_values,
    basis_type,
    nuisance_elapsed_seconds = NA_real_) {
  method_start_time <- proc.time()[["elapsed"]]

  transformed_covariates <- make_transformed_covariates(
    covariates = covariates,
    J_values = J_values,
    basis_type = basis_type
  )
  W <- gp_test(
    psedo_outcome = psedo_outcome,
    transformed_covariates = transformed_covariates
  )

  rows <- list()
  standardized_p_values <- numeric(length(J_values))
  unstandardized_p_values <- numeric(length(J_values))

  for (j in seq_along(J_values)) {
    standardized_p_values[j] <- stats::pnorm(
      W$stat_series[j],
      lower.tail = FALSE
    )
    unstandardized_p_values[j] <- weighted_chisq_pvalue(
      Sn = W$S_vec[j],
      M = W$M_list[[j]]
    )

    rows[[length(rows) + 1]] <- make_test_row(
      example = example,
      sample_size = sample_size,
      scenario_index = scenario_index,
      parameter_1 = parameter_1,
      parameter_2 = parameter_2,
      method = "gp",
      basis_type = basis_type,
      test_type = "standardized",
      J = J_values[j],
      combined = FALSE,
      statistic = W$stat_series[j],
      p_value = standardized_p_values[j],
      nuisance_elapsed_seconds = nuisance_elapsed_seconds
    )

    rows[[length(rows) + 1]] <- make_test_row(
      example = example,
      sample_size = sample_size,
      scenario_index = scenario_index,
      parameter_1 = parameter_1,
      parameter_2 = parameter_2,
      method = "gp",
      basis_type = basis_type,
      test_type = "unstandardized",
      J = J_values[j],
      combined = FALSE,
      statistic = W$S_vec[j],
      p_value = unstandardized_p_values[j],
      nuisance_elapsed_seconds = nuisance_elapsed_seconds
    )
  }

  rows[[length(rows) + 1]] <- make_test_row(
    example = example,
    sample_size = sample_size,
    scenario_index = scenario_index,
    parameter_1 = parameter_1,
    parameter_2 = parameter_2,
    method = "gp",
    basis_type = basis_type,
    test_type = "standardized",
    J = paste(J_values, collapse = ";"),
    combined = TRUE,
    statistic = NA_real_,
    p_value = bonferroni_pvalue(standardized_p_values),
    nuisance_elapsed_seconds = nuisance_elapsed_seconds
  )

  rows[[length(rows) + 1]] <- make_test_row(
    example = example,
    sample_size = sample_size,
    scenario_index = scenario_index,
    parameter_1 = parameter_1,
    parameter_2 = parameter_2,
    method = "gp",
    basis_type = basis_type,
    test_type = "unstandardized",
    J = paste(J_values, collapse = ";"),
    combined = TRUE,
    statistic = NA_real_,
    p_value = bonferroni_pvalue(unstandardized_p_values),
    nuisance_elapsed_seconds = nuisance_elapsed_seconds
  )

  out <- do.call(rbind, rows)
  out$method_elapsed_seconds <- as.numeric(proc.time()[["elapsed"]] - method_start_time)
  out
}

run_example_1_replicate <- function() {
  alpha_list <- mean_exchangeability_alpha_list()
  rows <- list()

  for (n_name in names(example_1_j_map)) {
    n <- as.integer(n_name)
    J_values <- example_1_j_map[[n_name]]

    for (scenario_index in seq_along(alpha_list)) {
      alpha <- alpha_list[[scenario_index]]
      dataset <- dgp_mean_exchangeability(n = n, alpha = alpha)
      psedo_timing <- time_expression(
        estimate_mean_exchangeability_pseudo_outcome(dataset, V = 2)
      )
      psedo_dat <- psedo_timing$value
      nuisance_elapsed_seconds <- psedo_timing$elapsed_seconds
      psedo_outcome <- psedo_dat$psedo_outcome
      covariates <- psedo_dat$x_vec

      projection_timing <- time_expression(
        projection_test(psedo_outcome, covariates)
      )
      projection <- projection_timing$value
      rows[[length(rows) + 1]] <- make_test_row(
        example = "example_1",
        sample_size = n,
        scenario_index = scenario_index,
        parameter_1 = alpha[1],
        parameter_2 = alpha[2],
        method = "projection",
        basis_type = NA_character_,
        test_type = "wald",
        J = NA_character_,
        combined = FALSE,
        statistic = projection$statistic,
        p_value = projection$p_value,
        df = projection$df,
        nuisance_elapsed_seconds = nuisance_elapsed_seconds,
        method_elapsed_seconds = projection_timing$elapsed_seconds
      )

      for (basis_type in c("legendre", "fourier")) {
        rows[[length(rows) + 1]] <- gp_test_rows(
          example = "example_1",
          sample_size = n,
          scenario_index = scenario_index,
          parameter_1 = alpha[1],
          parameter_2 = alpha[2],
          psedo_outcome = psedo_outcome,
          covariates = covariates,
          J_values = J_values,
          basis_type = basis_type,
          nuisance_elapsed_seconds = nuisance_elapsed_seconds
        )
      }

      racine_timing <- time_expression(
        tryCatch(
          test_conditional_independence_np_example_1(dataset),
          error = function(e) e
        )
      )
      racine_result <- racine_timing$value
      if (inherits(racine_result, "error")) {
        racine_p <- NA_real_
        racine_stat <- NA_real_
        racine_error <- conditionMessage(racine_result)
      } else {
        racine_component <- tryCatch(
          extract_np_significance_component(racine_result, variable_name = "factor(s)"),
          error = function(e) e
        )
        if (inherits(racine_component, "error")) {
          racine_p <- NA_real_
          racine_stat <- NA_real_
          racine_error <- conditionMessage(racine_component)
        } else {
          racine_p <- racine_component$p_value
          racine_stat <- racine_component$statistic
          racine_error <- NA_character_
        }
      }

      rows[[length(rows) + 1]] <- make_test_row(
        example = "example_1",
        sample_size = n,
        scenario_index = scenario_index,
        parameter_1 = alpha[1],
        parameter_2 = alpha[2],
        method = "racine",
        basis_type = NA_character_,
        test_type = "conditional_independence",
        J = NA_character_,
        combined = FALSE,
        statistic = racine_stat,
        p_value = racine_p,
        nuisance_elapsed_seconds = nuisance_elapsed_seconds,
        method_elapsed_seconds = racine_timing$elapsed_seconds,
        error_message = racine_error
      )
    }
  }

  do.call(rbind, rows)
}

run_example_2_replicate <- function() {
  beta_list <- iv_compatibility_beta_list()
  rows <- list()

  for (n_name in names(example_2_j_map)) {
    n <- as.integer(n_name)
    J_values <- example_2_j_map[[n_name]]

    for (scenario_index in seq_along(beta_list)) {
      beta <- beta_list[[scenario_index]]
      dataset <- dgp_iv_compatibility(n = n, beta = beta)
      psedo_timing <- time_expression(
        estimate_iv_compatibility_pseudo_outcome(dataset, V = 2)
      )
      psedo_dat <- psedo_timing$value
      nuisance_elapsed_seconds <- psedo_timing$elapsed_seconds
      psedo_outcome <- psedo_dat$psedo_outcome
      covariates <- psedo_dat$x_vec

      projection_timing <- time_expression(
        projection_test(psedo_outcome, covariates)
      )
      projection <- projection_timing$value
      rows[[length(rows) + 1]] <- make_test_row(
        example = "example_2",
        sample_size = n,
        scenario_index = scenario_index,
        parameter_1 = beta[1],
        parameter_2 = beta[2],
        method = "projection",
        basis_type = NA_character_,
        test_type = "wald",
        J = NA_character_,
        combined = FALSE,
        statistic = projection$statistic,
        p_value = projection$p_value,
        df = projection$df,
        nuisance_elapsed_seconds = nuisance_elapsed_seconds,
        method_elapsed_seconds = projection_timing$elapsed_seconds
      )

      for (basis_type in c("legendre", "fourier")) {
        rows[[length(rows) + 1]] <- gp_test_rows(
          example = "example_2",
          sample_size = n,
          scenario_index = scenario_index,
          parameter_1 = beta[1],
          parameter_2 = beta[2],
          psedo_outcome = psedo_outcome,
          covariates = covariates,
          J_values = J_values,
          basis_type = basis_type,
          nuisance_elapsed_seconds = nuisance_elapsed_seconds
        )
      }
    }
  }

  do.call(rbind, rows)
}

example_1_results <- run_example_1_replicate()
example_2_results <- run_example_2_replicate()

write.csv(
  example_1_results,
  file = paste0("output1/sim_", touse, ".csv"),
  row.names = FALSE
)

write.csv(
  example_2_results,
  file = paste0("output2/sim_", touse, ".csv"),
  row.names = FALSE
)
