## Shared DGPs and testing functions for the two simulation examples.
## Sourcing this file defines functions only; it does not run simulations.

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but is not installed.", package),
         call. = FALSE)
  }
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
  log(x / (1 - x))
}

legendre_orthonormal_matrix <- function(x, k) {
  stopifnot(is.numeric(x), k >= 0, k == as.integer(k))

  n <- length(x)
  M <- matrix(0, n, k + 1)
  M[, 1] <- 1

  if (k >= 1) {
    M[, 2] <- x
  }

  if (k >= 2) {
    for (m in 2:k) {
      M[, m + 1] <- ((2 * m - 1) * x * M[, m] -
                       (m - 1) * M[, m - 1]) / m
    }
  }

  scales <- sqrt((2 * (0:k) + 1) / 2)
  M <- sweep(M, 2, scales, FUN = "*")
  colnames(M) <- paste0("phi", 0:k)
  M
}

legendre_transformed_covariates <- function(
    covariates,
    k,
    include_intercept = TRUE,
    intercept_value = 1 / sqrt(2)) {
  x_vec <- as.matrix(covariates)
  bs_basis <- NULL

  for (j in seq_len(ncol(x_vec))) {
    bsj_basis <- legendre_orthonormal_matrix(x_vec[, j], k)
    if (ncol(bsj_basis) > 1) {
      bs_basis <- cbind(bs_basis, bsj_basis[, -1, drop = FALSE])
    }
  }

  if (include_intercept) {
    bs_basis <- cbind(rep(intercept_value, nrow(x_vec)), bs_basis)
    colnames(bs_basis)[1] <- ""
  }

  bs_basis
}

legendre_gp_transformed_covariates <- function(
    covariates,
    k_vec = c(2, 3, 5),
    include_default = TRUE) {
  x_vec <- as.matrix(covariates)
  N <- nrow(x_vec)
  k_all <- k_vec

  if (include_default) {
    default_k <- ceiling(N^(1 / 3) / log10(N))
    k_all <- c(default_k, k_vec)
  }

  series_covariates <- lapply(
    k_all,
    function(k) legendre_transformed_covariates(x_vec, k = k)
  )
  names(series_covariates) <- paste0("k_", k_all)

  structure(
    list(
      wald_covariates = x_vec,
      series_covariates = series_covariates,
      basis_degrees = k_all
    ),
    class = "gp_transformed_covariates"
  )
}




fit_super_learner <- function(y, x, family, sl_library, cv_folds = 5L) {
  require_package("SuperLearner")

  SuperLearner::SuperLearner(
    Y = y,
    X = as.data.frame(x),
    family = family,
    SL.library = sl_library,
    cvControl = list(V = as.integer(cv_folds), shuffle = TRUE, validRows = NULL)
  )
}

predict_super_learner <- function(fit, newdata) {
  as.numeric(stats::predict(fit, newdata = as.data.frame(newdata), onlySL = TRUE)$pred)
}

make_cv_folds <- function(n, V) {
  require_package("caret")
  caret::createFolds(seq_len(n), V)
}

## ---------------------------------------------------------------------------
## Example 1: mean exchangeability in data fusion
## ---------------------------------------------------------------------------

dgp_mean_exchangeability <- function(n = 1000, alpha = c(0, 0)) {
  x1 <- stats::runif(n, min = -1, max = 1)
  x2 <- stats::runif(n, min = -1, max = 1)

  ps <- expit(x1 - x2)
  s <- stats::rbinom(n, 1, ps)

  pa <- (s == 1) * expit(1.5 * x1 - 0.5 * x2) +
    (s == 0) * expit(x1 + 0.5 * x2)
  a <- stats::rbinom(n, 1, pa)

  y0 <- (s == 0) * (x1 + x2 + expit(x1)) +
    (s == 1) * (
      x1 + x2 + expit(x1) +
        alpha[1] * (cos(pi * x1) + cos(pi * x2)) +
        alpha[2] * (x1 + x2)
    ) +
    0.5 * stats::rnorm(n, 0, 1)
  y1 <- y0 + 2 * (x1 - x2)
  y <- a * y1 + (1 - a) * y0

  data.frame(a, s, x1, x2,y, y0, y1)
}



dgp_mean_exchangeability_high_dim <- function(n = 1000, alpha = c(0, 0)) {
  x1 <- stats::runif(n, min = -1, max = 1)
  x2 <- stats::runif(n, min = -1, max = 1)
  x3 <- stats::runif(n, min = -1, max = 1)
  x4 <- stats::runif(n, min = -1, max = 1)
  x5 <- stats::runif(n, min = -1, max = 1)
  
  x6 <- rbinom(n, size = 1, prob = 0.5)
  x7 <- rbinom(n, size = 1, prob = 0.5)
  x8 <- rbinom(n, size = 1, prob = 0.5)
  x9 <- rbinom(n, size = 1, prob = 0.5)
  x10 <- rbinom(n, size = 1, prob = 0.5)
  
  ps <- expit(x1 - x2 + 0.3*x7 -0.3*x9)
  s <- stats::rbinom(n, 1, ps)
  
  pa <- (s == 1) * expit(1.5 * x1 - 0.5 * x2) +
    (s == 0) * expit(x1 + 0.5 * x2)
  a <- stats::rbinom(n, 1, pa)
  
  y0 <- (s == 0) * (x1 + x2 + expit(x1) + x6 - x7) +
    (s == 1) * (
      x1 + x2 + expit(x1) + x6 - x7 +
        alpha[1] * (cos(pi * x3) + cos(pi * x4)) +
        alpha[2] * (x1 + x2)
    ) +
    0.5 * stats::rnorm(n, 0, 1)
  y1 <- y0 + 2 * (x1 - x2)
  y <- a * y1 + (1 - a) * y0
  
  data.frame(a, s, x1, x2,x3, x4,x5, x6,x7, x8,x9, x10, y, y0, y1)
}

estimate_mean_exchangeability_pseudo_outcome <- function(
    dataset,
    V = 2,
    treatment = 0,
    covariate_names = c("x1", "x2"),
    outcome_library = c("SL.glm", "SL.randomForest"),
    propensity_library = c("SL.glm"),
    cv_folds = 5L) {
  required_names <- c("a", "s", covariate_names, "y")
  missing_names <- setdiff(required_names, names(dataset))
  if (length(missing_names) > 0) {
    stop("dataset is missing required columns: ",
         paste(missing_names, collapse = ", "), call. = FALSE)
  }

  index <- make_cv_folds(nrow(dataset), V)
  pseudo_outcome <- NULL
  x_vec <- NULL

  for (v in seq_len(V)) {
    trainingset <- dataset[-index[[v]], , drop = FALSE]
    predictset <- dataset[index[[v]], , drop = FALSE]
    predict_x <- predictset[, covariate_names, drop = FALSE]

    x_vec <- rbind(x_vec, predict_x)
    # drop = FALSE tells R to preserve the result as a data frame/matrix structure
    trainingset_s0_a <- trainingset[trainingset$s == 0 & trainingset$a == treatment, ,
                                    drop = FALSE]
    trainingset_s1_a <- trainingset[trainingset$s == 1 & trainingset$a == treatment, ,
                                    drop = FALSE]

    sl_y_s0_a <- fit_super_learner(
      y = trainingset_s0_a$y,
      x = trainingset_s0_a[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = outcome_library,
      cv_folds = cv_folds
    )
    pred_y_s0_a <- predict_super_learner(sl_y_s0_a, predict_x)

    sl_y_s1_a <- fit_super_learner(
      y = trainingset_s1_a$y,
      x = trainingset_s1_a[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = outcome_library,
      cv_folds = cv_folds
    )
    pred_y_s1_a <- predict_super_learner(sl_y_s1_a, predict_x)

    sl_s1 <- fit_super_learner(
      y = trainingset$s,
      x = trainingset[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = propensity_library,
      cv_folds = cv_folds
    )
    pred_s1 <- predict_super_learner(sl_s1, predict_x)
    pred_s0 <- 1 - pred_s1

    trainingset_s0 <- trainingset[trainingset$s == 0, , drop = FALSE]
    trainingset_s1 <- trainingset[trainingset$s == 1, , drop = FALSE]

    sl_a1_s1 <- fit_super_learner(
      y = trainingset_s1$a,
      x = trainingset_s1[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = propensity_library,
      cv_folds = cv_folds
    )
    pred_a1_s1 <- predict_super_learner(sl_a1_s1, predict_x)

    sl_a1_s0 <- fit_super_learner(
      y = trainingset_s0$a,
      x = trainingset_s0[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = propensity_library,
      cv_folds = cv_folds
    )
    pred_a1_s0 <- predict_super_learner(sl_a1_s0, predict_x)

    pred_a_s1 <- if (treatment == 1) pred_a1_s1 else 1 - pred_a1_s1
    pred_a_s0 <- if (treatment == 1) pred_a1_s0 else 1 - pred_a1_s0
    pred_as1 <- pred_a_s1 * pred_s1
    pred_as0 <- pred_a_s0 * pred_s0

    ind_as1 <- as.numeric(predictset$a == treatment & predictset$s == 1)
    ind_as0 <- as.numeric(predictset$a == treatment & predictset$s == 0)

    pseudo_a <- (
      ind_as1 / pred_as1 * (predictset$y - pred_y_s1_a) + pred_y_s1_a
    ) - (
      ind_as0 / pred_as0 * (predictset$y - pred_y_s0_a) + pred_y_s0_a
    )
    pseudo_outcome <- c(pseudo_outcome, pseudo_a)
  }

  list(psedo_outcome = pseudo_outcome,
       pseudo_outcome = pseudo_outcome,
       x_vec = x_vec)
}

run_mean_exchangeability_gp_test <- function(
    n = 1000,
    alpha = c(0, 0),
    V = 2,
    k_vec = c(3, 5, 10, 15, 20),
    treatment = 0,
    ...) {
  dataset <- dgp_mean_exchangeability(n = n, alpha = alpha)
  pseudo_dat <- estimate_mean_exchangeability_pseudo_outcome(
    dataset = dataset,
    V = V,
    treatment = treatment,
    ...
  )
  transformed_covariates <- legendre_gp_transformed_covariates(
    covariates = pseudo_dat$x_vec,
    k_vec = k_vec
  )
  gp_test(psedo_outcome = pseudo_dat$psedo_outcome,
          transformed_covariates = transformed_covariates)
}

# The simulation setting for alpha
mean_exchangeability_alpha_list <- function() {
  list(c(0, 0), c(0.2, 0), c(0, 0.2), c(0.4, 0.2), c(0.2, 0.4))
}

mean_exchangeability_sample_size_list <- function() {
  c(250, 500, 1000, 1500, 3000)
}

## Backward-compatible names from example_1_new.R.
dat_gen_1 <- dgp_mean_exchangeability
test_a0_1 <- function(dataset, V = 2) {
  estimate_mean_exchangeability_pseudo_outcome(dataset = dataset, V = V, treatment = 0)
}

## ---------------------------------------------------------------------------
## Example 2: compatibility condition in IV problems
## ---------------------------------------------------------------------------

dgp_iv_compatibility_high_dim <- function(n = 1000, beta = c(0, 0)) {
  x1 <- stats::runif(n, min = -1, max = 1)
  x2 <- stats::runif(n, min = -1, max = 1)
  x3 <- stats::runif(n, min = -1, max = 1)
  x4 <- stats::runif(n, min = -1, max = 1)
  x5 <- stats::runif(n, min = -1, max = 1)
  
  x6 <- rbinom(n, size = 1, prob = 0.5)
  x7 <- rbinom(n, size = 1, prob = 0.5)
  x8 <- rbinom(n, size = 1, prob = 0.5)
  x9 <- rbinom(n, size = 1, prob = 0.5)
  x10 <- rbinom(n, size = 1, prob = 0.5)

  x <- cbind(x1, x2)
  x_star <- 1 * (x > 0)

  pz1 <- expit(0.5 + 0.5 * (x1>0) - 0.5 * (x2>0)+ 0.3*x7 -0.3*x9)
  pz2 <- expit(0.5 + 0.5 * (x1>0) + 0.5 * (x2>0)+ 0.3*x7 -0.3*x9)
  z1 <- stats::rbinom(n, 1, pz1)
  z2 <- stats::rbinom(n, 1, pz2)

  u <- stats::rnorm(n, -0.3, 0.3)

  par_ant <- c(0, -1)
  par_co1 <- c(1, 2)
  par_co2 <- c(1, 2)
  par_rco <- c(2, 2)
  par_eco <- c(2, 1)

  num_ant <- c(exp(1 + x_star %*% par_ant + 0.3 * (u > 0)))
  num_co1 <- c(exp(3.5 + 0.5 * x_star %*% par_co1 + 0.3 * (u > 0)))
  num_co2 <- c(exp(3.5 + 0.5 * x_star %*% par_co2 + 0.3 * (u > 0)))
  num_rco <- c(exp(2 + 0.5 * x_star %*% par_rco + 0.3 * (u > 0)))
  num_eco <- c(exp(2 + 0.5 * x_star %*% par_eco + 0.3 * (u > 0)))

  dem <- num_ant + num_co1 + num_co2 + num_rco + num_eco
  p_matrix <- cbind(
    p_ant = num_ant / dem,
    p_co1 = num_co1 / dem,
    p_co2 = num_co2 / dem,
    p_rco = num_rco / dem,
    p_eco = num_eco / dem
  )

  s <- character(n)
  strata <- c("ANT", "CO1", "CO2", "RCO", "ECO")
  for (i in seq_len(n)) {
    dat_i <- c(stats::rmultinom(1, 1, p_matrix[i, ]))
    s[i] <- strata[which(dat_i == 1)]
  }

  d <- numeric(n)
  d[s == "ANT"] <- 0
  d[s == "CO1" & z1 == 1] <- 1
  d[s == "CO1" & z1 == 0] <- 0
  d[s == "CO2" & z2 == 1] <- 1
  d[s == "CO2" & z2 == 0] <- 0
  d[s == "RCO"] <- 0
  d[s == "RCO" & z1 == 1 & z2 == 1] <- 1
  d[s == "ECO"] <- 0
  d[s == "ECO" & (z1 == 1 | z2 == 1)] <- 1

  epsilon <- stats::rnorm(n, 0, 1)
  y_0 <- 1 * (s == "ANT") * (1 + x1 + x2 + 0.3*x6 -0.3*x7+u) +
    1 * (s == "CO1") * (1 + x1 + x2 + 0.3*x6 -0.3*x7 +  u) +
    1 * (s == "CO2") * (1 + x1 + x2+ 0.3*x6 -0.3*x7 +u) +
    1 * (s == "RCO") * (1 + x1 + x2 + 0.3*x6 -0.3*x7 +u) +
    1 * (s == "ECO") * (1 + x1 + x2 + 0.3*x6 -0.3*x7 + u) +
    epsilon

  y_1 <- 1 * (s == "ANT") * (1 + x1+ x2 + u) +
    1 * (s == "CO1") * (1 - x1+ x2 + u) +
    1 * (s == "CO2") * (
      1 - x1+ x2 +
        beta[1] * (cos(pi * x3) + cos(pi * x4)) +
        beta[2] * (x1 + x2) +
        u
    ) +
    1 * (s == "RCO") * (1 - x1+ x2 + u) +
    1 * (s == "ECO") * (1 - x1+ x2 + u) +
    epsilon

  y <- (1 - d) * y_0 + d * y_1
  ite <- y_1 - y_0

  data.frame(x1, x2, x3, x4,x5, x6,x7, x8,x9, x10, d, z2, z1, y, y_0, y_1, s, u, ite)
}



dgp_iv_compatibility <- function(n = 1000, beta = c(0, 0)) {
  x1 <- stats::runif(n, min = -1, max = 1)
  x2 <- stats::runif(n, min = -1, max = 1)

  
  x <- cbind(x1, x2)
  x_star <- 1 * (x > 0)
  
  pz1 <- expit(0.5 + 0.5 * (x1>0) - 0.5 * (x2>0) )
  pz2 <- expit(0.5 + 0.5 * (x1>0) + 0.5 * (x2>0))
  z1 <- stats::rbinom(n, 1, pz1)
  z2 <- stats::rbinom(n, 1, pz2)
  
  u <- stats::rnorm(n, -0.3, 0.3)
  
  par_ant <- c(0, -1)
  par_co1 <- c(1, 2)
  par_co2 <- c(1, 2)
  par_rco <- c(2, 2)
  par_eco <- c(2, 1)
  
  num_ant <- c(exp(1 + x_star %*% par_ant + 0.3 * (u > 0)))
  num_co1 <- c(exp(3.5 + 0.5 * x_star %*% par_co1 + 0.3 * (u > 0)))
  num_co2 <- c(exp(3.5 + 0.5 * x_star %*% par_co2 + 0.3 * (u > 0)))
  num_rco <- c(exp(2 + 0.5 * x_star %*% par_rco + 0.3 * (u > 0)))
  num_eco <- c(exp(2 + 0.5 * x_star %*% par_eco + 0.3 * (u > 0)))
  
  dem <- num_ant + num_co1 + num_co2 + num_rco + num_eco
  p_matrix <- cbind(
    p_ant = num_ant / dem,
    p_co1 = num_co1 / dem,
    p_co2 = num_co2 / dem,
    p_rco = num_rco / dem,
    p_eco = num_eco / dem
  )
  
  s <- character(n)
  strata <- c("ANT", "CO1", "CO2", "RCO", "ECO")
  for (i in seq_len(n)) {
    dat_i <- c(stats::rmultinom(1, 1, p_matrix[i, ]))
    s[i] <- strata[which(dat_i == 1)]
  }
  
  d <- numeric(n)
  d[s == "ANT"] <- 0
  d[s == "CO1" & z1 == 1] <- 1
  d[s == "CO1" & z1 == 0] <- 0
  d[s == "CO2" & z2 == 1] <- 1
  d[s == "CO2" & z2 == 0] <- 0
  d[s == "RCO"] <- 0
  d[s == "RCO" & z1 == 1 & z2 == 1] <- 1
  d[s == "ECO"] <- 0
  d[s == "ECO" & (z1 == 1 | z2 == 1)] <- 1
  
  epsilon <- stats::rnorm(n, 0, 1)
  y_0 <- 1 * (s == "ANT") * (1 + x1+ x2 + u) +
    1 * (s == "CO1") * (1 + x1+ x2 + u) +
    1 * (s == "CO2") * (1 + x1+ x2 + u) +
    1 * (s == "RCO") * (1 + x1+ x2 + u) +
    1 * (s == "ECO") * (1 + x1+ x2 + u) +
    epsilon
  
  y_1 <- 1 * (s == "ANT") * (1 + x1+ x2 + u) +
    1 * (s == "CO1") * (1 - x1+ x2 + u) +
    1 * (s == "CO2") * (
      1 - x1+ x2 +
        beta[1] * (cos(pi * x1) + cos(pi * x2)) +
        beta[2] * (x1 + x2) +
        u
    ) +
    1 * (s == "RCO") * (1 - x1+ x2 + u) +
    1 * (s == "ECO") * (1 - x1 + x2 + u) +
    epsilon
  
  y <- (1 - d) * y_0 + d * y_1
  ite <- y_1 - y_0
  
  data.frame(x1, x2, d, z2, z1, y, y_0, y_1, s, u, ite)
}

estimate_iv_compatibility_pseudo_outcome <- function(
    dataset,
    V = 2,
    covariate_names = c("x1", "x2"),
    sl_library = NULL,
    z_library = c("SL.glm", "SL.randomForest"),
    d_library = c("SL.glm", "SL.randomForest"),
    y_library = c("SL.glm", "SL.randomForest", "SL.gam"),
    cv_folds = 5L) {
  if (!is.null(sl_library)) {
    z_library <- sl_library
    d_library <- sl_library
    y_library <- sl_library
  }

  required_names <- c(covariate_names, "d", "z1", "z2", "y")
  missing_names <- setdiff(required_names, names(dataset))
  if (length(missing_names) > 0) {
    stop("dataset is missing required columns: ",
         paste(missing_names, collapse = ", "), call. = FALSE)
  }

  index <- make_cv_folds(nrow(dataset), V)
  pseudo_outcome <- NULL
  x_vec <- NULL

  for (v in seq_len(V)) {
    trainingset <- dataset[-index[[v]], , drop = FALSE]
    predictset <- dataset[index[[v]], , drop = FALSE]
    predict_x <- predictset[, covariate_names, drop = FALSE]
    x_vec <- rbind(x_vec, predict_x)

    sl_z1 <- fit_super_learner(
      y = trainingset$z1,
      x = trainingset[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = z_library,
      cv_folds = cv_folds
    )
    
    sl_z2 <- fit_super_learner(
      y = trainingset$z2,
      x = trainingset[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = z_library,
      cv_folds = cv_folds
    )
    
    # trainingset_z1_1 <- trainingset[trainingset$z1 == 1, , drop = FALSE]
    # trainingset_z1_0 <- trainingset[trainingset$z1 == 0, , drop = FALSE]
    # 
    # sl_z2_z1_1 <- fit_super_learner(
    #   y = trainingset_z1_1$z2,
    #   x = trainingset_z1_1[, covariate_names, drop = FALSE],
    #   family = stats::binomial(),
    #   sl_library = z_library,
    #   cv_folds = cv_folds
    # )
    # sl_z2_z1_0 <- fit_super_learner(
    #   y = trainingset_z1_0$z2,
    #   x = trainingset_z1_0[, covariate_names, drop = FALSE],
    #   family = stats::binomial(),
    #   sl_library = z_library,
    #   cv_folds = cv_folds
    # )
    
    z1_1_pred <- (predict_super_learner(sl_z1, predict_x))
    z1_0_pred <- (1 - z1_1_pred)
    # z2_1_z1_1_pred <- (predict_super_learner(sl_z2_z1_1, predict_x))
    # z2_1_z1_0_pred <- (predict_super_learner(sl_z2_z1_0, predict_x))
    # 
    # z2_z1_1_1_pred <- z1_1_pred * z2_1_z1_1_pred
    # z2_z1_1_0_pred <- z1_0_pred * z2_1_z1_0_pred
    # z2_z1_0_1_pred <- z1_1_pred * (1 - z2_1_z1_1_pred)
    # z2_z1_0_0_pred <- z1_0_pred * (1 - z2_1_z1_0_pred)
    
    z2_1_pred <- (predict_super_learner(sl_z2, predict_x))
    z2_0_pred <- (1 - z2_1_pred)

    trainingset_co1_1 <- trainingset[trainingset$z1 == 1, , drop = FALSE]
    trainingset_co1_0 <- trainingset[trainingset$z1 == 0, , drop = FALSE]
    trainingset_co2_1 <- trainingset[trainingset$z2 == 1, , drop = FALSE]
    trainingset_co2_0 <- trainingset[trainingset$z2 == 0, , drop = FALSE]

    sl_d_co1_1 <- fit_super_learner(
      y = trainingset_co1_1$d,
      x = trainingset_co1_1[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = d_library,
      cv_folds = cv_folds
    )
    sl_d_co1_0 <- fit_super_learner(
      y = trainingset_co1_0$d,
      x = trainingset_co1_0[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = d_library,
      cv_folds = cv_folds
    )
    sl_d_co2_1 <- fit_super_learner(
      y = trainingset_co2_1$d,
      x = trainingset_co2_1[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = d_library,
      cv_folds = cv_folds
    )
    sl_d_co2_0 <- fit_super_learner(
      y = trainingset_co2_0$d,
      x = trainingset_co2_0[, covariate_names, drop = FALSE],
      family = stats::binomial(),
      sl_library = d_library,
      cv_folds = cv_folds
    )

    sl_y_co1_1 <- fit_super_learner(
      y = trainingset_co1_1$y,
      x = trainingset_co1_1[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = y_library,
      cv_folds = cv_folds
    )
    sl_y_co1_0 <- fit_super_learner(
      y = trainingset_co1_0$y,
      x = trainingset_co1_0[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = y_library,
      cv_folds = cv_folds
    )
    sl_y_co2_1 <- fit_super_learner(
      y = trainingset_co2_1$y,
      x = trainingset_co2_1[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = y_library,
      cv_folds = cv_folds
    )
    sl_y_co2_0 <- fit_super_learner(
      y = trainingset_co2_0$y,
      x = trainingset_co2_0[, covariate_names, drop = FALSE],
      family = stats::gaussian(),
      sl_library = y_library,
      cv_folds = cv_folds
    )

    d_co1_1_pred <- predict_super_learner(sl_d_co1_1, predict_x)
    d_co1_0_pred <- predict_super_learner(sl_d_co1_0, predict_x)
    d_co2_1_pred <- predict_super_learner(sl_d_co2_1, predict_x)
    d_co2_0_pred <- predict_super_learner(sl_d_co2_0, predict_x)

    y_co1_1_pred <- predict_super_learner(sl_y_co1_1, predict_x)
    y_co1_0_pred <- predict_super_learner(sl_y_co1_0, predict_x)
    y_co2_1_pred <- predict_super_learner(sl_y_co2_1, predict_x)
    y_co2_0_pred <- predict_super_learner(sl_y_co2_0, predict_x)

    psedo_outcome_co1 <- 1 / (d_co1_1_pred - d_co1_0_pred) * (
      as.numeric(predictset$z1 == 1) / z1_1_pred * (predictset$y - y_co1_1_pred) -
        as.numeric(predictset$z1 == 0) / z1_0_pred * (predictset$y - y_co1_0_pred)
    ) -
      (y_co1_1_pred - y_co1_0_pred) / (d_co1_1_pred - d_co1_0_pred)^2 * (
        as.numeric(predictset$z1 == 1) / z1_1_pred * (predictset$d - d_co1_1_pred) -
          as.numeric(predictset$z1 == 0) / z1_0_pred * (predictset$d - d_co1_0_pred)
      ) +
      (y_co1_1_pred - y_co1_0_pred) / (d_co1_1_pred - d_co1_0_pred)

    psedo_outcome_co2 <- 1 / (d_co2_1_pred - d_co2_0_pred) * (
      as.numeric(predictset$z2 == 1) / z2_1_pred * (predictset$y - y_co2_1_pred) -
        as.numeric(predictset$z2 == 0) / z2_0_pred * (predictset$y - y_co2_0_pred)
    ) -
      (y_co2_1_pred - y_co2_0_pred) / (d_co2_1_pred - d_co2_0_pred)^2 * (
        as.numeric(predictset$z2 == 1) / z2_1_pred * (predictset$d - d_co2_1_pred) -
          as.numeric(predictset$z2 == 0) / z2_0_pred * (predictset$d - d_co2_0_pred)
      ) +
      (y_co2_1_pred - y_co2_0_pred) / (d_co2_1_pred - d_co2_0_pred)

    psedo_outcome_diff <- psedo_outcome_co1 - psedo_outcome_co2

    pseudo_outcome <- c(pseudo_outcome, psedo_outcome_diff)
  }

  rownames(x_vec) <- NULL
  list(psedo_outcome = pseudo_outcome,
       pseudo_outcome = pseudo_outcome,
       x_vec = x_vec)
}

run_iv_compatibility_gp_test <- function(
    n = 1000,
    beta = c(0, 0),
    V = 2,
    k_vec = c(3, 5, 10, 15, 20),
    ...) {
  dataset <- dgp_iv_compatibility(n = n, beta = beta)
  pseudo_dat <- estimate_iv_compatibility_pseudo_outcome(dataset = dataset, V = V, ...)
  transformed_covariates <- legendre_gp_transformed_covariates(
    covariates = pseudo_dat$x_vec,
    k_vec = k_vec
  )
  gp_test(psedo_outcome = pseudo_dat$psedo_outcome,
          transformed_covariates = transformed_covariates)
}

iv_compatibility_beta_list <- function() {
  list(c(0, 0), c(0.3, 0), c(0, 0.3), c(0.6, 0.3), c(0.3, 0.6))
}

iv_compatibility_sample_size_list <- function() {
  c(1000, 2000, 4000, 8000)
}

## Backward-compatible names from example_2_new.R.
dat_gen <- dgp_iv_compatibility
nptest <- function(dataset, V = 2) {
  estimate_iv_compatibility_pseudo_outcome(dataset = dataset, V = V)
}

## ---------------------------------------------------------------------------
## Projection tests and summaries
## ---------------------------------------------------------------------------


# a helper function that converts an object into a numeric matrix 
# and checks that it has the expected number of rows
as_numeric_matrix <- function(x, n, argument_name) {
  matrix_x <- as.matrix(x)
  if (nrow(matrix_x) != n) {
    stop(sprintf(
      "%s has %s rows, but psedo_outcome has length %s.",
      argument_name,
      nrow(matrix_x),
      n
    ), call. = FALSE)
  }
  storage.mode(matrix_x) <- "numeric"
  matrix_x
}


# standardization / input-cleaning function for gp_test
# accepts transformed_covariates in several possible formats and converts them into one consistent structure
# it does not statistically normalize or standardize the covariates.
normalize_gp_transformed_covariates <- function(transformed_covariates, n) {
  wald_covariates <- NULL
  basis_degrees <- NULL

  if (inherits(transformed_covariates, "gp_transformed_covariates") ||
      (!is.data.frame(transformed_covariates) &&
       is.list(transformed_covariates) &&
       ("series_covariates" %in% names(transformed_covariates)))) {
    series_covariates <- transformed_covariates$series_covariates
    wald_covariates <- transformed_covariates$wald_covariates
    basis_degrees <- transformed_covariates$basis_degrees
  } else if (!is.data.frame(transformed_covariates) &&
             is.list(transformed_covariates)) {
    series_covariates <- transformed_covariates
  } else {
    series_covariates <- list(transformed_covariates)
  }

  if (!is.list(series_covariates) || length(series_covariates) == 0) {
    stop("transformed_covariates must contain at least one transformed covariate matrix.",
         call. = FALSE)
  }

  series_covariates <- lapply(
    seq_along(series_covariates),
    function(j) as_numeric_matrix(series_covariates[[j]], n,
                                  sprintf("transformed_covariates[[%s]]", j))
  )

  if (!is.null(wald_covariates)) {
    wald_covariates <- as_numeric_matrix(wald_covariates, n, "wald_covariates")
  }

  list(
    wald_covariates = wald_covariates,
    series_covariates = series_covariates,
    basis_degrees = basis_degrees
  )
}

gp_test <- function(
    psedo_outcome,
    transformed_covariates,
    wald_covariates = NULL,
    wald_include_intercept = TRUE,
    basis_degrees = NULL) {
  psedo_outcome <- as.numeric(psedo_outcome)
  N <- length(psedo_outcome)
  transformed <- normalize_gp_transformed_covariates(transformed_covariates, N)

  if (is.null(wald_covariates)) {
    wald_covariates <- transformed$wald_covariates
  } else {
    wald_covariates <- as_numeric_matrix(wald_covariates, N, "wald_covariates")
  }

  if (is.null(basis_degrees)) {
    basis_degrees <- transformed$basis_degrees
  }

  if (!is.null(wald_covariates)) {
    b_vec <- if (wald_include_intercept) cbind(1, wald_covariates) else wald_covariates
  } else {
    b_vec <- transformed$series_covariates[[1]]
  }

  m0 <- stats::lm(psedo_outcome ~ b_vec - 1)
  coef_0 <- stats::coef(m0)
  residual_0 <- as.numeric(stats::residuals(m0))
  SXX <- crossprod(b_vec) / N
  weighted_b_vec <- sweep(b_vec, 1, residual_0, FUN = "*")
  meat <- crossprod(weighted_b_vec) / N
  p_0 <- ncol(b_vec)

  W_wald <- N * t(coef_0) %*%
    solve(solve(SXX) %*% meat %*% solve(SXX)) %*%
    coef_0

  W_vec <- NULL
  S_vec <- NULL
  M_list <- list()

  for (bs_basis in transformed$series_covariates) {
    projection <- as.numeric(t(psedo_outcome) %*% bs_basis / N)
    sn <- N * crossprod(projection)
    M <- crossprod(bs_basis * as.numeric(psedo_outcome)) / N

    an <- sum(diag(M))
    bn <- norm(M, type = "F")
    W_vec <- c(W_vec, as.numeric((sn - an) / (sqrt(2) * bn)))
    S_vec <- c(S_vec, as.numeric(sn))
    M_list <- append(M_list, list(M))
  }

  list(stat_Wald = W_wald,
       stat_series = W_vec,
       S_vec = S_vec,
       M_list = M_list,
       wald_df = p_0,
       basis_degrees = basis_degrees)
}

test_summary_wald <- function(W, p) {
  as.numeric(c(W) > stats::qchisq(seq(0, 1, 0.05), p))
}

test_summary_series <- function(W) {
  as.numeric(c(W) > stats::qnorm(seq(0, 1, 0.05)))
}

test_summary_series2 <- function(Sn, M, simu_time = 1000) {
  e_values <- eigen(M, only.values = TRUE)$values
  chi_square_sample <- replicate(
    simu_time,
    sum(stats::rchisq(length(e_values), df = 1) * e_values)
  )
  as.numeric(Sn > stats::quantile(chi_square_sample, seq(0, 1, 0.05)))
}

test_conditional_independence_np_example_1 <- function(dataset) {
  require_package("np")

  dt_0 <- dataset[dataset$a == 0, , drop = FALSE]
  bw <- np::npregbw(
    formula = y ~ factor(s) + x1 + x2,
    data = dt_0,
    regtype = "ll",
    bwmethod = "cv.aic"
  )
  np::npsigtest(bws = bw, index = 1)
}

extract_np_significance_component <- function(sigtest_result, variable_name) {
  xnames <- sigtest_result$bws$xnames
  if (is.null(xnames)) {
    stop("The np significance-test result does not contain bandwidth variable names.",
         call. = FALSE)
  }

  variable_index <- which(xnames == variable_name)
  if (length(variable_index) != 1) {
    stop(sprintf(
      "Expected exactly one '%s' component in npsigtest output, found %s. Available variables: %s.",
      variable_name,
      length(variable_index),
      paste(xnames, collapse = ", ")
    ), call. = FALSE)
  }

  list(
    variable = variable_name,
    index = variable_index,
    statistic = as.numeric(sigtest_result$In[variable_index]),
    p_value = as.numeric(sigtest_result$P[variable_index])
  )
}
