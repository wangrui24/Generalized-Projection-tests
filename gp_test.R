

#### A function for generating legendre orthonormal basis
legendre_orthonormal_matrix <- function(x, k) {
  stopifnot(is.numeric(x), k >= 0, k == as.integer(k))
  
  n  <- length(x)
  M  <- matrix(0, n, k + 1)
  M[, 1] <- 1                      # P0
  
  if (k >= 1) M[, 2] <- x          # P1
  
  if (k >= 2) {
    for (m in 2:k) {               # recurrence for Pm
      M[, m + 1] <- ((2*m + 1) * x * M[, m] - (m) * M[, m - 1]) / (m+1)
    }
  }
  
  # # rescale columns to orthonormal basis
  scales <- sqrt((2*(0:k) + 1) / 2)
  M <- sweep(M, 2, scales, FUN = "*")
  
  colnames(M) <- paste0("phi", 0:k)  # φ₀, φ₁, …
  M
}


### The function for generalized projection test
### Including the projection test, chi-square approximation, and the normal approximation

gp_test <- function(psedo_outcome, covariates, k_vec = c(2,3,5)){
  
  ####
  #### The projection test
  ####
  
  N <- length(psedo_outcome)
  
  x_vec <- as.matrix(covariates)  # The covariate matrix
  m0 <- lm(psedo_outcome~x_vec) # Do the linear regression model
  coef_0 <- coef(m0)
  pred_0 <- predict(m0)
  
  b_vec <- cbind(1,x_vec) # include the intercept term
  SXX <- t(b_vec)%*%b_vec/N
  meat <- t(b_vec)%*%diag((psedo_outcome)^2)%*%b_vec/N
  p_0 <- ncol(b_vec) # degree of freedom
  
  W_wald <-  N*t(coef_0)%*%solve(solve(SXX)%*%meat%*%solve(SXX))%*%coef_0 # The Wald statistic
  
  ####
  #### The projection test, with the default number of basis
  ####
  
  k <- ceiling(((length(psedo_outcome)))^{1/3}/(log10(length(psedo_outcome))))
  W_vec <- c()
  
  bs_basis <- NULL
  for (j in c(1:ncol(x_vec))) {
    bsj_basis <- legendre_orthonormal_matrix(x_vec[,j],k)[,-1]
    bs_basis <- cbind(bs_basis, bsj_basis)
  }

  bs_basis <- cbind(1/sqrt(2), bs_basis)
  sn <- N*(t(psedo_outcome)%*%bs_basis/N)%*%t(t(psedo_outcome)%*%bs_basis/N) # S_n without scaling by n
  M <-  t(bs_basis)%*%diag(psedo_outcome^2)%*%(bs_basis)/N
  
  an <- sum(diag(M))
  bn <- norm(M, type = "F")
  W_vec <- c(W_vec, ((sn-an)/(sqrt(2)*bn)))
  S_vec <- c(sn)
  M_list <- list(M)
  
  
  for (k in k_vec) {
    bs_basis <- NULL
    for (j in c(1:ncol(x_vec))) {
      bsj_basis <- legendre_orthonormal_matrix(x_vec[,j],k)[,-1]
      bs_basis <- cbind(bs_basis, bsj_basis)
    }
    
    bs_basis <- cbind(1/sqrt(2), bs_basis)
    sn <- N*(t(psedo_outcome)%*%bs_basis/N)%*%t(t(psedo_outcome)%*%bs_basis/N) # S_n without scaling by n
    M <-  t(bs_basis)%*%diag(psedo_outcome^2)%*%(bs_basis)/N
    
    an <- sum(diag(M))
    bn <- norm(M, type = "F")
    W_vec <- c(W_vec, ((sn-an)/(sqrt(2)*bn)))
    S_vec <- c(S_vec,sn)
    M_list <- append(M_list,list(M))
  }
  return(list(stat_Wald = W_wald, stat_series = W_vec, S_vec = S_vec, M_list = M_list))
}




################ Summary function for projection test
test_summary_wald <- function(W,p){
  t <- 1*(c(W)>qchisq(seq(0,1,0.05),p))
  return(c(t))
}

################ Summary function for generalized projection test with normal approximation
test_summary_series <- function(W){
  t <- 1*(c(W)>qnorm(seq(0,1,0.05)))
  return(c(t))
}


################ Summary function for generalized projection test with chi-square approximation
test_summary_series2 <- function(Sn, M, simu_time = 1000){
  # Eigen values of M
  
  
  e_values <- eigen(M)$values
  sample <- c()
  for (j in c(1:simu_time)) {
    sample <- c(sample,sum(rchisq(length(e_values),df = 1)*e_values))
  }
  t <- c(1*(Sn>quantile(sample, seq(0,1,0.05))))
  
  return(t)
}



