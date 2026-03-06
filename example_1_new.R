touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(touse+1000)

source("gp_test.R")
library("MASS")
library("SuperLearner")
library("caret")
library("fda")
library("np")

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

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

# legendre_orthonormal_matrix(c(1),3)
# legendre_orthonormal_matrix(c(-1,-0.5,0,0.5,1),3)

dat_gen_1 <- function(n = 1000, alpha = c(0,0)){
  
  x1 <- runif(n, min = -1, max = 1)
  x2 <- runif(n, min = -1, max = 1)
  
  ps <- expit(x1-x2)
  
  s <- rbinom(n, 1, ps)
  
  pa <- 1*(s==1)*expit(1.5*x1-0.5*x2)+1*(s==0)*expit(1*x1+0.5*x2)
  
  a <- rbinom(n, 1, pa)
  
  y0 <- 1*(s==0)*(x1+x2+expit(x1)) + 1*(s==1)*(x1+x2+expit(x1)+alpha[1]*(cos(pi*x1)+cos(pi*x2))+alpha[2]*(x1+x2))+ 0.5*rnorm(n, 0, 1)
  y1 <- y0+2*(x1-x2)
  y <- a*y1+(1-a)*y0
  
  dat <- data.frame(a,s,x1,x2,y,y0,y1)
  return(dat)
}


test_a0_1 <-  function(dataset,V=2){
  # dataset <-  dat_gen_1(n = 3000)
  dataset$a0s1 <- 1*(dataset$a == 0)*(dataset$s == 1)
  dataset$a0s0 <- 1*(dataset$a == 0)*(dataset$s == 0)
  dataset$a1s1 <- 1*(dataset$a == 1)*(dataset$s == 1)
  dataset$a1s0 <- 1*(dataset$a == 1)*(dataset$s == 0)
  
  # V=2
  
  
  
  N <- nrow(dataset)
  index <- createFolds(1:nrow(dataset),V)
  if_vec <- NULL
  
  x_vec <- NULL
  psedo_a0_vec <- NULL
  
  for(v in 1:V)
  {
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    x_vec <- rbind(x_vec,predictset[,c("x1","x2")])
    
    trainingset_s0_a0 <- trainingset[which((trainingset$s==0)&(trainingset$a==0)),]
    trainingset_s1_a0 <- trainingset[which((trainingset$s==1)&(trainingset$a==0)),]
    
    
    
    
    sl_y_s0_a0 <- SuperLearner(Y = trainingset_s0_a0$y,
                               X = trainingset_s0_a0[,c("x1","x2")],
                               family = gaussian(),
                               SL.library = c("SL.glm", "SL.randomForest"),
                               cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_y_s0_a0 <- predict(sl_y_s0_a0, newdata = predictset, onlySL = T)$pred
    
    sl_y_s1_a0 <- SuperLearner(Y = trainingset_s1_a0$y,
                               X = trainingset_s1_a0[,c("x1","x2")],
                               family = gaussian(),
                               SL.library = c("SL.glm", "SL.randomForest"),
                               cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_y_s1_a0 <- predict(sl_y_s1_a0, newdata = predictset, onlySL = T)$pred
    
    sl_s1 <- SuperLearner(Y = trainingset$s,
                          X = trainingset[,c("x1","x2")],
                          family = binomial(),
                          SL.library = c("SL.glm"),
                          cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_s1 <- predict(sl_s1, newdata = predictset, onlySL = T)$pred
    pred_s0 <- 1-pred_s1
    
    
    trainingset_s0 <- trainingset[which((trainingset$s==0)),]
    trainingset_s1 <- trainingset[which((trainingset$s==1)),]
    sl_a1s1 <- SuperLearner(Y = trainingset_s1$a,
                            X = trainingset_s1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_a1_s1 <- predict(sl_a1s1, newdata = predictset, onlySL = T)$pred
    pred_a0_s1 <- 1 - pred_a1_s1
    
    sl_a1s0 <- SuperLearner(Y = trainingset_s0$a,
                            X = trainingset_s0[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_a1_s0 <- predict(sl_a1s0, newdata = predictset, onlySL = T)$pred
    pred_a0_s0 <- 1-pred_a1_s0
    
    pred_a0s1 <- pred_a0_s1*pred_s1
    pred_a0s0 <- pred_a0_s0*pred_s0
    
    #### Calculate psedo-outcome
    psedo_a0 <- (predictset$a0s1/pred_a0s1*(predictset$y-pred_y_s1_a0)+pred_y_s1_a0)-
      (predictset$a0s0/pred_a0s0*(predictset$y-pred_y_s0_a0)+pred_y_s0_a0)
    psedo_a0_vec <- c(psedo_a0_vec, psedo_a0)
  }
  
  return(list(psedo_outcome = psedo_a0_vec, x_vec = x_vec))
}


alpha_list <- list(c(0,0), c(0.2,0), c(0.2,0.2), c(0.2,0.4))
sample_size_list <- c(250,500,1000,1500)
simu_res = NULL
for (n in c(sample_size_list)){
  for (i in c(1:length(alpha_list))) {
    dt = dat_gen_1(n, alpha = alpha_list[[i]])
    psedo_dat <- test_a0_1(dt)
    W = gp_test(psedo_outcome = psedo_dat$psedo_outcome, covariates = psedo_dat$x_vec,
                k_vec <- ceiling(((n))^{1/5}/(log10(n))))
    W_wald <- W$stat_Wald
    W_series <- W$stat_series
    W_Sn <- W$S_vec
    W_M <- W$M_list
    
    simu_res <- rbind(simu_res,c(test_summary_wald(W_wald,3),n,i,0))
    for (k in c(1:length(W_series))) {
      simu_res <- rbind(simu_res,c(test_summary_series(W_series[k]),n,i,(k)))
      # simu_res <- rbind(simu_res,c(test_summary_series2(W_Sn[k],W_M[[k]]),n,i,(10+k)))
    }
    
    dt_0 <- dt[dt$a==0,]
    
    bw <- npregbw(formula=dt_0$y~factor(dt_0$s)+dt_0$x1+dt_0$x2,
                  regtype="ll",bwmethod="cv.aic")
    res_np <- npsigtest(bws=bw)
    
    simu_res <- rbind(simu_res,c(1*(res_np$P[1]<(1-seq(0,1,0.05))),n,i,(99)))
  }
}


FileSave <- paste0("/home/bzhang3/Wang_Rui/np_test/output1/sim_",touse,".csv")

write.table(simu_res, file=FileSave, row.names = FALSE)
