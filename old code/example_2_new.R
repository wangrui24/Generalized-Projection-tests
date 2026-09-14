touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(touse+2000)

library(SuperLearner)
library(caret)
library(MASS)
source("gp_test.R")


logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}


dat_gen <- function(n=1000,beta = c(0,0)){
  # n = 1000
  x1 <- runif(n, min = -1, max = 1)
  x2 <- runif(n, min = -1, max = 1)
  
  x <- cbind(x1,x2)
  
  x_star <- 1*(x>0)
  
  pz1 <- expit(0.5+0.5*x1-0.5*x2)
  pz2 <- expit(0.5+0.5*x2+0.5*x2)
  z1 <- rbinom(n,1,pz1)
  z2 <- rbinom(n,1,pz2)
  
  # Generate unmeasured confounder
  u <- rnorm(n, -0.3, 0.3)
  
  # Generate principal stratum
  par_ant <- c(0,-1)
  par_co1 <- c(1,2)
  par_co2 <- c(1,2)
  par_rco <- c(2,2)
  par_eco <- c(2,1)
  # par_aco <- c(0.5,1)
  
  
  # set the parameter
  num_ant <- c(exp(1+x_star%*%par_ant+0.3*(u>0)))
  num_co1 <- c(exp(3.5+0.5*x_star%*%par_co1+0.3*(u>0)))
  num_co2 <- c(exp(3.5+0.5*x_star%*%par_co2+0.3*(u>0)))
  num_rco <- c(exp(2+0.5*x_star%*%par_rco+0.3*(u>0)))
  num_eco <- c(exp(2+0.5*x_star%*%par_eco+0.3*(u>0)))
  # num_aco <- c(exp(2.5+x_star%*%par_aco+0.1*(u>0)))
  
  
  # num_ant <- 0
  # num_swi_1 <- c(exp(alpha[1]+alpha[2]*x%*%par_swi_1+alpha[3]*u))
  # num_swi_2 <- 0
  # num_atnt <- 0
  # num_aco <- c(exp(3+x%*%par_aco+0.5*u))
  # num_ntat <- 0
  # num_aat <- 0
  
  dem <- num_ant + num_co1 + num_co2 + num_rco + num_eco# + num_aco
  p_ant <- num_ant/dem
  p_co1 <- num_co1/dem
  p_co2 <- num_co2/dem
  p_rco <- num_rco/dem
  p_eco <- num_eco/dem
 # p_aco <- num_aco/dem
  
  p_matrix <- cbind(p_ant, p_co1, p_co2,p_rco,p_eco)
  
  s <- c()
  strata <- c("ANT", "CO1", "CO2","RCO","ECO")
  for (i in c(1:n)) {
    p_i <- c(p_matrix[i,])
    dat_i <- c(rmultinom(1,1,p_i))
    s[i] <- strata[which(dat_i==1)]
  }
  
  d <- c()
  d[(s=="ANT")]<-0
  d[(s=="CO1")&(z1==1)]<-1
  d[(s=="CO1")&(z1==0)]<-0
  
  d[(s=="CO2")&(z2==1)]<-1
  d[(s=="CO2")&(z2==0)]<-0
  
  d[(s=="RCO")]<-0
  d[(s=="RCO")&((z1==1)&(z2==1))]<-1
 
  
  d[(s=="ECO")]<-0
  d[(s=="ECO")&((z1==1)|(z2==1))]<-1
  
  
  epsilon <- rnorm(n,0,1)
  # Generate potential outcomes of Y
  y_0 <- 1*(s=="ANT")*(1+x1+x2+u)+
    1*(s=="CO1")*(1+x1+x2+u)+
    1*(s=="CO2")*(1+x1+x2+u)+
    1*(s=="RCO")*(1+x1+x2+u)+
    1*(s=="ECO")*(1+x1+x2+u)+epsilon
  
  
  
  
  # y_1 <- 1*(s=="ANT"|s=="AAT")*(1+x1+2*x2+2*x3+u)+
  #   1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+2*x3+u)+
  #   1*(s=="SWI1"|s=="SWI12")*(beta[1]+beta[2]*x1+2*x2+x3+beta[3]*x3^2+u)+
  #   1*(s=="ACO")*(1+x1+2*x2+0.2*x2^2+x3+u)
  
  y_1 <- 1*(s=="ANT")*(1+x1+x2+u)+
    1*(s=="CO1")*(1-x1+x2+u)+
    1*(s=="CO2")*(1-x1+x2+beta[1]*(cos(pi*x1)+cos(pi*x2))+beta[2]*(x1+x2)+u)+
    1*(s=="RCO")*(1-x1+x2+u)+
    1*(s=="ECO")*(1-x1+x2+u)+epsilon
  
  y <- (1-d)*y_0+(d)*y_1
  
  ite <- y_1-y_0
  
  dat <- data.frame(x1,x2,d,z2,z1,y,y_0, y_1,s,u,ite)
  return(dat)
}

dataset <- dat_gen(n = 5000)
table(dataset$s)/5000



nptest <- function(dataset,V=2)
{
  ####
  # V = 2
  ####
  index <- createFolds(1:nrow(dataset),V)
  
  psedo_outcome <- c()
  x1_vec <- c()
  x2_vec <- c()
  
  N <- nrow(dataset)
  for(v in 1:V)
  {
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    
    # Here we need to model the joint distribution of (Z_1,Z_2)
    # Model for P(Z_1,Z_2|X) and P(Z_2|X)
    sl_z1 <- SuperLearner(Y = trainingset$z1,
                            X = trainingset[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    trainingset_z1_1 <- trainingset[trainingset$z1==1,]
    trainingset_z1_0 <- trainingset[trainingset$z1==0,]
    
    sl_z2_z1_1 <- SuperLearner(Y = trainingset_z1_1$z2,
                          X = trainingset_z1_1[,c("x1","x2")],
                          family = binomial(),
                          SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                          cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_z2_z1_0 <- SuperLearner(Y = trainingset_z1_1$z2,
                               X = trainingset_z1_1[,c("x1","x2")],
                               family = binomial(),
                               SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                               cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    ##### Prediction for P(Z1,Z2|X)
    z1_1_pred <- predict(sl_z1, newdata = predictset, onlySL = T)$pred
    z1_0_pred <- 1-z1_1_pred
    
    z2_1_z1_1_pred <- predict(sl_z2_z1_1, newdata = predictset, onlySL = T)$pred
    z2_1_z1_0_pred <- predict(sl_z2_z1_0, newdata = predictset, onlySL = T)$pred
    
    z2_z1_1_1_pred <- z1_1_pred*z2_1_z1_1_pred
    z2_z1_1_0_pred <- z1_0_pred*z2_1_z1_0_pred
    z2_z1_0_1_pred <- z1_1_pred*(1-z2_1_z1_1_pred)
    z2_z1_0_0_pred <- z1_0_pred*(1-z2_1_z1_0_pred)
    
    z2_1_pred <- z2_z1_1_1_pred + z2_z1_1_0_pred
    z2_0_pred <- z2_z1_0_1_pred + z2_z1_0_0_pred
    # Model for P(D|Z1, X=x) and P(D|Z2, X=x)
    # It should be noted that P(D|Z1, X=x) and P(D|Z2, X=x) are variationally dependent, therefore we need to model this part
    # in a more careful way
    trainingset_CO1_1 <- trainingset[trainingset$z1==1,]
    trainingset_CO1_0 <- trainingset[trainingset$z1==0,]
    trainingset_CO2_1 <- trainingset[trainingset$z2==1,]
    trainingset_CO2_0 <- trainingset[trainingset$z2==0,]

    sl_d_CO1_1 <- SuperLearner(Y = trainingset_CO1_1$d,
                            X = trainingset_CO1_1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_d_CO1_0 <- SuperLearner(Y = trainingset_CO1_0$d,
                              X = trainingset_CO1_0[,c("x1","x2")],
                              family = binomial(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_d_CO2_1 <- SuperLearner(Y = trainingset_CO2_1$d,
                              X = trainingset_CO2_1[,c("x1","x2")],
                              family = binomial(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_d_CO2_0 <- SuperLearner(Y = trainingset_CO2_0$d,
                              X = trainingset_CO2_0[,c("x1","x2")],
                              family = binomial(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(Y|Z,X=x)
    sl_y_CO1_1 <- SuperLearner(Y = trainingset_CO1_1$y,
                              X = trainingset_CO1_1[,c("x1","x2")],
                              family = gaussian(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_y_CO1_0 <- SuperLearner(Y = trainingset_CO1_0$y,
                              X = trainingset_CO1_0[,c("x1","x2")],
                              family = gaussian(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_y_CO2_1 <- SuperLearner(Y = trainingset_CO2_1$y,
                              X = trainingset_CO2_1[,c("x1","x2")],
                              family = gaussian(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    sl_y_CO2_0 <- SuperLearner(Y = trainingset_CO2_0$y,
                              X = trainingset_CO2_0[,c("x1","x2")],
                              family = gaussian(),
                              SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                              cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    
    # prediction for P(D=1|Z1,X), P(D=1|Z2,X)
    d_CO1_1_pred <- predict(sl_d_CO1_1, newdata = predictset, onlySL = T)$pred
    d_CO1_0_pred <- predict(sl_d_CO1_0, newdata = predictset, onlySL = T)$pred
    d_CO2_1_pred <- predict(sl_d_CO2_1, newdata = predictset, onlySL = T)$pred
    d_CO2_0_pred <- predict(sl_d_CO2_0, newdata = predictset, onlySL = T)$pred
    
    
    # prediction for P(Y|Z1,X), P(Y|Z2,X)
    y_CO1_1_pred <- predict(sl_y_CO1_1, newdata = predictset, onlySL = T)$pred
    y_CO1_0_pred <- predict(sl_y_CO1_0, newdata = predictset, onlySL = T)$pred
    y_CO2_1_pred <- predict(sl_y_CO2_1, newdata = predictset, onlySL = T)$pred
    y_CO2_0_pred <- predict(sl_y_CO2_0, newdata = predictset, onlySL = T)$pred
    
    ####################################
    psedo_outcome_CO1 <- 1/(d_CO1_1_pred-d_CO1_0_pred)*(1*(predictset$z1==1)/(z1_1_pred)*(predictset$y-y_CO1_1_pred)-
                                                       1*(predictset$z1==0)/(z1_0_pred)*(predictset$y-y_CO1_0_pred))-
      (y_CO1_1_pred-y_CO1_0_pred)/(d_CO1_1_pred-d_CO1_0_pred)^2*(1*(predictset$z1==1)/(z1_1_pred)*(predictset$d-d_CO1_1_pred)-
                                                               1*(predictset$z1==0)/(z1_0_pred)*(predictset$d-d_CO1_0_pred))+
      (y_CO1_1_pred-y_CO1_0_pred)/(d_CO1_1_pred-d_CO1_0_pred)
    
    psedo_outcome_CO2 <- 1/(d_CO2_1_pred-d_CO2_0_pred)*(1*(predictset$z2==1)/(z2_1_pred)*(predictset$y-y_CO2_1_pred)-
                                                       1*(predictset$z2==0)/(z2_0_pred)*(predictset$y-y_CO2_0_pred))-
      (y_CO2_1_pred-y_CO2_0_pred)/(d_CO2_1_pred-d_CO2_0_pred)^2*(1*(predictset$z2==1)/(z2_1_pred)*(predictset$d-d_CO2_1_pred)-
                                                               1*(predictset$z2==0)/(z2_0_pred)*(predictset$d-d_CO2_0_pred))+
      (y_CO2_1_pred-y_CO2_0_pred)/(d_CO2_1_pred-d_CO2_0_pred)
    
    
    psedo_outcome_diff <- psedo_outcome_CO1-psedo_outcome_CO2
    
    psedo_outcome <- c(psedo_outcome, psedo_outcome_diff)
    x1_vec <- c(x1_vec, predictset$x1)
    x2_vec <- c(x2_vec, predictset$x2)
  }
  return(list(psedo_outcome = psedo_outcome, x_vec = cbind(x1_vec,x2_vec)))
}




beta_list <- list(c(0,0),c(0.5,0),c(0.5,0.3),c(0.5,0.5))
sample_size_list <- c(1000,2000,3000,5000)
simu_res = NULL
for (n in c(sample_size_list)){
  for (i in c(1:length(beta_list))) {
    dt = dat_gen(n, beta = beta_list[[i]])
    psedo_dat <- nptest(dt)
    W = gp_test(psedo_outcome = psedo_dat$psedo_outcome,
                covariates = psedo_dat$x_vec, k_vec <- ceiling(((n))^{1/5}/(log10(n))))
    W_wald <- W$stat_Wald
    W_series <- W$stat_series
    W_Sn <- W$S_vec
    W_M <- W$M_list
    
    simu_res <- rbind(simu_res,c(test_summary_wald(W_wald,3),n,i,0))
    for (k in c(1:length(W_series))) {
      simu_res <- rbind(simu_res,c(test_summary_series(W_series[k]),n,i,(k)))
      simu_res <- rbind(simu_res,c(test_summary_series2(W_Sn[k],W_M[[k]]),n,i,(10+k)))
    }
    
   
  }
}

FileSave <- paste0("/home/bzhang3/Wang_Rui/np_test/output2/sim_",touse,".csv")

write.table(simu_res, file=FileSave, row.names = FALSE)
