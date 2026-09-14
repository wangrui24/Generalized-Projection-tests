

library("MASS")
library("SuperLearner")
library("caret")
library("fda")
library("np")
library("tidyverse")


library('survminer')
library('survival')
# library('rigr')
library('corrplot')
library("ggpubr")

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


test_a0_1 <-  function(dat,x_con_name,x_bin_name,a_name,y_name,s_name,K=3,V=2){
  
  dataset <- dat[,c(x_con_name,x_bin_name,a_name,y_name,s_name)]
  colnames(dataset) <- c(x_con_name,x_bin_name,"a","y","s")
  dataset[, x_con_name] <- lapply(dataset[, x_con_name, drop = FALSE], function(col) {
    col / max(abs(col), na.rm = TRUE)
  })  # scale the x variables
  
  dataset <- as.data.frame(dataset)
  
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
    
    x_vec <- rbind(x_vec,predictset[,c(x_con_name,x_bin_name)])
    
    trainingset_s0_a0 <- trainingset[which((trainingset$s==0)&(trainingset$a==0)),]
    trainingset_s1_a0 <- trainingset[which((trainingset$s==1)&(trainingset$a==0)),]
    
    
    
    
    sl_y_s0_a0 <- SuperLearner(Y = trainingset_s0_a0$y,
                               X = trainingset_s0_a0[,c(x_con_name,x_bin_name)],
                               family = binomial(),
                               SL.library = c("SL.randomForest"),
                               cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_y_s0_a0 <- predict(sl_y_s0_a0, newdata = predictset, onlySL = T)$pred
    
    sl_y_s1_a0 <- SuperLearner(Y = trainingset_s1_a0$y,
                               X = trainingset_s1_a0[,c(x_con_name,x_bin_name)],
                               family = binomial(),
                               SL.library = c("SL.randomForest"),
                               cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_y_s1_a0 <- predict(sl_y_s1_a0, newdata = predictset, onlySL = T)$pred
    
    sl_s1 <- SuperLearner(Y = trainingset$s,
                          X = trainingset[,c(x_con_name,x_bin_name)],
                          family = binomial(),
                          SL.library = c("SL.glm", "SL.randomForest"),
                          cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    pred_s1 <- predict(sl_s1, newdata = predictset, onlySL = T)$pred
    pred_s0 <- 1-pred_s1
    
    
    
    pred_a0_s1 <- 1 
    
    
    pred_a0_s0 <- 1
    
    pred_a0s1 <- pred_a0_s1*pred_s1
    pred_a0s0 <- pred_a0_s0*pred_s0
    
    #### Calculate psedo-outcome
    psedo_a0 <- (predictset$a0s1/pred_a0s1*(predictset$y-pred_y_s1_a0)+pred_y_s1_a0)-
      (predictset$a0s0/pred_a0s0*(predictset$y-pred_y_s0_a0)+pred_y_s0_a0)
    psedo_a0_vec <- c(psedo_a0_vec, psedo_a0)
    print(which(is.na(psedo_a0)))
  }
  
  
  #### Test 0, using the orginal covariates 
  x_vec <- as.matrix(x_vec)
  m0 <- lm(psedo_a0_vec~x_vec)
  coef_0 <- coef(m0)
  pred_0 <- predict(m0)
  
  b_vec <- cbind(1,x_vec)
  SXX <- t(b_vec)%*%b_vec/N
  meat <- t(b_vec)%*%diag((psedo_a0_vec-pred_0)^2)%*%b_vec/N
  p_0 <- ncol(b_vec)
  
  W_wald <-  N*t(coef_0)%*%solve(solve(SXX)%*%meat%*%solve(SXX))%*%coef_0
  #### Test using the series
  # k <- ceiling(((nrow(dataset)))^{1/3}/(log10(nrow(dataset))))
  W_vec <- c()
  
  bs_basis <- NULL
  
  for (i in c(1:length(x_con_name))) {
    bs_i_basis <- legendre_orthonormal_matrix(x_vec[,i],K)[,-1]
    bs_basis <- cbind(bs_basis,bs_i_basis)
  }
  bs_basis <- cbind(1/sqrt(2), bs_basis)
  
  bs_basis_bin <- x_vec[,c((length(x_con_name)+1):ncol(x_vec))]
  
  bs_basis <- cbind(bs_basis,bs_basis_bin)
  
  sn <- (t(psedo_a0_vec)%*%bs_basis/N)%*%t(t(psedo_a0_vec)%*%bs_basis/N)
  M <-  t(bs_basis)%*%diag(psedo_a0_vec^2)%*%(bs_basis)/N
  
  an <- sum(diag(M))
  bn <- norm(M, type = "F")
  W_series <- (N*sn-an)/(sqrt(2)*bn)

  

  return(list(stat_Wald = W_wald, stat_series = W_series))
}



test_summary_wald <- function(W,p){
  t <- 1*(c(W)>qchisq(seq(0,1,0.05),p))
  return(c(t))
}


test_summary_series <- function(W){
  t <- 1*(c(W)>qnorm(seq(0,1,0.05)))
  return(c(t))
}




data <- read.csv("covail_data_processed_20240313.csv")
p_values_wald <- c()
p_values_series <- c()


################### 
################### KM curve
################### 


#### Moderna

data_km <- data[data$arm %in% c(1,2,4,5,6),c("COVIDtimeD22toD181","COVIDIndD22toD181","arm")]

data_km <- data_km[complete.cases(data_use),]
data_km <- data_km %>% 
  mutate(s = case_when(
    arm %in% c(1) ~ "Moderna (prototype)",
    arm %in% c(2,4,5,6) ~ "Moderna (Omicron)"
  ), s = fct_relevel(s, "Moderna (prototype)", "Moderna (Omicron)"))

table(data_km$s)

# KM curves
model_fit <- survfit(Surv((COVIDtimeD22toD181), COVIDIndD22toD181) ~ s, data = data_km)

km1 <- ggsurvplot(
  model_fit,
  data = data_km,
  fun = "event",
  size = 1.2,                 # change line size
  palette =
    c("#FF8C00", "#0066FF", "#BA55D3", "#009999"),# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = FALSE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.title = "Vaccine",
  legend.labs =
    c("Moderna (prototype)","Moderna (Omicron)"),     # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  censor = FALSE,
  title = "Moderna",
  xlab = "\n days since baseline",
  ylab = "Cumulative incidence of COVID-19 \n",
  ylim = c(0, 0.4),
  legend = "bottom"
)
km1
ggsave("KM1.pdf", width = 5, height = 4)

#### Pfizer

data_km <- data[data$arm %in% c(7,8,9,12),c("COVIDtimeD22toD181","COVIDIndD22toD181","arm")]

data_km <- data_km[complete.cases(data_use),]
data_km <- data_km %>% 
  mutate(s = case_when(
    arm %in% c(7) ~ "Pfizer (prototype)",
    arm %in% c(8,9,12) ~ "Pfizer (Omicron)"
  ), s = fct_relevel(s, "Pfizer (prototype)", "Pfizer (Omicron)"))

table(data_km$s)

# KM curves
model_fit <- survfit(Surv((COVIDtimeD22toD181), COVIDIndD22toD181) ~ s, data = data_km)

km2 <- ggsurvplot(
  model_fit,
  data = data_km,
  fun = "event",
  size = 1.2,                 # change line size
  palette =
    c("#BA55D3", "#009999"),# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = FALSE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.title = "Vaccine",
  legend.labs =
    c("Pfizer (prototype)", "Pfizer (Omicron)"),     # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  censor = FALSE,
  title = "Pfizer",
  xlab = "\n days since baseline",
  ylab = "Cumulative incidence of COVID-19 \n",
  ylim = c(0, 0.4),
  legend = "bottom"
)
km2



ggsave("KM2.pdf", width = 5, height = 4)



######## ALL
data_km <- data[data$arm %in% c(1,2,4,5,6,7,8,9,12),c("COVIDtimeD22toD181","COVIDIndD22toD181","arm")]

data_km <- data_km[complete.cases(data_use),]
data_km <- data_km %>% 
  mutate(s = case_when(
    arm %in% c(1) ~ "Moderna (prototype)",
    arm %in% c(2,4,5,6) ~ "Moderna (Omicron)",
    arm %in% c(7) ~ "Pfizer (prototype)",
    arm %in% c(8,9,12) ~ "Pfizer (Omicron)"
  ), s = fct_relevel(s, "Moderna (prototype)", "Moderna (Omicron)",
                     "Pfizer (prototype)", "Pfizer (Omicron)"))

table(data_km$s)

# KM curves
model_fit <- survfit(Surv((COVIDtimeD22toD181), COVIDIndD22toD181) ~ s, data = data_km)

km_moderna <- ggsurvplot(
  model_fit,
  data = data_km,
  fun = "event",
  size = 1.2,                 # change line size
  palette =
    c("#FF8C00", "#0066FF", "#BA55D3", "#009999"),# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = FALSE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.title = "Vaccine platform",
  legend.labs =
    c("Moderna (prototype)","Moderna (Omicron)","Pfizer (prototype)" ,"Pfizer (Omicron)"),     # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  censor = FALSE,
  xlab = "\n days since baseline",
  ylab = "Cumulative incidence of COVID-19 \n",
  legend = "right"
)
km1
ggsave("KM.pdf", width = 5, height = 4)



set.seed("1234")
################### 
################### Analysis 1
################### 

data_use <- data[data$arm %in% c(2,4,5,6,8,9,12),c("naive","risk_score","FOIstandardized","Bpseudoneutid50_BA.1",
                                        "arm","COVIDIndD22toD181")]

# FOI standardized between stage

data_use <- data_use[complete.cases(data_use),]
data_use <- data_use %>% 
  mutate(s = case_when(
    arm %in% c(2,4,5,6) ~ 0,
    arm %in% c(8,9,12) ~ 1
  )) %>% mutate(
    a = 0
  )


dat = data_use
x_con_name = c("risk_score","FOIstandardized","Day15pseudoneutid50_BA.1")
x_bin_name = c("naive")
a_name = c("a")
y_name = c("COVIDIndD22toD181")
s_name = c("s")
K = 3
V = 3

res1 <- test_a0_1(dat = data_use, x_con_name = c("risk_score","FOIstandardized","Bpseudoneutid50_BA.1"),
          x_bin_name = c("naive"), a_name = c("a"),y_name = c("COVIDIndD22toD181"),
          s_name = c("s"), K = 3, V = 2)

qchisq(0.96,5)
qnorm(0.95)


p_values_wald <- c(p_values_wald, (1-pchisq(res1$stat_Wald, df = 5)))
p_values_series <- c(p_values_series, (1-pnorm(res1$stat_series)))

################### 
################### Analysis 2
################### 

data_use <- data[data$arm %in% c(2,4,5,6,1),c("naive","risk_score","FOIstandardized","Bpseudoneutid50_BA.1",
                                                   "Bpseudoneutid50_BA.1",
                                                   "arm","COVIDIndD22toD181")]

# FOI standardized between stage

data_use <- data_use[complete.cases(data_use),]
data_use <- data_use %>% 
  mutate(s = case_when(
    arm %in% c(2,4,5,6) ~ 0,
    arm %in% c(1) ~ 1
  )) %>% mutate(
    a = 0
  )






res2 <- test_a0_1(dat = data_use, x_con_name = c("risk_score","FOIstandardized","Bpseudoneutid50_BA.1"),
          x_bin_name = c("naive"), a_name = c("a"),y_name = c("COVIDIndD22toD181"),
          s_name = c("s"), K = 3, V =2)

p_values_wald <- c(p_values_wald, (1-pchisq(res2$stat_Wald, df = 5)))
p_values_series <- c(p_values_series, (1-pnorm(res2$stat_series)))



################### 
################### Analysis 3
################### 

data_use <- data[data$arm %in% c(8,9,12,7),c("naive","risk_score","FOIstandardized","Bpseudoneutid50_BA.1",
                                              "Bpseudoneutid50_BA.1",
                                              "arm","COVIDIndD22toD181")]

# FOI standardized between stage

data_use <- data_use[complete.cases(data_use),]
data_use <- data_use %>% 
  mutate(s = case_when(
    arm %in% c(8,9,12) ~ 0,
    arm %in% c(7) ~ 1
  )) %>% mutate(
    a = 0
  )



res3 <- test_a0_1(dat = data_use, x_con_name = c("risk_score","FOIstandardized","Bpseudoneutid50_BA.1"),
          x_bin_name = c("naive"), a_name = c("a"),y_name = c("COVIDIndD22toD181"),
          s_name = c("s"), K = 3, V = 2)

qchisq(0.96,5)
qnorm(0.95)

p_values_wald <- c(p_values_wald, (1-pchisq(res3$stat_Wald, df = 5)))
p_values_series <- c(p_values_series, (1-pnorm(res3$stat_series)))


p_values_wald
p_values_series
# 
# 
# x_con_name <- c("risk_score","FOIstandardized","Day15pseudoneutid50_BA.1")
# x_bin_name <- c("Age65C","naive")
# a_name <- c("a")
# y_name <- c("COVIDIndD22toD181")
# s_name <- c("s")
# K=3
# V=2
# 
# test_a0_1(dataset = data_use, x_con_name = c("risk_score","FOIstandardized"),
#           x_bin_name = c("Age65C","Sex"), a_name = c("a"),y_name = c("COVIDIndD22toD181"),
#           s_name = c("s"), K = 3, V = 5)
# 
# qchisq(0.96,5)
# qnorm(0.95)