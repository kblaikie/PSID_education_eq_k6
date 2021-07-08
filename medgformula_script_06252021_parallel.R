#### Mediational g-formula script - marginal models - same timepoint Y #####
#Author - Kieran Blaikie
#Date - 25 June 2021

#Changes - Compared to medgformula_script_currY_05052021, this script: 
#          1 - Computes the mediational g-formula across 10 imputed datasets for each intervention length
#          2 - Runs all bootstrap repetitions in parallel, saving computational time
#          3 - Performs Rubin's Rules to obtain pooled mean effect estimates, SE, 95%CI

#Overview - This script creates 3-7 time-point marginal models with 95% CI
#         - K-6 is treated as binary, with the following assumed DAG
#           Base C > HS > rep[EQ > SRH > MAR > K-6 ] > ... > T
#         - Baseline C 1-9: 1) year, 2) age, 3) sex, 4) race, 5) nativity, 
#                           6) disability, 7) parental wealth, 8) region, 9) occupation

#Loading libraries
library(doParallel)
library(tidyverse)

#Setting directory
directory <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"

#### Creating mediational g-formula functions ####
medgf_3 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction$fit, prediction$residual.scale)
    }
    return(noisy)
  }
  
  #### 3 time-point marginal model #### 
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                       C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                   A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a1     <- noisy_est("continuous", pM_0a1)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a1     <- noisy_est("continuous", pM_1a1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a1     <- noisy_est("continuous", pM_2a1)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1, M_1 = pM_1a1, M_2 = pM_2a1)
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a0     <- noisy_est("continuous", pM_0a0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a0     <- noisy_est("continuous", pM_1a0)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a0     <- noisy_est("continuous", pM_2a0)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0, M_1 = pM_1a0, M_2 = pM_2a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    for (a in 0:1) {
      #G changes permuted M as K moves from 0-30 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
      }
      
      #Q(a,a) or Q(a,a*)
      pM_0a1g     <- G$M_0
      pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
      pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
      pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a1g     <- noisy_est("binary", pY_0a1g)
      
      pM_1a1g     <- G$M_1
      pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
      pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
      pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a1g     <- noisy_est("binary", pY_1a1g)
      
      pM_2a1g     <- G$M_2
      pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
      pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
      
      #Q(a*,a) or Q(a*,a*)
      pM_0a0g     <- G$M_0
      pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
      pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
      pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a0g     <- noisy_est("binary", pY_0a0g)
      
      pM_1a0g     <- G$M_1
      pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
      pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
      pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a0g     <- noisy_est("binary", pY_1a0g)
      
      pM_2a0g     <- G$M_2
      pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
      pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
      
      #3b - Simulate outcome for each i and save mean outcome
      pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a1g     <- noisy_est("binary", pY_2a1g)
      pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a0g     <- noisy_est("binary", pY_2a0g)
      
      MC$pY_2a1g <- pY_2a1g
      MC$pY_2a0g <- pY_2a0g
      
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_2a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_2a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_2a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_2a0g)
      }
    }
  }
  
  #3c - Calculating causal effects
  outcomes$TE <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #3e - Taking the mean of each effect estimate
  results <- data.frame(rep    = boot,
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),      TE_m   = mean(outcomes$TE_m),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
}
medgf_4 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction$fit, prediction$residual.scale)
    }
    return(noisy)
  }
  
  #### 4 time-point marginal model #### 
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                       C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                         data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                   A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                          data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a1     <- noisy_est("continuous", pM_0a1)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a1     <- noisy_est("continuous", pM_1a1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a1     <- noisy_est("continuous", pM_2a1)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a1     <- noisy_est("continuous", pM_3a1)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1, M_1 = pM_1a1, M_2 = pM_2a1, M_3 = pM_3a1)
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a0     <- noisy_est("continuous", pM_0a0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a0     <- noisy_est("continuous", pM_1a0)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a0     <- noisy_est("continuous", pM_2a0)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a0     <- noisy_est("continuous", pM_3a0)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,  M_1 = pM_1a0, M_2 = pM_2a0, M_3 = pM_3a0)
  
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    for (a in 0:1) {
      #G changes permuted M as K moves from 0-30 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
      }
      
      #Q(a,a) or Q(a,a*)
      pM_0a1g     <- G$M_0
      pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
      pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
      pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a1g     <- noisy_est("binary", pY_0a1g)
      
      pM_1a1g     <- G$M_1
      pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
      pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
      pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a1g     <- noisy_est("binary", pY_1a1g)
      
      pM_2a1g     <- G$M_2
      pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
      pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
      pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a1g     <- noisy_est("binary", pY_2a1g)
      
      pM_3a1g     <- G$M_3
      pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
      pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
      
      #Q(a*,a) or Q(a*,a*)
      pM_0a0g     <- G$M_0
      pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
      pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
      pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a0g     <- noisy_est("binary", pY_0a0g)
      
      pM_1a0g     <- G$M_1
      pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
      pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
      pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a0g     <- noisy_est("binary", pY_1a0g)
      
      pM_2a0g     <- G$M_2
      pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
      pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
      pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a0g     <- noisy_est("binary", pY_2a0g)
      
      pM_3a0g     <- G$M_3
      pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
      pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
      
      #3b - Simulate outcome for each i and save mean outcome
      pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a1g     <- noisy_est("binary", pY_3a1g)
      pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a0g     <- noisy_est("binary", pY_3a0g)
      
      MC$pY_3a1g <- pY_3a1g
      MC$pY_3a0g <- pY_3a0g
      
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_3a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_3a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_3a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_3a0g)
      }
    }
  }
  
  #3c - Calculating causal effects
  outcomes$TE <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #3e - Taking the mean of each effect estimate
  results <- data.frame(rep    = boot,
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),      TE_m   = mean(outcomes$TE_m),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
} 
medgf_5 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction$fit, prediction$residual.scale)
    }
    return(noisy)
  }
  
  #### 5 time-point marginal model #### 
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                       C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                         data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                   A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                          data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                       Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                 M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~         TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_4   <- glm(Y_4   ~ TV2_4 + TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a1     <- noisy_est("continuous", pM_0a1)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a1     <- noisy_est("continuous", pM_1a1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a1     <- noisy_est("continuous", pM_2a1)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a1     <- noisy_est("continuous", pM_3a1)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a1     <- noisy_est("continuous", pM_4a1)
  pTV1_4a1   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a1   <- noisy_est("binary", pTV1_4a1)
  pTV2_4a1   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a1   <- noisy_est("binary", pTV2_4a1)
  pY_4a1     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a1     <- noisy_est("binary", pY_4a1)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1, M_1 = pM_1a1, M_2 = pM_2a1, M_3 = pM_3a1, M_4 = pM_4a1)
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a0     <- noisy_est("continuous", pM_0a0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a0     <- noisy_est("continuous", pM_1a0)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a0     <- noisy_est("continuous", pM_2a0)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a0     <- noisy_est("continuous", pM_3a0)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a0     <- noisy_est("continuous", pM_4a0)
  pTV1_4a0   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a0   <- noisy_est("binary", pTV1_4a0)
  pTV2_4a0   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a0   <- noisy_est("binary", pTV2_4a0)
  pY_4a0     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a0     <- noisy_est("binary", pY_4a0)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  #   - Below is the predicted mediator for each i under A=0 at T=0 - T=4 (M_0 to M_4)
  #   - If predicted M for T=0-4 are taken from the same prediction, below 30 rows are selected for Step 3 K repeats
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0, M_1 = pM_1a0, M_2 = pM_2a0, M_3 = pM_3a0, M_4 = pM_4a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    for (a in 0:1) {
      #G changes permuted M as K moves from 0-30 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
      }
      
      #Q(a,a) or Q(a,a*)
      pM_0a1g     <- G$M_0
      pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
      pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
      pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a1g     <- noisy_est("binary", pY_0a1g)
      
      pM_1a1g     <- G$M_1
      pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
      pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
      pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a1g     <- noisy_est("binary", pY_1a1g)
      
      pM_2a1g     <- G$M_2
      pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
      pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
      pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a1g     <- noisy_est("binary", pY_2a1g)
      
      pM_3a1g     <- G$M_3
      pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
      pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
      pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a1g     <- noisy_est("binary", pY_3a1g)
      
      pM_4a1g     <- G$M_4
      pTV1_4a1g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a1g   <- noisy_est("binary", pTV1_4a1g)
      pTV2_4a1g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a1g   <- noisy_est("binary", pTV2_4a1g)
      
      #Q(a*,a) or Q(a*,a*)
      pM_0a0g     <- G$M_0
      pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
      pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
      pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a0g     <- noisy_est("binary", pY_0a0g)
      
      pM_1a0g     <- G$M_1
      pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
      pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
      pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a0g     <- noisy_est("binary", pY_1a0g)
      
      pM_2a0g     <- G$M_2
      pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
      pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
      pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a0g     <- noisy_est("binary", pY_2a0g)
      
      pM_3a0g     <- G$M_3
      pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
      pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
      pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a0g     <- noisy_est("binary", pY_3a0g)
      
      pM_4a0g     <- G$M_4
      pTV1_4a0g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a0g   <- noisy_est("binary", pTV1_4a0g)
      pTV2_4a0g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a0g   <- noisy_est("binary", pTV2_4a0g)
      
      #3b - Simulate outcome for each i and save mean outcome
      pY_4a1g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a1g     <- noisy_est("binary", pY_4a1g)
      pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a0g     <- noisy_est("binary", pY_4a0g)
      
      MC$pY_4a1g <- pY_4a1g
      MC$pY_4a0g <- pY_4a0g
      
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_4a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_4a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_4a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_4a0g)
      }
    }
  }
  
  #3c - Calculating causal effects
  outcomes$TE <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #3e - Taking the mean of each effect estimate
  results <- data.frame(rep    = boot,
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),      TE_m   = mean(outcomes$TE_m),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
}
medgf_6 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction$fit, prediction$residual.scale)
    }
    return(noisy)
  }
  
  #### 6 time-point marginal model #### 
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                       C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                         data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                   A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                          data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                       Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                 M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~         TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mY_4   <- glm(Y_4   ~ TV2_4 + TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mM_5   <- glm(M_5   ~                       Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = gaussian("identity"))
  mTV1_5 <- glm(TV1_5 ~                 M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mTV2_5 <- glm(TV2_5 ~         TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_5   <- glm(Y_5   ~ TV2_5 + TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a1     <- noisy_est("continuous", pM_0a1)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a1     <- noisy_est("continuous", pM_1a1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a1     <- noisy_est("continuous", pM_2a1)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a1     <- noisy_est("continuous", pM_3a1)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a1     <- noisy_est("continuous", pM_4a1)
  pTV1_4a1   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a1   <- noisy_est("binary", pTV1_4a1)
  pTV2_4a1   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a1   <- noisy_est("binary", pTV2_4a1)
  pY_4a1     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a1     <- noisy_est("binary", pY_4a1)
  
  pM_5a1     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_5a1     <- noisy_est("continuous", pM_5a1)
  pTV1_5a1   <- predict(mTV1_5, newdata = data.frame(                                      M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5a1   <- noisy_est("binary", pTV1_5a1)
  pTV2_5a1   <- predict(mTV2_5, newdata = data.frame(                    TV1_5 = pTV1_5a1, M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5a1   <- noisy_est("binary", pTV2_5a1)
  pY_5a1     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_5a1     <- noisy_est("binary", pY_5a1)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1, M_1 = pM_1a1, M_2 = pM_2a1, M_3 = pM_3a1, M_4 = pM_4a1, M_5 = pM_5a1)
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a0     <- noisy_est("continuous", pM_0a0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a0     <- noisy_est("continuous", pM_1a0)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a0     <- noisy_est("continuous", pM_2a0)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a0     <- noisy_est("continuous", pM_3a0)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a0     <- noisy_est("continuous", pM_4a0)
  pTV1_4a0   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a0   <- noisy_est("binary", pTV1_4a0)
  pTV2_4a0   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a0   <- noisy_est("binary", pTV2_4a0)
  pY_4a0     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a0     <- noisy_est("binary", pY_4a0)
  
  pM_5a0     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_5a0     <- noisy_est("continuous", pM_5a0)
  pTV1_5a0   <- predict(mTV1_5, newdata = data.frame(                                      M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5a0   <- noisy_est("binary", pTV1_5a0)
  pTV2_5a0   <- predict(mTV2_5, newdata = data.frame(                    TV1_5 = pTV1_5a0, M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5a0   <- noisy_est("binary", pTV2_5a0)
  pY_5a0     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_5a0     <- noisy_est("binary", pY_5a0)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0, M_1 = pM_1a0, M_2 = pM_2a0, M_3 = pM_3a0, M_4 = pM_4a0, M_5 = pM_5a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    for (a in 0:1) {
      #G changes permuted M as K moves from 0-30 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
      }
      
      #Q(a,a) or Q(a,a*)
      pM_0a1g     <- G$M_0
      pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
      pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
      pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a1g     <- noisy_est("binary", pY_0a1g)
      
      pM_1a1g     <- G$M_1
      pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
      pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
      pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a1g     <- noisy_est("binary", pY_1a1g)
      
      pM_2a1g     <- G$M_2
      pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
      pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
      pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a1g     <- noisy_est("binary", pY_2a1g)
      
      pM_3a1g     <- G$M_3
      pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
      pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
      pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a1g     <- noisy_est("binary", pY_3a1g)
      
      pM_4a1g     <- G$M_4
      pTV1_4a1g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a1g   <- noisy_est("binary", pTV1_4a1g)
      pTV2_4a1g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a1g   <- noisy_est("binary", pTV2_4a1g)
      pY_4a1g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a1g     <- noisy_est("binary", pY_4a1g)
      
      pM_5a1g     <- G$M_5
      pTV1_5a1g   <- predict(mTV1_5, newdata = data.frame(                                        M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_5a1g   <- noisy_est("binary", pTV1_5a1g)
      pTV2_5a1g   <- predict(mTV2_5, newdata = data.frame(                     TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_5a1g   <- noisy_est("binary", pTV2_5a1g)
      
      #Q(a*,a) or Q(a*,a*)
      pM_0a0g     <- G$M_0
      pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
      pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
      pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a0g     <- noisy_est("binary", pY_0a0g)
      
      pM_1a0g     <- G$M_1
      pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
      pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
      pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a0g     <- noisy_est("binary", pY_1a0g)
      
      pM_2a0g     <- G$M_2
      pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
      pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
      pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a0g     <- noisy_est("binary", pY_2a0g)
      
      pM_3a0g     <- G$M_3
      pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
      pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
      pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a0g     <- noisy_est("binary", pY_3a0g)
      
      pM_4a0g     <- G$M_4
      pTV1_4a0g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a0g   <- noisy_est("binary", pTV1_4a0g)
      pTV2_4a0g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a0g   <- noisy_est("binary", pTV2_4a0g)
      pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a0g     <- noisy_est("binary", pY_4a0g)
      
      pM_5a0g     <- G$M_5
      pTV1_5a0g   <- predict(mTV1_5, newdata = data.frame(                                        M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_5a0g   <- noisy_est("binary", pTV1_5a0g)
      pTV2_5a0g   <- predict(mTV2_5, newdata = data.frame(                     TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_5a0g   <- noisy_est("binary", pTV2_5a0g)
      
      #3b - Simulate outcome for each i and save mean outcome
      pY_5a1g     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_5a1g     <- noisy_est("binary", pY_5a1g)
      pY_5a0g     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_5a0g     <- noisy_est("binary", pY_5a0g)
      
      MC$pY_5a1g <- pY_5a1g
      MC$pY_5a0g <- pY_5a0g
      
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_5a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_5a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_5a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_5a0g)
      }
    }
  }
  
  #3c - Calculating causal effects
  outcomes$TE <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #3e - Taking the mean of each effect estimate
  results <- data.frame(rep    = boot,
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),      TE_m   = mean(outcomes$TE_m),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
} 
medgf_7 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction$fit, prediction$residual.scale)
    }
    return(noisy)
  }
  
  #### 7 time-point marginal model #### 
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                       C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                         data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                   A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                          data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0, data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                       Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                 M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~         TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mY_4   <- glm(Y_4   ~ TV2_4 + TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mM_5   <- glm(M_5   ~                       Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = gaussian("identity"))
  mTV1_5 <- glm(TV1_5 ~                 M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mTV2_5 <- glm(TV2_5 ~         TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mY_5   <- glm(Y_5   ~ TV2_5 + TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mM_6   <- glm(M_6   ~                       Y_5 + TV2_5 + TV1_5 + M_5 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = gaussian("identity"))
  mTV1_6 <- glm(TV1_6 ~                 M_6 + Y_5 + TV2_5 + TV1_5 + M_5 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_6 + C3*M_6 + C4*M_6, data = data_boot, family = binomial("logit"))
  mTV2_6 <- glm(TV2_6 ~         TV1_6 + M_6 + Y_5 + TV2_5 + TV1_5 + M_5 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_6 + C3*M_6 + C4*M_6, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_6   <- glm(Y_6   ~ TV2_6 + TV1_6 + M_6 + Y_5 + TV2_5 + TV1_5 + M_5 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_6 + C3*M_6 + C4*M_6, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a1     <- noisy_est("continuous", pM_0a1)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a1     <- noisy_est("continuous", pM_1a1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a1     <- noisy_est("continuous", pM_2a1)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a1     <- noisy_est("continuous", pM_3a1)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a1     <- noisy_est("continuous", pM_4a1)
  pTV1_4a1   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a1   <- noisy_est("binary", pTV1_4a1)
  pTV2_4a1   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a1   <- noisy_est("binary", pTV2_4a1)
  pY_4a1     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a1     <- noisy_est("binary", pY_4a1)
  
  pM_5a1     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_5a1     <- noisy_est("continuous", pM_5a1)
  pTV1_5a1   <- predict(mTV1_5, newdata = data.frame(                                      M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5a1   <- noisy_est("binary", pTV1_5a1)
  pTV2_5a1   <- predict(mTV2_5, newdata = data.frame(                    TV1_5 = pTV1_5a1, M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5a1   <- noisy_est("binary", pTV2_5a1)
  pY_5a1     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_5a1     <- noisy_est("binary", pY_5a1)
  
  pM_6a1     <- predict(mM_6, newdata = data.frame(                                                      Y_5 = pY_5a1, TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_6a1     <- noisy_est("continuous", pM_6a1)
  pTV1_6a1   <- predict(mTV1_6, newdata = data.frame(                                      M_6 = pM_6a1, Y_5 = pY_5a1, TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_6a1   <- noisy_est("binary", pTV1_6a1)
  pTV2_6a1   <- predict(mTV2_6, newdata = data.frame(                    TV1_6 = pTV1_6a1, M_6 = pM_6a1, Y_5 = pY_5a1, TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_6a1   <- noisy_est("binary", pTV2_6a1)
  pY_6a1     <- predict(mY_6, newdata = data.frame(    TV2_6 = pTV2_6a1, TV1_6 = pTV1_6a1, M_6 = pM_6a1, Y_5 = pY_5a1, TV2_5 = pTV2_5a1, TV1_5 = pTV1_5a1, M_5 = pM_5a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_6a1     <- noisy_est("binary", pY_6a1)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1, M_1 = pM_1a1, M_2 = pM_2a1, M_3 = pM_3a1, M_4 = pM_4a1, M_5 = pM_5a1, M_6 = pM_6a1)
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_0a0     <- noisy_est("continuous", pM_0a0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_1a0     <- noisy_est("continuous", pM_1a0)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_2a0     <- noisy_est("continuous", pM_2a0)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_3a0     <- noisy_est("continuous", pM_3a0)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_4a0     <- noisy_est("continuous", pM_4a0)
  pTV1_4a0   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a0   <- noisy_est("binary", pTV1_4a0)
  pTV2_4a0   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a0   <- noisy_est("binary", pTV2_4a0)
  pY_4a0     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a0     <- noisy_est("binary", pY_4a0)
  
  pM_5a0     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_5a0     <- noisy_est("continuous", pM_5a0)
  pTV1_5a0   <- predict(mTV1_5, newdata = data.frame(                                      M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5a0   <- noisy_est("binary", pTV1_5a0)
  pTV2_5a0   <- predict(mTV2_5, newdata = data.frame(                    TV1_5 = pTV1_5a0, M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5a0   <- noisy_est("binary", pTV2_5a0)
  pY_5a0     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_5a0     <- noisy_est("binary", pY_5a0)
  
  pM_6a0     <- predict(mM_6, newdata = data.frame(                                                      Y_5 = pY_5a0, TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), se.fit = T, interval = "prediction")
  pM_6a0     <- noisy_est("continuous", pM_6a0)
  pTV1_6a0   <- predict(mTV1_6, newdata = data.frame(                                      M_6 = pM_6a0, Y_5 = pY_5a0, TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_6a0   <- noisy_est("binary", pTV1_6a0)
  pTV2_6a0   <- predict(mTV2_6, newdata = data.frame(                    TV1_6 = pTV1_6a0, M_6 = pM_6a0, Y_5 = pY_5a0, TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_6a0   <- noisy_est("binary", pTV2_6a0)
  pY_6a0     <- predict(mY_6, newdata = data.frame(    TV2_6 = pTV2_6a0, TV1_6 = pTV1_6a0, M_6 = pM_6a0, Y_5 = pY_5a0, TV2_5 = pTV2_5a0, TV1_5 = pTV1_5a0, M_5 = pM_5a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_6a0     <- noisy_est("binary", pY_6a0)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0, M_1 = pM_1a0, M_2 = pM_2a0, M_3 = pM_3a0, M_4 = pM_4a0, M_5 = pM_5a0, M_6 = pM_6a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    for (a in 0:1) {
      #G changes permuted M as K moves from 0-30 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
      }
      
      #Q(a,a) or Q(a,a*)
      pM_0a1g     <- G$M_0
      pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
      pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
      pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a1g     <- noisy_est("binary", pY_0a1g)
      
      pM_1a1g     <- G$M_1
      pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
      pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
      pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a1g     <- noisy_est("binary", pY_1a1g)
      
      pM_2a1g     <- G$M_2
      pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
      pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
      pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a1g     <- noisy_est("binary", pY_2a1g)
      
      pM_3a1g     <- G$M_3
      pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
      pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
      pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a1g     <- noisy_est("binary", pY_3a1g)
      
      pM_4a1g     <- G$M_4
      pTV1_4a1g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a1g   <- noisy_est("binary", pTV1_4a1g)
      pTV2_4a1g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a1g   <- noisy_est("binary", pTV2_4a1g)
      pY_4a1g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a1g     <- noisy_est("binary", pY_4a1g)
      
      pM_5a1g     <- G$M_5
      pTV1_5a1g   <- predict(mTV1_5, newdata = data.frame(                                        M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_5a1g   <- noisy_est("binary", pTV1_5a1g)
      pTV2_5a1g   <- predict(mTV2_5, newdata = data.frame(                     TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_5a1g   <- noisy_est("binary", pTV2_5a1g)
      pY_5a1g     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_5a1g     <- noisy_est("binary", pY_5a1g)
      
      pM_6a1g     <- G$M_6
      pTV1_6a1g   <- predict(mTV1_6, newdata = data.frame(                                        M_6 = pM_6a1g, Y_5 = pY_5a1g, TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_6a1g   <- noisy_est("binary", pTV1_6a1g)
      pTV2_6a1g   <- predict(mTV2_6, newdata = data.frame(                     TV1_6 = pTV1_6a1g, M_6 = pM_6a1g, Y_5 = pY_5a1g, TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_6a1g   <- noisy_est("binary", pTV2_6a1g)
      
      #Q(a*,a) or Q(a*,a*)
      pM_0a0g     <- G$M_0
      pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
      pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
      pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_0a0g     <- noisy_est("binary", pY_0a0g)
      
      pM_1a0g     <- G$M_1
      pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
      pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
      pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_1a0g     <- noisy_est("binary", pY_1a0g)
      
      pM_2a0g     <- G$M_2
      pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
      pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
      pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_2a0g     <- noisy_est("binary", pY_2a0g)
      
      pM_3a0g     <- G$M_3
      pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
      pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
      pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_3a0g     <- noisy_est("binary", pY_3a0g)
      
      pM_4a0g     <- G$M_4
      pTV1_4a0g   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_4a0g   <- noisy_est("binary", pTV1_4a0g)
      pTV2_4a0g   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_4a0g   <- noisy_est("binary", pTV2_4a0g)
      pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_4a0g     <- noisy_est("binary", pY_4a0g)
      
      pM_5a0g     <- G$M_5
      pTV1_5a0g   <- predict(mTV1_5, newdata = data.frame(                                        M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_5a0g   <- noisy_est("binary", pTV1_5a0g)
      pTV2_5a0g   <- predict(mTV2_5, newdata = data.frame(                     TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_5a0g   <- noisy_est("binary", pTV2_5a0g)
      pY_5a0g     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_5a0g     <- noisy_est("binary", pY_5a0g)
      
      pM_6a0g     <- G$M_6
      pTV1_6a0g   <- predict(mTV1_6, newdata = data.frame(                                        M_6 = pM_6a0g, Y_5 = pY_5a0g, TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV1_6a0g   <- noisy_est("binary", pTV1_6a0g)
      pTV2_6a0g   <- predict(mTV2_6, newdata = data.frame(                     TV1_6 = pTV1_6a0g, M_6 = pM_6a0g, Y_5 = pY_5a0g, TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pTV2_6a0g   <- noisy_est("binary", pTV2_6a0g)
      
      #3b - Simulate outcome for each i and save mean outcome
      pY_6a1g     <- predict(mY_6, newdata = data.frame(    TV2_6 = pTV2_6a1g, TV1_6 = pTV1_6a1g, M_6 = pM_6a1g, Y_5 = pY_5a1g, TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_6a1g     <- noisy_est("binary", pY_6a1g)
      pY_6a0g     <- predict(mY_6, newdata = data.frame(    TV2_6 = pTV2_6a0g, TV1_6 = pTV1_6a0g, M_6 = pM_6a0g, Y_5 = pY_5a0g, TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
      pY_6a0g     <- noisy_est("binary", pY_6a0g)
      
      MC$pY_6a1g <- pY_6a1g
      MC$pY_6a0g <- pY_6a0g
      
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_6a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_6a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_6a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_6a0g)
      }
    }
  }
  
  #3c - Calculating causal effects
  outcomes$TE <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #3e - Taking the mean of each effect estimate
  results <- data.frame(rep    = boot,
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),      TE_m   = mean(outcomes$TE_m),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
} 

#### Setting seed for replicability ####
set.seed(2021)

#### Performing mediational g-formula for each of 10 multiply imputed datasets for 3-7 period interventions ####
#Running the mediational g-formula with 200 bootstrap resamples using 40 cores takes ~3-7 minutes per dataset
#Parallelising the process speeds up the formula itself, but when combined with multiple imputation and multiple
#intervention periods (each requiring their own function to be supplied, created above), the process still 
#takes several hours to fully complete.
#The below loop performs the mediational g-formula for each imputed dataset (mi) for each intervention length (int)
for (int in 3:7) {
  for (mi in 1:10) {
    start <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " Start time: ", start))
    
    #Loading multiply imputed dataset
    medgf_wide <- read.csv(paste0(directory, "MI_data/r_dataset_06242021_wide_", int, "_mi_", mi, ".csv"))
    
    #Subsetting to necessary variables for analysis
    medgf_wide %>% select(unique_id, base_hs_or_less, base_years_since_2001, 
                          base_age_mean_centred, female, poc_hisp, foreign, base_disabled_work, 
                          parents_poor, base_Region, base_occ_category, 
                          starts_with("srh_vgood_exc"), starts_with("married"),
                          starts_with("eq_score_trim_"),
                          starts_with("k6_bin")) -> dataset
    rm(medgf_wide)
    
    #Simplifying naming in working dataset
    if (int == 3) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", 
                          "TV2_0", "TV2_1", "TV2_2", 
                          "M_0", "M_1", "M_2",  
                          "Y_0", "Y_1", "Y_2")
    }
    if (int == 4) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3",
                          "M_0", "M_1", "M_2", "M_3", 
                          "Y_0", "Y_1", "Y_2", "Y_3")
    }
    if (int == 5) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4",
                          "M_0", "M_1", "M_2", "M_3", "M_4", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4")
    }
    if (int == 6) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", "TV1_5", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4", "TV2_5",
                          "M_0", "M_1", "M_2", "M_3", "M_4", "M_5", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5")
    }
    if (int == 7) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", "TV1_5", "TV1_6", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4", "TV2_5", "TV2_6",
                          "M_0", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6",
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5", "Y_6")
    }
    
    #Performing mediational g-formula in parallel and exporting output
    clusters <- parallel::makeCluster(50)
    doParallel::registerDoParallel(clusters)
    if (int == 3) {
      assign(paste0("int_3_mi_", mi), foreach(boot = 1:200, .combine = rbind) %dopar% medgf_3(dataset))
    }
    if (int == 4) {
      assign(paste0("int_4_mi_", mi), foreach(boot = 1:200, .combine = rbind) %dopar% medgf_4(dataset))
    }
    if (int == 5) {
      assign(paste0("int_5_mi_", mi), foreach(boot = 1:200, .combine = rbind) %dopar% medgf_5(dataset))
    }
    if (int == 6) {
      assign(paste0("int_6_mi_", mi), foreach(boot = 1:200, .combine = rbind) %dopar% medgf_6(dataset))
    }
    if (int == 7) {
      assign(paste0("int_7_mi_", mi), foreach(boot = 1:200, .combine = rbind) %dopar% medgf_7(dataset))
    }
    parallel::stopCluster(clusters)
    write.csv(eval(parse(text = paste0("int_", int, "_mi_", mi))), 
              paste0(directory, "MI_data/boot_results_06252021_", int, "_mi_", mi, ".csv"), row.names = FALSE, na = "")
    end <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " completed. Time elapsed: ", (end - start)))
    print("")
    rm(dataset, start, end, clusters)
  }
}
rm(int, mi)

#### Combining imputation-specific findings to give summary findings ####
#Following Rubin's rules and the MI-Boot process recommended by Schomaker & Heumann (2019),
#the below takes mean estimates and their bootstrap SE across each MI dataset, 
#then follows Rubin's rules to pool these and give summary estimates.
for (int in 3:7) {
  int_results <- data.frame(mi_set = NA,
                            Ya0g0  = NA,  Ya0g1  = NA,  Ya1g0  = NA,  Ya1g1  = NA,
                            TE     = NA,  TE_m   = NA,
                            pNDE   = NA,  tNDE   = NA,  pNIE   = NA,  tNIE   = NA,
                            pNDE_m = NA,  tNDE_m = NA,  pNIE_m = NA,  tNIE_m = NA)
  int_results <- int_results[0, ]
  for (mi in 1:10) {
    data <- eval(parse(text = paste0("int_", int, "_mi_", mi)))
    results <- data.frame(mi_set = mi,
                          Ya0g0  = mean(data$Ya0g0),   Ya0g1  = mean(data$Ya0g1),
                          Ya1g0  = mean(data$Ya1g0),   Ya1g1  = mean(data$Ya1g1),
                          TE     = mean(data$TE),      TE_m   = mean(data$TE_m),
                          pNDE   = mean(data$pNDE),    tNDE   = mean(data$tNDE),
                          pNIE   = mean(data$pNIE),    tNIE   = mean(data$tNIE),
                          pNDE_m = mean(data$pNDE_m),  tNDE_m = mean(data$tNDE_m),
                          pNIE_m = mean(data$pNIE_m),  tNIE_m = mean(data$tNIE_m))
    results_se <- data.frame(mi_set = mi,
                             Ya0g0  = sd(data$Ya0g0),   Ya0g1  = sd(data$Ya0g1),
                             Ya1g0  = sd(data$Ya1g0),   Ya1g1  = sd(data$Ya1g1),
                             TE     = sd(data$TE),      TE_m   = sd(data$TE_m),
                             pNDE   = sd(data$pNDE),    tNDE   = sd(data$tNDE),
                             pNIE   = sd(data$pNIE),    tNIE   = sd(data$tNIE),
                             pNDE_m = sd(data$pNDE_m),  tNDE_m = sd(data$tNDE_m),
                             pNIE_m = sd(data$pNIE_m),  tNIE_m = sd(data$tNIE_m))
    results <- rbind(results, results_se)
    int_results <- rbind(int_results, results)
    rm(results, results_se)
  }
  assign(paste0("int_", int, "_results"), int_results)
  int_point <- int_results[c(1,3,5,7,9,11,13,15,17,19), c(2:15)]
  int_se    <- int_results[c(2,4,6,8,10,12,14,16,18,20), c(2:15)]
  int_point_summary <- data.frame(Ya0g0  = mean(int_point$Ya0g0),   Ya0g1  = mean(int_point$Ya0g1),
                                  Ya1g0  = mean(int_point$Ya1g0),   Ya1g1  = mean(int_point$Ya1g1),
                                  TE     = mean(int_point$TE),      TE_m   = mean(int_point$TE_m),
                                  pNDE   = mean(int_point$pNDE),    tNDE   = mean(int_point$tNDE),
                                  pNIE   = mean(int_point$pNIE),    tNIE   = mean(int_point$tNIE),
                                  pNDE_m = mean(int_point$pNDE_m),  tNDE_m = mean(int_point$tNDE_m),
                                  pNIE_m = mean(int_point$pNIE_m),  tNIE_m = mean(int_point$tNIE_m))
  int_within_var_summary <- data.frame(Ya0g0  = mean((int_se$Ya0g0^2)),   Ya0g1  = mean((int_se$Ya0g1^2)),
                                       Ya1g0  = mean((int_se$Ya1g0^2)),   Ya1g1  = mean((int_se$Ya1g1^2)),
                                       TE     = mean((int_se$TE^2)),      TE_m   = mean((int_se$TE_m^2)),
                                       pNDE   = mean((int_se$pNDE^2)),    tNDE   = mean((int_se$tNDE^2)),
                                       pNIE   = mean((int_se$pNIE^2)),    tNIE   = mean((int_se$tNIE^2)),
                                       pNDE_m = mean((int_se$pNDE_m^2)),  tNDE_m = mean((int_se$tNDE_m^2)),
                                       pNIE_m = mean((int_se$pNIE_m^2)),  tNIE_m = mean((int_se$tNIE_m^2)))
  int_between_var_summary <- data.frame(Ya0g0 = sum((int_point$Ya0g0 - int_point_summary$Ya0g0)^2)/9,
                                        Ya0g1 = sum((int_point$Ya0g1 - int_point_summary$Ya0g1)^2)/9,
                                        Ya1g0 = sum((int_point$Ya1g0 - int_point_summary$Ya1g0)^2)/9,
                                        Ya1g1 = sum((int_point$Ya1g1 - int_point_summary$Ya1g1)^2)/9,
                                        TE   = sum((int_point$TE - int_point_summary$TE)^2)/9,
                                        TE_m = sum((int_point$TE_m - int_point_summary$TE_m)^2)/9,
                                        pNDE = sum((int_point$pNDE - int_point_summary$pNDE)^2)/9,
                                        tNDE = sum((int_point$tNDE - int_point_summary$tNDE)^2)/9,
                                        pNIE = sum((int_point$pNIE - int_point_summary$pNIE)^2)/9,
                                        tNIE = sum((int_point$tNIE - int_point_summary$tNIE)^2)/9,
                                        pNDE_m = sum((int_point$pNDE_m - int_point_summary$pNDE_m)^2)/9,
                                        tNDE_m = sum((int_point$tNDE_m - int_point_summary$tNDE_m)^2)/9,
                                        pNIE_m = sum((int_point$pNIE_m - int_point_summary$pNIE_m)^2)/9,
                                        tNIE_m = sum((int_point$tNIE_m - int_point_summary$tNIE_m)^2)/9)
  int_se_summary <- data.frame(Ya0g0  = sqrt(sum(int_within_var_summary$Ya0g0 + int_between_var_summary$Ya0g0 + (int_between_var_summary$Ya0g0/10))),
                               Ya0g1  = sqrt(sum(int_within_var_summary$Ya0g1 + int_between_var_summary$Ya0g1 + (int_between_var_summary$Ya0g1/10))),
                               Ya1g0  = sqrt(sum(int_within_var_summary$Ya1g0 + int_between_var_summary$Ya1g0 + (int_between_var_summary$Ya1g0/10))),
                               Ya1g1  = sqrt(sum(int_within_var_summary$Ya1g1 + int_between_var_summary$Ya1g1 + (int_between_var_summary$Ya1g1/10))),
                               TE     = sqrt(sum(int_within_var_summary$TE + int_between_var_summary$TE + (int_between_var_summary$TE/10))),
                               TE_m   = sqrt(sum(int_within_var_summary$TE_m + int_between_var_summary$TE_m + (int_between_var_summary$TE_m/10))),
                               pNDE   = sqrt(sum(int_within_var_summary$pNDE + int_between_var_summary$pNDE + (int_between_var_summary$pNDE/10))),
                               tNDE   = sqrt(sum(int_within_var_summary$tNDE + int_between_var_summary$tNDE + (int_between_var_summary$tNDE/10))),
                               pNIE   = sqrt(sum(int_within_var_summary$pNIE + int_between_var_summary$pNIE + (int_between_var_summary$pNIE/10))),
                               tNIE   = sqrt(sum(int_within_var_summary$tNIE + int_between_var_summary$tNIE + (int_between_var_summary$tNIE/10))),
                               pNDE_m = sqrt(sum(int_within_var_summary$pNDE_m + int_between_var_summary$pNDE_m + (int_between_var_summary$pNDE_m/10))),
                               tNDE_m = sqrt(sum(int_within_var_summary$tNDE_m + int_between_var_summary$tNDE_m + (int_between_var_summary$tNDE_m/10))),
                               pNIE_m = sqrt(sum(int_within_var_summary$pNIE_m + int_between_var_summary$pNIE_m + (int_between_var_summary$pNIE_m/10))),
                               tNIE_m = sqrt(sum(int_within_var_summary$tNIE_m + int_between_var_summary$tNIE_m + (int_between_var_summary$tNIE_m/10))))
  int_lci_summary <- data.frame(Ya0g0 = int_point_summary$Ya0g0 - 1.96*int_se_summary$Ya0g0,
                                Ya0g1  = int_point_summary$Ya0g1 - 1.96*int_se_summary$Ya0g1,
                                Ya1g0  = int_point_summary$Ya1g0 - 1.96*int_se_summary$Ya1g0,
                                Ya1g1  = int_point_summary$Ya1g1 - 1.96*int_se_summary$Ya1g1,
                                TE     = int_point_summary$TE - 1.96*int_se_summary$TE,
                                TE_m   = int_point_summary$TE_m - 1.96*int_se_summary$TE_m,
                                pNDE   = int_point_summary$pNDE - 1.96*int_se_summary$pNDE,
                                tNDE   = int_point_summary$tNDE - 1.96*int_se_summary$tNDE,
                                pNIE   = int_point_summary$pNIE - 1.96*int_se_summary$pNIE,
                                tNIE   = int_point_summary$tNIE - 1.96*int_se_summary$tNIE,
                                pNDE_m = int_point_summary$pNDE_m - 1.96*int_se_summary$pNDE_m,
                                tNDE_m = int_point_summary$tNDE_m - 1.96*int_se_summary$tNDE_m,
                                pNIE_m = int_point_summary$pNIE_m - 1.96*int_se_summary$pNIE_m,
                                tNIE_m = int_point_summary$tNIE_m - 1.96*int_se_summary$tNIE_m)
  int_uci_summary <- data.frame(Ya0g0  = int_point_summary$Ya0g0 + 1.96*int_se_summary$Ya0g0,
                                Ya0g1  = int_point_summary$Ya0g1 + 1.96*int_se_summary$Ya0g1,
                                Ya1g0  = int_point_summary$Ya1g0 + 1.96*int_se_summary$Ya1g0,
                                Ya1g1  = int_point_summary$Ya1g1 + 1.96*int_se_summary$Ya1g1,
                                TE     = int_point_summary$TE + 1.96*int_se_summary$TE,
                                TE_m   = int_point_summary$TE_m + 1.96*int_se_summary$TE_m,
                                pNDE   = int_point_summary$pNDE + 1.96*int_se_summary$pNDE,
                                tNDE   = int_point_summary$tNDE + 1.96*int_se_summary$tNDE,
                                pNIE   = int_point_summary$pNIE + 1.96*int_se_summary$pNIE,
                                tNIE   = int_point_summary$tNIE + 1.96*int_se_summary$tNIE,
                                pNDE_m = int_point_summary$pNDE_m + 1.96*int_se_summary$pNDE_m,
                                tNDE_m = int_point_summary$tNDE_m + 1.96*int_se_summary$tNDE_m,
                                pNIE_m = int_point_summary$pNIE_m + 1.96*int_se_summary$pNIE_m,
                                tNIE_m = int_point_summary$tNIE_m + 1.96*int_se_summary$tNIE_m)
  int_summary <- rbind(int_point_summary, int_se_summary, int_lci_summary, int_uci_summary)
  int_summary <- cbind(data.frame(Measure = c("Point", "SE", "LCI", "UCI")), int_summary)
  assign(paste0("int_", int, "_summary"), int_summary)
  rm(int_results, int_point, int_se, 
     int_point_summary, int_within_var_summary, int_between_var_summary, 
     int_se_summary, int_lci_summary, int_uci_summary)
}
rm(int, mi)

#### Creating and saving summary dataset with all intervention period findings ####
summary_all <- data.frame(Period = NA, Point = NA, SE = NA, LCI = NA, UCI = NA)
summary_all <- summary_all[0, ]
for (int in 3:7) {
  data <- eval(parse(text = paste0("int_", int, "_summary")))
  data <- as.data.frame(t(data))
  names(data) <- data[1,]
  data <- data[2:15, ]
  data$Period <- int*2
  data$Measure <- row.names(data)
  data <- data[,c("Period", "Measure", "Point", "SE", "LCI", "UCI")]
  assign(paste0("int_", int, "_summary"), data)
  write.csv(eval(parse(text = paste0("int_", int, "_summary"))), paste0(directory, "MI_data/pooled_summary_06252022_int_", int, ".csv"), row.names = FALSE, na = "")
  summary_all <- rbind(summary_all, data)
  write.csv(summary_all, paste0(directory, "MI_data/pooled_summary_06252022.csv"), row.names = FALSE, na = "")
}
rm(data, int)
