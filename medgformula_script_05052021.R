#### Mediational g-formula script - marginal models - same timepoint Y #####
#Author - Kieran Blaikie
#Date - 5 May 2021

#Changes - Compared to the 05012021.R marginal model file, this script includes
#          marital status as a time-varying confounder (TV2)

#Overview - This script creates 3-6 time-point marginal models with 95% CI
#         - K-6 is treated as binary, with the following assumed DAG
#           Base C > HS > rep[EQ > SRH > MAR > K-6 ] > ... > T
#         - Baseline C 1-9: 1) year, 2) age, 3) sex, 4) race, 5) nativity, 
#                           6) disability, 7) parental wealth, 8) region, 9) occupation

library(tidyverse)

directory <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"
medgf_wide <- read.csv(paste0(directory, "r_dataset_wide_04302021.csv"))

#Subsetting to necessary variables for analysis
medgf_wide %>% select(unique_id, hs_or_less_majority, base_year, 
                      base_age, female, poc_hisp, foreign, base_disabled_work, 
                      parents_poor, base_Region, occ1990dd_broad_0, 
                      starts_with("srh_vgood_exc"), starts_with("married"),
                      starts_with("eq_score_imp_pct_age_trim"),
                      starts_with("k6_impute_mod")) -> medgf_wide

#Simplifies naming in working dataset
#  _N indicates time N
# A = hs_or_less_majority, M = EQ score, Y = K-6
# C1-8 = base year, age, female, poc_hisp, foreign, disabled, parents poor, region
# TV1-4 = employment type, low_eq_occ, SRH, married
dataset <- medgf_wide
names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                    "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", "TV1_5", "TV1_6", "TV1_7", "TV1_8",
                    "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4", "TV2_5", "TV2_6", "TV2_7", "TV2_8",
                    "M_0", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8",
                    "Y_0", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5", "Y_6", "Y_7", "Y_8")

#Function to add appropriate noise to predictions for each i
#   - Via predict(), each i within covariate strata are assigned the same predicted mean
#   - Within each MC sample, the variable distr. roughly equals that observed. Across MC samples, each i would roughly equal their observed value
#   - To reflect the actual variable variance, need to add in noise based on prediction uncertainty
#   - For binary vars, due to mean-variance relationship, we can do this by sampling from a binomial distribution
#   - For cont. vars, need to estimate prediction SE, then add noise to point prediction from N(0, SE) 
#   - Given noisy_est is used here only with MC resamples each 10000 in size, changed 1:NROW(prediction) to 1:10000
noisy_est <- function(type, prediction) {
  if (type == "binary") {
    noisy <- sapply(1:10000, function(x) {rbinom(1,1,prediction[x])})
  } else {
    noisy <- sapply(1:10000, function(x) {rnorm(1,prediction$fit[x], prediction$residual.scale)})
  }
  return(noisy)
}

#Setting seed for replicability
set.seed(2021)

#### 3 time-point marginal model #### 
#Restricting to those with 3 time-points
dataset_3 <- dataset[!is.na(dataset$Y_2), ]

#Creating summary dataset of results for obtaining 95% CI
results <- data.frame(bootstrap  = rep(NA, 30),
                      Ya0g0  = rep(NA, 30),  Ya0g1  = rep(NA, 30),  Ya1g0  = rep(NA, 30),
                      Ya1g1  = rep(NA, 30),  pNDE   = rep(NA, 30),  tNDE   = rep(NA, 30),
                      pNIE   = rep(NA, 30),  tNIE   = rep(NA, 30),  TE     = rep(NA, 30),
                      pNDE_m = rep(NA, 30),  tNDE_m = rep(NA, 30),  pNIE_m = rep(NA, 30),
                      tNIE_m = rep(NA, 30),  TE_m   = rep(NA, 30))

#Creating progress bar for below (very long)
pb <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)

for (boot in 1:30) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Bootstrap number
  results$bootstrap[boot] <- boot
  
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset_3), NROW(dataset_3), replace = T)
  data_boot <- dataset_3[index_boot, ]
  
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
  
  #1b - parametric model for final Y given measured past
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1 + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #Creating a dataset to store the estimated outcomes for each k repetition in 3d
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #2ai - Generate time t=0-(T-1) M & tv-C based on Step 1 model coefficients 
  #2aii - Assign time t A under intervention A=1
  #     - In this scenario if A is assumed to be fixed before all M, tv-C, Y, we can assign a or a* at baseline only
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
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1)
  
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
  #   - Below is the predicted mediator for each i under A=0 at T=0 - T=4 (M_0 to M_4)
  #   - If predicted M for T=0-4 are taken from the same prediction, below 30 rows are selected for Step 3 K repeats
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    #3a/d - The below simulates Y at T=5 K times, saving the mean Y for each K in the above 'outcomes' dataset where A=0/1 and M=Ga based on M|A=1
    #     - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
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
    #Updating progress bar
    setTxtProgressBar(pk,k)
    
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
  #Adding above summary measures to 'results' dataset for 95% CI estimation
  results$Ya0g0[boot] <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot] <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot] <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot] <- mean(outcomes$Ya1g1)
  results$TE[boot]    <- mean(outcomes$TE)
  results$pNDE[boot]  <- mean(outcomes$pNDE)
  results$tNDE[boot]  <- mean(outcomes$tNDE)
  results$pNIE[boot]  <- mean(outcomes$pNIE)
  results$tNIE[boot]  <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Saving bootstrap results repeatedly to prevent crashing
  write.csv(results, paste0(directory, "boot_results_05052021_3_currY.csv"), row.names = FALSE, na = "")
} 
results_3 <- results

rm(pb, boot, index_boot, data_boot, results, 
   mA, mM_0, mM_1, mM_2, mTV1_0, mTV1_1, mTV1_2, mTV2_0, mTV2_1, mTV2_2,
   mY_0, mY_1, mY_2, outcomes, pk, k, index, MC, 
   pM_0a1, pM_1a1, pM_2a1, pTV1_0a1, pTV1_1a1, pTV1_2a1,
   pTV2_0a1, pTV2_1a1, pTV2_2a1, pY_0a1, pY_1a1, pY_2a1, perm_Mt_a1,
   pM_0a0, pM_1a0, pM_2a0, pTV1_0a0, pTV1_1a0, pTV1_2a0,
   pTV2_0a0, pTV2_1a0, pTV2_2a0, pY_0a0, pY_1a0, pY_2a0, perm_Mt_a0,
   a, G, pM_0a1g, pM_1a1g, pM_2a1g, pTV1_0a1g, pTV1_1a1g, pTV1_2a1g,
   pTV2_0a1g, pTV2_1a1g, pTV2_2a1g, pY_0a1g, pY_1a1g, pY_2a1g,
   pM_0a0g, pM_1a0g, pM_2a0g, pTV1_0a0g, pTV1_1a0g, pTV1_2a0g,
   pTV2_0a0g, pTV2_1a0g, pTV2_2a0g,pY_0a0g, pY_1a0g, pY_2a0g)

#### 4 time-point marginal model #### 
#Restricting to those with 3 time-points
dataset_4 <- dataset[!is.na(dataset$Y_3), ]

#Creating summary dataset of results for obtaining 95% CI
results <- data.frame(bootstrap  = rep(NA, 30),
                      Ya0g0  = rep(NA, 30),  Ya0g1  = rep(NA, 30),  Ya1g0  = rep(NA, 30),
                      Ya1g1  = rep(NA, 30),  pNDE   = rep(NA, 30),  tNDE   = rep(NA, 30),
                      pNIE   = rep(NA, 30),  tNIE   = rep(NA, 30),  TE     = rep(NA, 30),
                      pNDE_m = rep(NA, 30),  tNDE_m = rep(NA, 30),  pNIE_m = rep(NA, 30),
                      tNIE_m = rep(NA, 30),  TE_m   = rep(NA, 30))

#Creating progress bar for below (very long)
pb <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)

for (boot in 1:30) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Bootstrap number
  results$bootstrap[boot] <- boot
  
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset_4), NROW(dataset_4), replace = T)
  data_boot <- dataset_4[index_boot, ]
  
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
  #2a - Set baseline C to observed for each i at each time
  #2ai - Generate time t=0-(T-1) M & tv-C based on Step 1 model coefficients 
  #2aii - Assign time t A under intervention A=1
  #     - In this scenario if A is assumed to be fixed before all M, tv-C, Y, we can assign a or a* at baseline only
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
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1,
                           M_3 = pM_3a1)
  
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
  #   - Below is the predicted mediator for each i under A=0 at T=0 - T=4 (M_0 to M_4)
  #   - If predicted M for T=0-4 are taken from the same prediction, below 30 rows are selected for Step 3 K repeats
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0,
                           M_3 = pM_3a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    #3a/d - The below simulates Y at T=5 K times, saving the mean Y for each K in the above 'outcomes' dataset where A=0/1 and M=Ga based on M|A=1
    #     - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
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
    #Updating progress bar
    setTxtProgressBar(pk,k)
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
  #Adding above summary measures to 'results' dataset for 95% CI estimation
  results$Ya0g0[boot] <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot] <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot] <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot] <- mean(outcomes$Ya1g1)
  results$TE[boot]    <- mean(outcomes$TE)
  results$pNDE[boot]  <- mean(outcomes$pNDE)
  results$tNDE[boot]  <- mean(outcomes$tNDE)
  results$pNIE[boot]  <- mean(outcomes$pNIE)
  results$tNIE[boot]  <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Saving bootstrap results repeatedly to prevent crashing
  write.csv(results, paste0(directory, "boot_results_05052021_4_currY.csv"), row.names = FALSE, na = "")
} 
results_4 <- results

rm(pb, boot, index_boot, data_boot, results, 
   mA, mM_0, mM_1, mM_2, mM_3, mTV1_0, mTV1_1, mTV1_2, mTV1_3, 
   mTV2_0, mTV2_1, mTV2_2, mTV2_3,
   mY_0, mY_1, mY_2, mY_3, outcomes, pk, k, index, MC, 
   pM_0a1, pM_1a1, pM_2a1, pM_3a1, pTV1_0a1, pTV1_1a1, pTV1_2a1, pTV1_3a1,
   pTV2_0a1, pTV2_1a1, pTV2_2a1, pTV2_3a1, pY_0a1, pY_1a1, pY_2a1, pY_3a1, perm_Mt_a1,
   pM_0a0, pM_1a0, pM_2a0, pM_3a0, pTV1_0a0, pTV1_1a0, pTV1_2a0, pTV1_3a0,
   pTV2_0a0, pTV2_1a0, pTV2_2a0, pTV2_3a0, pY_0a0, pY_1a0, pY_2a0, pY_3a0, perm_Mt_a0,
   a, G, pM_0a1g, pM_1a1g, pM_2a1g, pM_3a1g, pTV1_0a1g, pTV1_1a1g, pTV1_2a1g, pTV1_3a1g,
   pTV2_0a1g, pTV2_1a1g, pTV2_2a1g, pTV2_3a1g, pY_0a1g, pY_1a1g, pY_2a1g, pY_3a1g,
   pM_0a0g, pM_1a0g, pM_2a0g, pM_3a0g, pTV1_0a0g, pTV1_1a0g, pTV1_2a0g, pTV1_3a0g,
   pTV2_0a0g, pTV2_1a0g, pTV2_2a0g, pTV2_3a0g, pY_0a0g, pY_1a0g, pY_2a0g, pY_3a0g)

#### 5 time-point marginal model #### 
#Restricting to those with 3 time-points
dataset_5 <- dataset[!is.na(dataset$Y_4), ]

#Creating summary dataset of results for obtaining 95% CI
results <- data.frame(bootstrap  = rep(NA, 30),
                      Ya0g0  = rep(NA, 30),  Ya0g1  = rep(NA, 30),  Ya1g0  = rep(NA, 30),
                      Ya1g1  = rep(NA, 30),  pNDE   = rep(NA, 30),  tNDE   = rep(NA, 30),
                      pNIE   = rep(NA, 30),  tNIE   = rep(NA, 30),  TE     = rep(NA, 30),
                      pNDE_m = rep(NA, 30),  tNDE_m = rep(NA, 30),  pNIE_m = rep(NA, 30),
                      tNIE_m = rep(NA, 30),  TE_m   = rep(NA, 30))

#Creating progress bar for below (very long)
pb <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)

for (boot in 1:30) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Bootstrap number
  results$bootstrap[boot] <- boot
  
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset_5), NROW(dataset_5), replace = T)
  data_boot <- dataset_5[index_boot, ]
  
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
  #2a - Set baseline C to observed for each i at each time
  #2ai - Generate time t=0-(T-1) M & tv-C based on Step 1 model coefficients 
  #2aii - Assign time t A under intervention A=1
  #     - In this scenario if A is assumed to be fixed before all M, tv-C, Y, we can assign a or a* at baseline only
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
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1,
                           M_3 = pM_3a1,
                           M_4 = pM_4a1)
  
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
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0,
                           M_3 = pM_3a0,
                           M_4 = pM_4a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    #3a/d - The below simulates Y at T=5 K times, saving the mean Y for each K in the above 'outcomes' dataset where A=0/1 and M=Ga based on M|A=1
    #     - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
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
    #Updating progress bar
    setTxtProgressBar(pk,k)
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
  #Adding above summary measures to 'results' dataset for 95% CI estimation
  results$Ya0g0[boot] <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot] <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot] <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot] <- mean(outcomes$Ya1g1)
  results$TE[boot]    <- mean(outcomes$TE)
  results$pNDE[boot]  <- mean(outcomes$pNDE)
  results$tNDE[boot]  <- mean(outcomes$tNDE)
  results$pNIE[boot]  <- mean(outcomes$pNIE)
  results$tNIE[boot]  <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Saving bootstrap results repeatedly to prevent crashing
  write.csv(results, paste0(directory, "boot_results_05052021_5_currY.csv"), row.names = FALSE, na = "")
} 
results_5 <- results

rm(pb, boot, index_boot, data_boot, results, 
   mA, mM_0, mM_1, mM_2, mM_3, mM_4, mTV1_0, mTV1_1, mTV1_2, mTV1_3, mTV1_4, 
   mTV2_0, mTV2_1, mTV2_2, mTV2_3, mTV2_4,
   mY_0, mY_1, mY_2, mY_3, mY_4, outcomes, pk, k, index, MC, 
   pM_0a1, pM_1a1, pM_2a1, pM_3a1, pM_4a1, pTV1_0a1, pTV1_1a1, pTV1_2a1, pTV1_3a1, pTV1_4a1,
   pTV2_0a1, pTV2_1a1, pTV2_2a1, pTV2_3a1, pTV2_4a1, pY_0a1, pY_1a1, pY_2a1, pY_3a1, pY_4a1, perm_Mt_a1,
   pM_0a0, pM_1a0, pM_2a0, pM_3a0, pM_4a0, pTV1_0a0, pTV1_1a0, pTV1_2a0, pTV1_3a0, pTV1_4a0,
   pTV2_0a0, pTV2_1a0, pTV2_2a0, pTV2_3a0, pTV2_4a0, pY_0a0, pY_1a0, pY_2a0, pY_3a0, pY_4a0, perm_Mt_a0,
   a, G, pM_0a1g, pM_1a1g, pM_2a1g, pM_3a1g, pM_4a1g, pTV1_0a1g, pTV1_1a1g, pTV1_2a1g, pTV1_3a1g, pTV1_4a1g,
   pTV2_0a1g, pTV2_1a1g, pTV2_2a1g, pTV2_3a1g, pTV2_4a1g, pY_0a1g, pY_1a1g, pY_2a1g, pY_3a1g, pY_4a1g,
   pM_0a0g, pM_1a0g, pM_2a0g, pM_3a0g, pM_4a0g, pTV1_0a0g, pTV1_1a0g, pTV1_2a0g, pTV1_3a0g, pTV1_4a0g,
   pTV2_0a0g, pTV2_1a0g, pTV2_2a0g, pTV2_3a0g, pTV2_4a0g, pY_0a0g, pY_1a0g, pY_2a0g, pY_3a0g, pY_4a0g)

#### 6 time-point marginal model #### 
#Restricting to those with 3 time-points
dataset_6 <- dataset[!is.na(dataset$Y_5), ]

#Creating summary dataset of results for obtaining 95% CI
results <- data.frame(bootstrap  = rep(NA, 30),
                      Ya0g0  = rep(NA, 30),  Ya0g1  = rep(NA, 30),  Ya1g0  = rep(NA, 30),
                      Ya1g1  = rep(NA, 30),  pNDE   = rep(NA, 30),  tNDE   = rep(NA, 30),
                      pNIE   = rep(NA, 30),  tNIE   = rep(NA, 30),  TE     = rep(NA, 30),
                      pNDE_m = rep(NA, 30),  tNDE_m = rep(NA, 30),  pNIE_m = rep(NA, 30),
                      tNIE_m = rep(NA, 30),  TE_m   = rep(NA, 30))

#Creating progress bar for below (very long)
pb <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)

for (boot in 1:30) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Bootstrap number
  results$bootstrap[boot] <- boot
  
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset_6), NROW(dataset_6), replace = T)
  data_boot <- dataset_6[index_boot, ]
  
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
  #2a - Set baseline C to observed for each i at each time
  #2ai - Generate time t=0-(T-1) M & tv-C based on Step 1 model coefficients 
  #2aii - Assign time t A under intervention A=1
  #     - In this scenario if A is assumed to be fixed before all M, tv-C, Y, we can assign a or a* at baseline only
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
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1,
                           M_3 = pM_3a1,
                           M_4 = pM_4a1,
                           M_5 = pM_5a1)
  
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
  #   - Below is the predicted mediator for each i under A=0 at T=0 - T=4 (M_0 to M_4)
  #   - If predicted M for T=0-4 are taken from the same prediction, below 30 rows are selected for Step 3 K repeats
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0,
                           M_3 = pM_3a0,
                           M_4 = pM_4a0,
                           M_5 = pM_5a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    #3a/d - The below simulates Y at T=5 K times, saving the mean Y for each K in the above 'outcomes' dataset where A=0/1 and M=Ga based on M|A=1
    #     - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
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
    #Updating progress bar
    setTxtProgressBar(pk,k)
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
  #Adding above summary measures to 'results' dataset for 95% CI estimation
  results$Ya0g0[boot] <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot] <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot] <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot] <- mean(outcomes$Ya1g1)
  results$TE[boot]    <- mean(outcomes$TE)
  results$pNDE[boot]  <- mean(outcomes$pNDE)
  results$tNDE[boot]  <- mean(outcomes$tNDE)
  results$pNIE[boot]  <- mean(outcomes$pNIE)
  results$tNIE[boot]  <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Saving bootstrap results repeatedly to prevent crashing
  write.csv(results, paste0(directory, "boot_results_05052021_6_currY.csv"), row.names = FALSE, na = "")
} 
results_6 <- results

rm(pb, boot, index_boot, data_boot, results, 
   mA, mM_0, mM_1, mM_2, mM_3, mM_4, mM_5, mTV1_0, mTV1_1, mTV1_2, mTV1_3, mTV1_4, mTV1_5,
   mTV2_0, mTV2_1, mTV2_2, mTV2_3, mTV2_4, mTV2_5,mY_0, mY_1, mY_2, mY_3, mY_4, mY_5, 
   outcomes, pk, k, index, MC, 
   pM_0a1, pM_1a1, pM_2a1, pM_3a1, pM_4a1, pM_5a1, pTV1_0a1, pTV1_1a1, pTV1_2a1, pTV1_3a1, pTV1_4a1, pTV1_5a1, 
   pTV2_0a1, pTV2_1a1, pTV2_2a1, pTV2_3a1, pTV2_4a1, pTV2_5a1, pY_0a1, pY_1a1, pY_2a1, pY_3a1, pY_4a1, pY_5a1, perm_Mt_a1,
   pM_0a0, pM_1a0, pM_2a0, pM_3a0, pM_4a0, pM_5a0, pTV1_0a0, pTV1_1a0, pTV1_2a0, pTV1_3a0, pTV1_4a0, pTV1_5a0,
   pTV2_0a0, pTV2_1a0, pTV2_2a0, pTV2_3a0, pTV2_4a0, pTV2_5a0, pY_0a0, pY_1a0, pY_2a0, pY_3a0, pY_4a0, pY_5a0, 
   perm_Mt_a0, a, G, pM_0a1g, pM_1a1g, pM_2a1g, pM_3a1g, pM_4a1g, pM_5a1g, pTV1_0a1g, pTV1_1a1g, pTV1_2a1g, pTV1_3a1g, pTV1_4a1g, pTV1_5a1g,
   pTV2_0a1g, pTV2_1a1g, pTV2_2a1g, pTV2_3a1g, pTV2_4a1g, pTV2_5a1g, pY_0a1g, pY_1a1g, pY_2a1g, pY_3a1g, pY_4a1g, pY_5a1g,
   pM_0a0g, pM_1a0g, pM_2a0g, pM_3a0g, pM_4a0g, pM_5a0g, pTV1_0a0g, pTV1_1a0g, pTV1_2a0g, pTV1_3a0g, pTV1_4a0g, pTV1_5a0g,
   pTV2_0a0g, pTV2_1a0g, pTV2_2a0g, pTV2_3a0g, pTV2_4a0g, pTV2_5a0g, pY_0a0g, pY_1a0g, pY_2a0g, pY_3a0g, pY_4a0g, pY_5a0g)

#### 7 time-point marginal model #### 
#Restricting to those with 3 time-points
dataset_7 <- dataset[!is.na(dataset$Y_6), ]

#Creating summary dataset of results for obtaining 95% CI
results <- data.frame(bootstrap  = rep(NA, 30),
                      Ya0g0  = rep(NA, 30),  Ya0g1  = rep(NA, 30),  Ya1g0  = rep(NA, 30),
                      Ya1g1  = rep(NA, 30),  pNDE   = rep(NA, 30),  tNDE   = rep(NA, 30),
                      pNIE   = rep(NA, 30),  tNIE   = rep(NA, 30),  TE     = rep(NA, 30),
                      pNDE_m = rep(NA, 30),  tNDE_m = rep(NA, 30),  pNIE_m = rep(NA, 30),
                      tNIE_m = rep(NA, 30),  TE_m   = rep(NA, 30))

#Creating progress bar for below (very long)
pb <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)

for (boot in 1:30) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Bootstrap number
  results$bootstrap[boot] <- boot
  
  #Creating bootstrap resample of original dataset
  index_boot <- sample(1:NROW(dataset_7), NROW(dataset_7), replace = T)
  data_boot <- dataset_7[index_boot, ]
  
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
  #2a - Set baseline C to observed for each i at each time
  #2ai - Generate time t=0-(T-1) M & tv-C based on Step 1 model coefficients 
  #2aii - Assign time t A under intervention A=1
  #     - In this scenario if A is assumed to be fixed before all M, tv-C, Y, we can assign a or a* at baseline only
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
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1,
                           M_3 = pM_3a1,
                           M_4 = pM_4a1,
                           M_5 = pM_5a1,
                           M_6 = pM_6a1)
  
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
  #   - Below is the predicted mediator for each i under A=0 at T=0 - T=4 (M_0 to M_4)
  #   - If predicted M for T=0-4 are taken from the same prediction, below 30 rows are selected for Step 3 K repeats
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0,
                           M_3 = pM_3a0,
                           M_4 = pM_4a0,
                           M_5 = pM_5a0,
                           M_6 = pM_6a0)
  
  #3d - Repeat 1-3 K (30) times using different permutations of M
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  for (k in 1:30) {
    #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) K times using different 2b/d M for Ga/ Ga*
    #3a/d - The below simulates Y at T=5 K times, saving the mean Y for each K in the above 'outcomes' dataset where A=0/1 and M=Ga based on M|A=1
    #     - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
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
    #Updating progress bar
    setTxtProgressBar(pk,k)
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
  #Adding above summary measures to 'results' dataset for 95% CI estimation
  results$Ya0g0[boot] <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot] <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot] <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot] <- mean(outcomes$Ya1g1)
  results$TE[boot]    <- mean(outcomes$TE)
  results$pNDE[boot]  <- mean(outcomes$pNDE)
  results$tNDE[boot]  <- mean(outcomes$tNDE)
  results$pNIE[boot]  <- mean(outcomes$pNIE)
  results$tNIE[boot]  <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Saving bootstrap results repeatedly to prevent crashing
  write.csv(results, paste0(directory, "boot_results_05052021_7_currY.csv"), row.names = FALSE, na = "")
} 
results_7 <- results

rm(pb, boot, index_boot, data_boot, results, 
   mA, mM_0, mM_1, mM_2, mM_3, mM_4, mM_5, mM_6, mTV1_0, mTV1_1, mTV1_2, mTV1_3, mTV1_4, mTV1_5, mTV1_6,
   mTV2_0, mTV2_1, mTV2_2, mTV2_3, mTV2_4, mTV2_5, mTV2_6, mY_0, mY_1, mY_2, mY_3, mY_4, mY_5, mY_6,
   outcomes, pk, k, index, MC, pM_0a1, pM_1a1, pM_2a1, pM_3a1, pM_4a1, pM_5a1, pM_6a1, 
   pTV1_0a1, pTV1_1a1, pTV1_2a1, pTV1_3a1, pTV1_4a1, pTV1_5a1, pTV1_6a1,
   pTV2_0a1, pTV2_1a1, pTV2_2a1, pTV2_3a1, pTV2_4a1, pTV2_5a1, pTV2_6a1, pY_0a1, pY_1a1, pY_2a1, pY_3a1, pY_4a1, pY_5a1, pY_6a1, perm_Mt_a1,
   pM_0a0, pM_1a0, pM_2a0, pM_3a0, pM_4a0, pM_5a0, pM_6a0, pTV1_0a0, pTV1_1a0, pTV1_2a0, pTV1_3a0, pTV1_4a0, pTV1_5a0, pTV1_6a0,
   pTV2_0a0, pTV2_1a0, pTV2_2a0, pTV2_3a0, pTV2_4a0, pTV2_5a0, pTV2_6a0, pY_0a0, pY_1a0, pY_2a0, pY_3a0, pY_4a0, pY_5a0, pY_6a0, 
   perm_Mt_a0, a, G, pM_0a1g, pM_1a1g, pM_2a1g, pM_3a1g, pM_4a1g, pM_5a1g, pM_6a1g, 
   pTV1_0a1g, pTV1_1a1g, pTV1_2a1g, pTV1_3a1g, pTV1_4a1g, pTV1_5a1g, pTV1_6a1g,
   pTV2_0a1g, pTV2_1a1g, pTV2_2a1g, pTV2_3a1g, pTV2_4a1g, pTV2_5a1g, pTV2_6a1g, pY_0a1g, pY_1a1g, pY_2a1g, pY_3a1g, pY_4a1g, pY_5a1g, pY_6a1g,
   pM_0a0g, pM_1a0g, pM_2a0g, pM_3a0g, pM_4a0g, pM_5a0g, pM_6a0g, pTV1_0a0g, pTV1_1a0g, pTV1_2a0g, pTV1_3a0g, pTV1_4a0g, pTV1_5a0g, pTV1_6a0g,
   pTV2_0a0g, pTV2_1a0g, pTV2_2a0g, pTV2_3a0g, pTV2_4a0g, pTV2_5a0g, pTV2_6a0g, pY_0a0g, pY_1a0g, pY_2a0g, pY_3a0g, pY_4a0g, pY_5a0g, pY_6a0g)
