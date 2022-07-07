#### Mediational g-formula functions script - marginal models - same timepoint Y #####
#Author - Kieran Blaikie
#Date - 11 Oct 2021

#Overview - This script creates 3-6 time-point marginal models with 95% CI for both
#           primary analyses and complete-case analyses
#            - Outcome Y is treated as binary
#            - Mediator M is treatd as continuous
#            - Exposure A is treated as binary
#            - 9 Baseline confounders C are included, named C1 - C9
#            - 3 Time-Varying confounders TV are included, named TV1-TV3 in primary analyses
#            - 2 Time-Varying confounders TV are included, named TV1-TV2 in the complete-case analyses
#            - The causal order is: C > HS > rep[M > TV1 > TV2 > TV3 > Y] > ... > T

#Changes - Compared to the functions created in the July 15th version of 'making_medgformula_functions',
#          the following changes are made:
#           1 - Incorporated estimates from the observed data (i.e. fitting each model on the observed
#               data instead of resampled data), given as the first row in the 'results' dataframe
#           2 - Permuted the mediator at each time-point instead of per-record, separating the 
#               correlated structure of each mediator over time
#           3 - Stored the mediator error term separately for each T model instead of using
#               predict(, se.fit = T), as this is computationally faster
#           4 - Removed the Step 3d K = 1:30 repetition, as through using a large MC sample this 
#               step isn't necessary
#           5 - For each final time-point outcome, the prediction is left unresolved, as this
#               should reduce the needed size of the MC sample
#           6 - Natural Course simulations are added assigning exposure as observed instead of
#               simulated. As exposure is static, simulating exposure isn't necessary.
#           7 - Based on the 'MC_sample_check' script, the MC size is increased from 10,000 to 50,000.
#               This shouldn't affect point-estimates, but should reduce simulation error in our estimates
#        - Compared to the functions created in the Sep 7th version of 'making_medgformula_functions',
#           1 - Based on the 'checking_model_assumptions_093021.R' script, a 2nd order polynomial
#               is included for M in models for TV1_0 (i.e. Time 0 Unemployed/NILF status) and 
#               for M_t-1 in models for M_t

#### Creating mediational g-formula functions ####
#### Primary Analysis - 3 time-point ####
medgf_3 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 3 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                                             poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV3_0 <- glm(TV3_0 ~                                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                                     TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                               Y_0 + TV3_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                         M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~                 TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV3_1 <- glm(TV3_1 ~         TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV3_1 + TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                               Y_1 + TV3_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                         M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~                 TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV3_2 <- glm(TV3_2 ~         TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_2   <- glm(Y_2   ~ TV3_2 + TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                                                                 MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                   A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pTV3_0nat   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0nat   <- noisy_est("binary", pTV3_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                                            Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pTV3_1nat   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1nat   <- noisy_est("binary", pTV3_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                                            Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pTV3_2nat   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2nat   <- noisy_est("binary", pTV3_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A =  pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                            MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                          M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                                                                      TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pTV3_0natA   <- predict(mTV3_0, newdata = data.frame(                                                                                                                  TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0natA   <- noisy_est("binary", pTV3_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                                                TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                                                Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                                              M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                                          TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pTV3_1natA   <- predict(mTV3_1, newdata = data.frame(                      TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1natA   <- noisy_est("binary", pTV3_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                                                Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                                              M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                                          TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pTV3_2natA   <- predict(mTV3_2, newdata = data.frame(                      TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2natA   <- noisy_est("binary", pTV3_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_2nat overall, by exposure
  MC$Y_nat <- pY_2nat
  MC$YA_nat <- pY_2natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pTV3_0a1   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a1   <- noisy_est("binary", pTV3_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pTV3_1a1   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a1   <- noisy_est("binary", pTV3_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pTV3_0a0   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a0   <- noisy_est("binary", pTV3_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pTV3_1a0   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a0   <- noisy_est("binary", pTV3_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
    }
    
    #Q(a,a) or Q(a,a*)
    pM_0a1g     <- G$M_0
    pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
    pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
    pTV3_0a1g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a1g   <- noisy_est("binary", pTV3_0a1g)
    pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a1g     <- noisy_est("binary", pY_0a1g)
    
    pM_1a1g     <- G$M_1
    pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
    pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
    pTV3_1a1g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a1g   <- noisy_est("binary", pTV3_1a1g)
    pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a1g     <- noisy_est("binary", pY_1a1g)
    
    pM_2a1g     <- G$M_2
    pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
    pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
    pTV3_2a1g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a1g   <- noisy_est("binary", pTV3_2a1g)
    
    #Q(a*,a) or Q(a*,a*)
    pM_0a0g     <- G$M_0
    pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
    pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
    pTV3_0a0g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a0g   <- noisy_est("binary", pTV3_0a0g)
    pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a0g     <- noisy_est("binary", pY_0a0g)
    
    pM_1a0g     <- G$M_1
    pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
    pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
    pTV3_1a0g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a0g   <- noisy_est("binary", pTV3_1a0g)
    pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a0g     <- noisy_est("binary", pY_1a0g)
    
    pM_2a0g     <- G$M_2
    pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
    pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
    pTV3_2a0g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a0g   <- noisy_est("binary", pTV3_2a0g)
    
    #3b - Simulate outcome for each i and save mean outcome
    pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_2a1g <- pY_2a1g
    MC$pY_2a0g <- pY_2a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_2a1g)
      outcomes$Ya0g1 <- mean(MC$pY_2a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_2a1g)
      outcomes$Ya0g0 <- mean(MC$pY_2a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}



#### Primary Analysis - 4 time-point ####
medgf_4 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 4 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                                             poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV3_0 <- glm(TV3_0 ~                                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                                     TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                               Y_0 + TV3_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                         M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~                 TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV3_1 <- glm(TV3_1 ~         TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV3_1 + TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                               Y_1 + TV3_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                         M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~                 TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV3_2 <- glm(TV3_2 ~         TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV3_2 + TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                               Y_2 + TV3_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                         M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~                 TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV3_3 <- glm(TV3_3 ~         TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_3   <- glm(Y_3   ~ TV3_3 + TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                                                                 MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                   A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pTV3_0nat   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0nat   <- noisy_est("binary", pTV3_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                                            Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pTV3_1nat   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1nat   <- noisy_est("binary", pTV3_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                                            Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pTV3_2nat   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2nat   <- noisy_est("binary", pTV3_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                                            Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pTV3_3nat   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3nat   <- noisy_est("binary", pTV3_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                            MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                          M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                                                                      TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pTV3_0natA   <- predict(mTV3_0, newdata = data.frame(                                                                                                                  TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0natA   <- noisy_est("binary", pTV3_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                                                TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                                                Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                                              M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                                          TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pTV3_1natA   <- predict(mTV3_1, newdata = data.frame(                      TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1natA   <- noisy_est("binary", pTV3_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                                                Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                                              M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                                          TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pTV3_2natA   <- predict(mTV3_2, newdata = data.frame(                      TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2natA   <- noisy_est("binary", pTV3_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                                                Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                                              M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                                          TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pTV3_3natA   <- predict(mTV3_3, newdata = data.frame(                      TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3natA   <- noisy_est("binary", pTV3_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_3nat overall, by exposure
  MC$Y_nat <- pY_3nat
  MC$YA_nat <- pY_3natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pTV3_0a1   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a1   <- noisy_est("binary", pTV3_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pTV3_1a1   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a1   <- noisy_est("binary", pTV3_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pTV3_2a1   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a1   <- noisy_est("binary", pTV3_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pTV3_0a0   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a0   <- noisy_est("binary", pTV3_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pTV3_1a0   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a0   <- noisy_est("binary", pTV3_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pTV3_2a0   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a0   <- noisy_est("binary", pTV3_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
    }
    
    #Q(a,a) or Q(a,a*)
    pM_0a1g     <- G$M_0
    pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
    pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
    pTV3_0a1g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a1g   <- noisy_est("binary", pTV3_0a1g)
    pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a1g     <- noisy_est("binary", pY_0a1g)
    
    pM_1a1g     <- G$M_1
    pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
    pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
    pTV3_1a1g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a1g   <- noisy_est("binary", pTV3_1a1g)
    pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a1g     <- noisy_est("binary", pY_1a1g)
    
    pM_2a1g     <- G$M_2
    pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
    pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
    pTV3_2a1g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a1g   <- noisy_est("binary", pTV3_2a1g)
    pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a1g     <- noisy_est("binary", pY_2a1g)
    
    pM_3a1g     <- G$M_3
    pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
    pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
    pTV3_3a1g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a1g   <- noisy_est("binary", pTV3_3a1g)
    
    #Q(a*,a) or Q(a*,a*)
    pM_0a0g     <- G$M_0
    pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
    pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
    pTV3_0a0g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a0g   <- noisy_est("binary", pTV3_0a0g)
    pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a0g     <- noisy_est("binary", pY_0a0g)
    
    pM_1a0g     <- G$M_1
    pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
    pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
    pTV3_1a0g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a0g   <- noisy_est("binary", pTV3_1a0g)
    pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a0g     <- noisy_est("binary", pY_1a0g)
    
    pM_2a0g     <- G$M_2
    pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
    pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
    pTV3_2a0g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a0g   <- noisy_est("binary", pTV3_2a0g)
    pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a0g     <- noisy_est("binary", pY_2a0g)
    
    pM_3a0g     <- G$M_3
    pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
    pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
    pTV3_3a0g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a0g   <- noisy_est("binary", pTV3_3a0g)
    
    #3b - Simulate outcome for each i and save mean outcome
    pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_3a1g <- pY_3a1g
    MC$pY_3a0g <- pY_3a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_3a1g)
      outcomes$Ya0g1 <- mean(MC$pY_3a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_3a1g)
      outcomes$Ya0g0 <- mean(MC$pY_3a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}


#### Primary Analysis - 5 time-point ####
medgf_5 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 5 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                                             poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV3_0 <- glm(TV3_0 ~                                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                                     TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                               Y_0 + TV3_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                         M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~                 TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV3_1 <- glm(TV3_1 ~         TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV3_1 + TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                               Y_1 + TV3_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                         M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~                 TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV3_2 <- glm(TV3_2 ~         TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV3_2 + TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                               Y_2 + TV3_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                         M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~                 TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV3_3 <- glm(TV3_3 ~         TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV3_3 + TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                               Y_3 + TV3_3 + TV2_3 + TV1_3 + poly(M_3, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3,                           data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                         M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~                 TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV3_4 <- glm(TV3_4 ~         TV2_4 + TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_4   <- glm(Y_4   ~ TV3_4 + TV2_4 + TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  eM_4 <- sqrt(summary(mM_4)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                                                                 MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                   A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pTV3_0nat   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0nat   <- noisy_est("binary", pTV3_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                                            Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pTV3_1nat   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1nat   <- noisy_est("binary", pTV3_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                                            Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pTV3_2nat   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2nat   <- noisy_est("binary", pTV3_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                                            Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pTV3_3nat   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3nat   <- noisy_est("binary", pTV3_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3nat     <- noisy_est("binary", pY_3nat)
  
  pM_4nat     <- predict(mM_4, newdata = data.frame(                                                                            Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4nat     <- noisy_est("continuous", pM_4nat, eM_4)
  pTV1_4nat   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4nat   <- noisy_est("binary", pTV1_4nat)
  pTV2_4nat   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4nat   <- noisy_est("binary", pTV2_4nat)
  pTV3_4nat   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4nat   <- noisy_est("binary", pTV3_4nat)
  pY_4nat     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                            MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                          M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                                                                      TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pTV3_0natA   <- predict(mTV3_0, newdata = data.frame(                                                                                                                  TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0natA   <- noisy_est("binary", pTV3_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                                                TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                                                Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                                              M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                                          TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pTV3_1natA   <- predict(mTV3_1, newdata = data.frame(                      TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1natA   <- noisy_est("binary", pTV3_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                                                Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                                              M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                                          TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pTV3_2natA   <- predict(mTV3_2, newdata = data.frame(                      TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2natA   <- noisy_est("binary", pTV3_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                                                Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                                              M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                                          TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pTV3_3natA   <- predict(mTV3_3, newdata = data.frame(                      TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3natA   <- noisy_est("binary", pTV3_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3natA     <- noisy_est("binary", pY_3natA)
  
  pM_4natA     <- predict(mM_4, newdata = data.frame(                                                                                Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4natA     <- noisy_est("continuous", pM_4natA, eM_4)
  pTV1_4natA   <- predict(mTV1_4, newdata = data.frame(                                                              M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4natA   <- noisy_est("binary", pTV1_4natA)
  pTV2_4natA   <- predict(mTV2_4, newdata = data.frame(                                          TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4natA   <- noisy_est("binary", pTV2_4natA)
  pTV3_4natA   <- predict(mTV3_4, newdata = data.frame(                      TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4natA   <- noisy_est("binary", pTV3_4natA)
  pY_4natA     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_4nat overall, by exposure
  MC$Y_nat <- pY_4nat
  MC$YA_nat <- pY_4natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pTV3_0a1   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a1   <- noisy_est("binary", pTV3_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pTV3_1a1   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a1   <- noisy_est("binary", pTV3_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pTV3_2a1   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a1   <- noisy_est("binary", pTV3_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                                        M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                                      TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pTV3_3a1   <- predict(mTV3_3, newdata = data.frame(                    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3a1   <- noisy_est("binary", pTV3_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                                        Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a1     <- noisy_est("continuous", pM_4a1, eM_4)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pTV3_0a0   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a0   <- noisy_est("binary", pTV3_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pTV3_1a0   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a0   <- noisy_est("binary", pTV3_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pTV3_2a0   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a0   <- noisy_est("binary", pTV3_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                                        M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                                      TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pTV3_3a0   <- predict(mTV3_3, newdata = data.frame(                    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3a0   <- noisy_est("binary", pTV3_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                                        Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a0     <- noisy_est("continuous", pM_4a0, eM_4)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
    }
    
    #Q(a,a) or Q(a,a*)
    pM_0a1g     <- G$M_0
    pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
    pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
    pTV3_0a1g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a1g   <- noisy_est("binary", pTV3_0a1g)
    pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a1g     <- noisy_est("binary", pY_0a1g)
    
    pM_1a1g     <- G$M_1
    pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
    pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
    pTV3_1a1g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a1g   <- noisy_est("binary", pTV3_1a1g)
    pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a1g     <- noisy_est("binary", pY_1a1g)
    
    pM_2a1g     <- G$M_2
    pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
    pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
    pTV3_2a1g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a1g   <- noisy_est("binary", pTV3_2a1g)
    pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a1g     <- noisy_est("binary", pY_2a1g)
    
    pM_3a1g     <- G$M_3
    pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
    pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
    pTV3_3a1g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a1g   <- noisy_est("binary", pTV3_3a1g)
    pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_3a1g     <- noisy_est("binary", pY_3a1g)
    
    pM_4a1g     <- G$M_4
    pTV1_4a1g   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_4a1g   <- noisy_est("binary", pTV1_4a1g)
    pTV2_4a1g   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_4a1g   <- noisy_est("binary", pTV2_4a1g)
    pTV3_4a1g   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_4a1g   <- noisy_est("binary", pTV3_4a1g)
    
    #Q(a*,a) or Q(a*,a*)
    pM_0a0g     <- G$M_0
    pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
    pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
    pTV3_0a0g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a0g   <- noisy_est("binary", pTV3_0a0g)
    pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a0g     <- noisy_est("binary", pY_0a0g)
    
    pM_1a0g     <- G$M_1
    pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
    pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
    pTV3_1a0g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a0g   <- noisy_est("binary", pTV3_1a0g)
    pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a0g     <- noisy_est("binary", pY_1a0g)
    
    pM_2a0g     <- G$M_2
    pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
    pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
    pTV3_2a0g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a0g   <- noisy_est("binary", pTV3_2a0g)
    pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a0g     <- noisy_est("binary", pY_2a0g)
    
    pM_3a0g     <- G$M_3
    pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
    pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
    pTV3_3a0g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a0g   <- noisy_est("binary", pTV3_3a0g)
    pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_3a0g     <- noisy_est("binary", pY_3a0g)
    
    pM_4a0g     <- G$M_4
    pTV1_4a0g   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_4a0g   <- noisy_est("binary", pTV1_4a0g)
    pTV2_4a0g   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_4a0g   <- noisy_est("binary", pTV2_4a0g)
    pTV3_4a0g   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_4a0g   <- noisy_est("binary", pTV3_4a0g)
    
    #3b - Simulate outcome for each i and save mean outcome
    pY_4a1g     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_4a1g <- pY_4a1g
    MC$pY_4a0g <- pY_4a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_4a1g)
      outcomes$Ya0g1 <- mean(MC$pY_4a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_4a1g)
      outcomes$Ya0g0 <- mean(MC$pY_4a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}

#### Primary Analysis - 6 time-point ####
medgf_6 <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 6 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                                             poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV3_0 <- glm(TV3_0 ~                                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                                     TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                               Y_0 + TV3_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                         M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~                 TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV3_1 <- glm(TV3_1 ~         TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV3_1 + TV2_1 + TV1_1 + M_1 + Y_0 + TV3_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                               Y_1 + TV3_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                         M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~                 TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV3_2 <- glm(TV3_2 ~         TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV3_2 + TV2_2 + TV1_2 + M_2 + Y_1 + TV3_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                               Y_2 + TV3_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                         M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~                 TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV3_3 <- glm(TV3_3 ~         TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV3_3 + TV2_3 + TV1_3 + M_3 + Y_2 + TV3_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                               Y_3 + TV3_3 + TV2_3 + TV1_3 + poly(M_3, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3,                           data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                         M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~                 TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV3_4 <- glm(TV3_4 ~         TV2_4 + TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mY_4   <- glm(Y_4   ~ TV3_4 + TV2_4 + TV1_4 + M_4 + Y_3 + TV3_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mM_5   <- glm(M_5   ~                               Y_4 + TV3_4 + TV2_4 + TV1_4 + poly(M_4, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4,                           data = data_boot, family = gaussian("identity"))
  mTV1_5 <- glm(TV1_5 ~                         M_5 + Y_4 + TV3_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mTV2_5 <- glm(TV2_5 ~                 TV1_5 + M_5 + Y_4 + TV3_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mTV3_5 <- glm(TV3_5 ~         TV2_5 + TV1_5 + M_5 + Y_4 + TV3_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_5   <- glm(Y_5   ~ TV3_5 + TV2_5 + TV1_5 + M_5 + Y_4 + TV3_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  eM_4 <- sqrt(summary(mM_4)$dispersion)
  eM_5 <- sqrt(summary(mM_5)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                                                                 MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                   A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pTV3_0nat   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0nat   <- noisy_est("binary", pTV3_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                                            Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pTV3_1nat   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1nat   <- noisy_est("binary", pTV3_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV3_0 = pTV3_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                                            Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pTV3_2nat   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2nat   <- noisy_est("binary", pTV3_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV3_1 = pTV3_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                                            Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pTV3_3nat   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3nat   <- noisy_est("binary", pTV3_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV3_2 = pTV3_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3nat     <- noisy_est("binary", pY_3nat)
  
  pM_4nat     <- predict(mM_4, newdata = data.frame(                                                                            Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4nat     <- noisy_est("continuous", pM_4nat, eM_4)
  pTV1_4nat   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4nat   <- noisy_est("binary", pTV1_4nat)
  pTV2_4nat   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4nat   <- noisy_est("binary", pTV2_4nat)
  pTV3_4nat   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4nat   <- noisy_est("binary", pTV3_4nat)
  pY_4nat     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV3_3 = pTV3_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4nat     <- noisy_est("binary", pY_4nat)
  
  pM_5nat     <- predict(mM_5, newdata = data.frame(                                                                            Y_4 = pY_4nat, TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5nat     <- noisy_est("continuous", pM_5nat, eM_5)
  pTV1_5nat   <- predict(mTV1_5, newdata = data.frame(                                                           M_5 = pM_5nat, Y_4 = pY_4nat, TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5nat   <- noisy_est("binary", pTV1_5nat)
  pTV2_5nat   <- predict(mTV2_5, newdata = data.frame(                                        TV1_5 = pTV1_5nat, M_5 = pM_5nat, Y_4 = pY_4nat, TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5nat   <- noisy_est("binary", pTV2_5nat)
  pTV3_5nat   <- predict(mTV3_5, newdata = data.frame(                     TV2_5 = pTV2_5nat, TV1_5 = pTV1_5nat, M_5 = pM_5nat, Y_4 = pY_4nat, TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_5nat   <- noisy_est("binary", pTV3_5nat)
  pY_5nat     <- predict(mY_5, newdata = data.frame(    TV3_5 = pTV3_5nat, TV2_5 = pTV2_5nat, TV1_5 = pTV1_5nat, M_5 = pM_5nat, Y_4 = pY_4nat, TV3_4 = pTV3_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                                                            MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                          M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                                                                      TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pTV3_0natA   <- predict(mTV3_0, newdata = data.frame(                                                                                                                  TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0natA   <- noisy_est("binary", pTV3_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                                                TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                                                Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                                              M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                                          TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pTV3_1natA   <- predict(mTV3_1, newdata = data.frame(                      TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1natA   <- noisy_est("binary", pTV3_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV3_0 = pTV3_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                                                Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                                              M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                                          TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pTV3_2natA   <- predict(mTV3_2, newdata = data.frame(                      TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2natA   <- noisy_est("binary", pTV3_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV3_1 = pTV3_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                                                Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                                              M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                                          TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pTV3_3natA   <- predict(mTV3_3, newdata = data.frame(                      TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3natA   <- noisy_est("binary", pTV3_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV3_2 = pTV3_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3natA     <- noisy_est("binary", pY_3natA)
  
  pM_4natA     <- predict(mM_4, newdata = data.frame(                                                                                Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4natA     <- noisy_est("continuous", pM_4natA, eM_4)
  pTV1_4natA   <- predict(mTV1_4, newdata = data.frame(                                                              M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4natA   <- noisy_est("binary", pTV1_4natA)
  pTV2_4natA   <- predict(mTV2_4, newdata = data.frame(                                          TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4natA   <- noisy_est("binary", pTV2_4natA)
  pTV3_4natA   <- predict(mTV3_4, newdata = data.frame(                      TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4natA   <- noisy_est("binary", pTV3_4natA)
  pY_4natA     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV3_3 = pTV3_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4natA     <- noisy_est("binary", pY_4natA)
  
  pM_5natA     <- predict(mM_5, newdata = data.frame(                                                                                Y_4 = pY_4natA, TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5natA     <- noisy_est("continuous", pM_5natA, eM_5)
  pTV1_5natA   <- predict(mTV1_5, newdata = data.frame(                                                              M_5 = pM_5natA, Y_4 = pY_4natA, TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5natA   <- noisy_est("binary", pTV1_5natA)
  pTV2_5natA   <- predict(mTV2_5, newdata = data.frame(                                          TV1_5 = pTV1_5natA, M_5 = pM_5natA, Y_4 = pY_4natA, TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5natA   <- noisy_est("binary", pTV2_5natA)
  pTV3_5natA   <- predict(mTV3_5, newdata = data.frame(                      TV2_5 = pTV2_5natA, TV1_5 = pTV1_5natA, M_5 = pM_5natA, Y_4 = pY_4natA, TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_5natA   <- noisy_est("binary", pTV3_5natA)
  pY_5natA     <- predict(mY_5, newdata = data.frame(    TV3_5 = pTV3_5natA, TV2_5 = pTV2_5natA, TV1_5 = pTV1_5natA, M_5 = pM_5natA, Y_4 = pY_4natA, TV3_4 = pTV3_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_5nat overall, by exposure
  MC$Y_nat <- pY_5nat
  MC$YA_nat <- pY_5natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pTV3_0a1   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a1   <- noisy_est("binary", pTV3_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pTV3_1a1   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a1   <- noisy_est("binary", pTV3_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV3_0 = pTV3_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pTV3_2a1   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a1   <- noisy_est("binary", pTV3_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV3_1 = pTV3_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                                        M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                                      TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pTV3_3a1   <- predict(mTV3_3, newdata = data.frame(                    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3a1   <- noisy_est("binary", pTV3_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV3_2 = pTV3_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                                        Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a1     <- noisy_est("continuous", pM_4a1, eM_4)
  pTV1_4a1   <- predict(mTV1_4, newdata = data.frame(                                                        M_4 = pM_4a1, Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a1   <- noisy_est("binary", pTV1_4a1)
  pTV2_4a1   <- predict(mTV2_4, newdata = data.frame(                                      TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a1   <- noisy_est("binary", pTV2_4a1)
  pTV3_4a1   <- predict(mTV3_4, newdata = data.frame(                    TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4a1   <- noisy_est("binary", pTV3_4a1)
  pY_4a1     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV3_3 = pTV3_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a1     <- noisy_est("binary", pY_4a1)
  
  pM_5a1     <- predict(mM_5, newdata = data.frame(                                                                        Y_4 = pY_4a1, TV3_4 = pTV3_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5a1     <- noisy_est("continuous", pM_5a1, eM_5)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_5 = pM_5a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                                                          A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                          M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                                                        TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pTV3_0a0   <- predict(mTV3_0, newdata = data.frame(                                                                                                      TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_0a0   <- noisy_est("binary", pTV3_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                                      TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                                        Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                                        M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                                      TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pTV3_1a0   <- predict(mTV3_1, newdata = data.frame(                    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_1a0   <- noisy_est("binary", pTV3_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV3_0 = pTV3_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                                        Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                                        M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                                      TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pTV3_2a0   <- predict(mTV3_2, newdata = data.frame(                    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_2a0   <- noisy_est("binary", pTV3_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV3_1 = pTV3_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                                        Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                                        M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                                      TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pTV3_3a0   <- predict(mTV3_3, newdata = data.frame(                    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_3a0   <- noisy_est("binary", pTV3_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV3_2 = pTV3_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                                        Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a0     <- noisy_est("continuous", pM_4a0, eM_4)
  pTV1_4a0   <- predict(mTV1_4, newdata = data.frame(                                                        M_4 = pM_4a0, Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a0   <- noisy_est("binary", pTV1_4a0)
  pTV2_4a0   <- predict(mTV2_4, newdata = data.frame(                                      TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a0   <- noisy_est("binary", pTV2_4a0)
  pTV3_4a0   <- predict(mTV3_4, newdata = data.frame(                    TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV3_4a0   <- noisy_est("binary", pTV3_4a0)
  pY_4a0     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV3_3 = pTV3_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a0     <- noisy_est("binary", pY_4a0)
  
  pM_5a0     <- predict(mM_5, newdata = data.frame(                                                                        Y_4 = pY_4a0, TV3_4 = pTV3_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5a0     <- noisy_est("continuous", pM_5a0, eM_5)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_5 = pM_5a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
    }
    
    #Q(a,a) or Q(a,a*)
    pM_0a1g     <- G$M_0
    pTV1_0a1g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a1g   <- noisy_est("binary", pTV1_0a1g)
    pTV2_0a1g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a1g   <- noisy_est("binary", pTV2_0a1g)
    pTV3_0a1g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a1g   <- noisy_est("binary", pTV3_0a1g)
    pY_0a1g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a1g     <- noisy_est("binary", pY_0a1g)
    
    pM_1a1g     <- G$M_1
    pTV1_1a1g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a1g   <- noisy_est("binary", pTV1_1a1g)
    pTV2_1a1g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a1g   <- noisy_est("binary", pTV2_1a1g)
    pTV3_1a1g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a1g   <- noisy_est("binary", pTV3_1a1g)
    pY_1a1g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, Y_0 = pY_0a1g, TV3_0 = pTV3_0a1g, TV2_0 = pTV2_0a1g, TV1_0 = pTV1_0a1g, M_0 = pM_0a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a1g     <- noisy_est("binary", pY_1a1g)
    
    pM_2a1g     <- G$M_2
    pTV1_2a1g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a1g   <- noisy_est("binary", pTV1_2a1g)
    pTV2_2a1g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a1g   <- noisy_est("binary", pTV2_2a1g)
    pTV3_2a1g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a1g   <- noisy_est("binary", pTV3_2a1g)
    pY_2a1g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, Y_1 = pY_1a1g, TV3_1 = pTV3_1a1g, TV2_1 = pTV2_1a1g, TV1_1 = pTV1_1a1g, M_1 = pM_1a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a1g     <- noisy_est("binary", pY_2a1g)
    
    pM_3a1g     <- G$M_3
    pTV1_3a1g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a1g   <- noisy_est("binary", pTV1_3a1g)
    pTV2_3a1g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a1g   <- noisy_est("binary", pTV2_3a1g)
    pTV3_3a1g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a1g   <- noisy_est("binary", pTV3_3a1g)
    pY_3a1g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, Y_2 = pY_2a1g, TV3_2 = pTV3_2a1g, TV2_2 = pTV2_2a1g, TV1_2 = pTV1_2a1g, M_2 = pM_2a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_3a1g     <- noisy_est("binary", pY_3a1g)
    
    pM_4a1g     <- G$M_4
    pTV1_4a1g   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_4a1g   <- noisy_est("binary", pTV1_4a1g)
    pTV2_4a1g   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_4a1g   <- noisy_est("binary", pTV2_4a1g)
    pTV3_4a1g   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_4a1g   <- noisy_est("binary", pTV3_4a1g)
    pY_4a1g     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, Y_3 = pY_3a1g, TV3_3 = pTV3_3a1g, TV2_3 = pTV2_3a1g, TV1_3 = pTV1_3a1g, M_3 = pM_3a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_4a1g     <- noisy_est("binary", pY_4a1g)
    
    pM_5a1g     <- G$M_5
    pTV1_5a1g   <- predict(mTV1_5, newdata = data.frame(                                                           M_5 = pM_5a1g, Y_4 = pY_4a1g, TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_5a1g   <- noisy_est("binary", pTV1_5a1g)
    pTV2_5a1g   <- predict(mTV2_5, newdata = data.frame(                                        TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_5a1g   <- noisy_est("binary", pTV2_5a1g)
    pTV3_5a1g   <- predict(mTV3_5, newdata = data.frame(                     TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_5a1g   <- noisy_est("binary", pTV3_5a1g)
    
    #Q(a*,a) or Q(a*,a*)
    pM_0a0g     <- G$M_0
    pTV1_0a0g   <- predict(mTV1_0, newdata = data.frame(                                                                                                                                                  M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_0a0g   <- noisy_est("binary", pTV1_0a0g)
    pTV2_0a0g   <- predict(mTV2_0, newdata = data.frame(                                                                                                                               TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_0a0g   <- noisy_est("binary", pTV2_0a0g)
    pTV3_0a0g   <- predict(mTV3_0, newdata = data.frame(                                                                                                            TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_0a0g   <- noisy_est("binary", pTV3_0a0g)
    pY_0a0g     <- predict(mY_0, newdata = data.frame(                                                                                           TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_0a0g     <- noisy_est("binary", pY_0a0g)
    
    pM_1a0g     <- G$M_1
    pTV1_1a0g   <- predict(mTV1_1, newdata = data.frame(                                                           M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_1a0g   <- noisy_est("binary", pTV1_1a0g)
    pTV2_1a0g   <- predict(mTV2_1, newdata = data.frame(                                        TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_1a0g   <- noisy_est("binary", pTV2_1a0g)
    pTV3_1a0g   <- predict(mTV3_1, newdata = data.frame(                     TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_1a0g   <- noisy_est("binary", pTV3_1a0g)
    pY_1a0g     <- predict(mY_1, newdata = data.frame(    TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, Y_0 = pY_0a0g, TV3_0 = pTV3_0a0g, TV2_0 = pTV2_0a0g, TV1_0 = pTV1_0a0g, M_0 = pM_0a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_1a0g     <- noisy_est("binary", pY_1a0g)
    
    pM_2a0g     <- G$M_2
    pTV1_2a0g   <- predict(mTV1_2, newdata = data.frame(                                                           M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_2a0g   <- noisy_est("binary", pTV1_2a0g)
    pTV2_2a0g   <- predict(mTV2_2, newdata = data.frame(                                        TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_2a0g   <- noisy_est("binary", pTV2_2a0g)
    pTV3_2a0g   <- predict(mTV3_2, newdata = data.frame(                     TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_2a0g   <- noisy_est("binary", pTV3_2a0g)
    pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV3_1 = pTV3_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_2a0g     <- noisy_est("binary", pY_2a0g)
    
    pM_3a0g     <- G$M_3
    pTV1_3a0g   <- predict(mTV1_3, newdata = data.frame(                                                           M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_3a0g   <- noisy_est("binary", pTV1_3a0g)
    pTV2_3a0g   <- predict(mTV2_3, newdata = data.frame(                                        TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_3a0g   <- noisy_est("binary", pTV2_3a0g)
    pTV3_3a0g   <- predict(mTV3_3, newdata = data.frame(                     TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_3a0g   <- noisy_est("binary", pTV3_3a0g)
    pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV3_2 = pTV3_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_3a0g     <- noisy_est("binary", pY_3a0g)
    
    pM_4a0g     <- G$M_4
    pTV1_4a0g   <- predict(mTV1_4, newdata = data.frame(                                                           M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_4a0g   <- noisy_est("binary", pTV1_4a0g)
    pTV2_4a0g   <- predict(mTV2_4, newdata = data.frame(                                        TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_4a0g   <- noisy_est("binary", pTV2_4a0g)
    pTV3_4a0g   <- predict(mTV3_4, newdata = data.frame(                     TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_4a0g   <- noisy_est("binary", pTV3_4a0g)
    pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV3_3 = pTV3_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_4a0g     <- noisy_est("binary", pY_4a0g)
    
    pM_5a0g     <- G$M_5
    pTV1_5a0g   <- predict(mTV1_5, newdata = data.frame(                                                           M_5 = pM_5a0g, Y_4 = pY_4a0g, TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV1_5a0g   <- noisy_est("binary", pTV1_5a0g)
    pTV2_5a0g   <- predict(mTV2_5, newdata = data.frame(                                        TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV2_5a0g   <- noisy_est("binary", pTV2_5a0g)
    pTV3_5a0g   <- predict(mTV3_5, newdata = data.frame(                     TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pTV3_5a0g   <- noisy_est("binary", pTV3_5a0g)
    
    #3b - Simulate outcome for each i and save mean outcome
    pY_5a1g     <- predict(mY_5, newdata = data.frame(    TV3_5 = pTV3_5a1g, TV2_5 = pTV2_5a1g, TV1_5 = pTV1_5a1g, M_5 = pM_5a1g, Y_4 = pY_4a1g, TV3_4 = pTV3_4a1g, TV2_4 = pTV2_4a1g, TV1_4 = pTV1_4a1g, M_4 = pM_4a1g, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    pY_5a0g     <- predict(mY_5, newdata = data.frame(    TV3_5 = pTV3_5a0g, TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV3_4 = pTV3_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_5a1g <- pY_5a1g
    MC$pY_5a0g <- pY_5a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_5a1g)
      outcomes$Ya0g1 <- mean(MC$pY_5a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_5a1g)
      outcomes$Ya0g0 <- mean(MC$pY_5a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}

#### Complete-Case Analysis - 3 time-point ####
medgf_3_cc <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 3 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                           MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                             A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                         Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                         Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                    MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                  M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                              TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                            TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                            Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                          M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                      TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                            Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                          M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                      TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_2nat overall, by exposure
  MC$Y_nat <- pY_2nat
  MC$YA_nat <- pY_2natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
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
    pY_2a0g     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, Y_1 = pY_1a0g, TV2_1 = pTV2_1a0g, TV1_1 = pTV1_1a0g, M_1 = pM_1a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_2a1g <- pY_2a1g
    MC$pY_2a0g <- pY_2a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_2a1g)
      outcomes$Ya0g1 <- mean(MC$pY_2a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_2a1g)
      outcomes$Ya0g0 <- mean(MC$pY_2a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}



#### Complete-Case Analysis - 4 time-point ####
medgf_4_cc <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 4 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                           MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                             A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                         Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                         Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                         Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                    MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                  M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                              TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                            TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                            Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                          M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                      TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                            Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                          M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                      TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                            Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                          M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                      TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_3nat overall, by exposure
  MC$Y_nat <- pY_3nat
  MC$YA_nat <- pY_3natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
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
    pY_3a0g     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, Y_2 = pY_2a0g, TV2_2 = pTV2_2a0g, TV1_2 = pTV1_2a0g, M_2 = pM_2a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_3a1g <- pY_3a1g
    MC$pY_3a0g <- pY_3a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_3a1g)
      outcomes$Ya0g1 <- mean(MC$pY_3a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_3a1g)
      outcomes$Ya0g0 <- mean(MC$pY_3a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}




#### Complete-Case Analysis - 5 time-point ####
medgf_5_cc <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 5 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                       Y_3 + TV2_3 + TV1_3 + poly(M_3, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3,                           data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                 M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~         TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_4   <- glm(Y_4   ~ TV2_4 + TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  eM_4 <- sqrt(summary(mM_4)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                           MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                             A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                         Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                         Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                         Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3nat     <- noisy_est("binary", pY_3nat)
  
  pM_4nat     <- predict(mM_4, newdata = data.frame(                                                         Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4nat     <- noisy_est("continuous", pM_4nat, eM_4)
  pTV1_4nat   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4nat   <- noisy_est("binary", pTV1_4nat)
  pTV2_4nat   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4nat   <- noisy_est("binary", pTV2_4nat)
  pY_4nat     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                    MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                  M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                              TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                            TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                            Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                          M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                      TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                            Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                          M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                      TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                            Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                          M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                      TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3natA     <- noisy_est("binary", pY_3natA)
  
  pM_4natA     <- predict(mM_4, newdata = data.frame(                                                            Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4natA     <- noisy_est("continuous", pM_4natA, eM_4)
  pTV1_4natA   <- predict(mTV1_4, newdata = data.frame(                                          M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4natA   <- noisy_est("binary", pTV1_4natA)
  pTV2_4natA   <- predict(mTV2_4, newdata = data.frame(                      TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4natA   <- noisy_est("binary", pTV2_4natA)
  pY_4natA     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_4nat overall, by exposure
  MC$Y_nat <- pY_4nat
  MC$YA_nat <- pY_4natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a1     <- noisy_est("continuous", pM_4a1, eM_4)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a0     <- noisy_est("continuous", pM_4a0, eM_4)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
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
    pY_4a0g     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, Y_3 = pY_3a0g, TV2_3 = pTV2_3a0g, TV1_3 = pTV1_3a0g, M_3 = pM_3a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_4a1g <- pY_4a1g
    MC$pY_4a0g <- pY_4a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_4a1g)
      outcomes$Ya0g1 <- mean(MC$pY_4a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_4a1g)
      outcomes$Ya0g0 <- mean(MC$pY_4a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}



#### Complete-Case Analysis - 6 time-point ####
medgf_6_cc <- function(dataset) {
  #Creating noisy function
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(50000, 1, prediction)
    } else {
      noisy <- rnorm(50000, prediction, error)
    }
    return(noisy)
  }
  
  #### 6 time-point marginal model #### 
  #Creating bootstrap resample of original dataset (or keep original sample where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #Creating Step 1 model coefficients based on observed sample
  mA     <- glm(A     ~                                                                C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4,                                                                   data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                                            A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 ,                                                    data = data_boot, family = gaussian("identity"))
  mTV1_0 <- glm(TV1_0 ~                                             M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mTV2_0 <- glm(TV2_0 ~                                     TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mY_0   <- glm(Y_0   ~                             TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~                       Y_0 + TV2_0 + TV1_0 + poly(M_0, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0,                           data = data_boot, family = gaussian("identity"))
  mTV1_1 <- glm(TV1_1 ~                 M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mTV2_1 <- glm(TV2_1 ~         TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mY_1   <- glm(Y_1   ~ TV2_1 + TV1_1 + M_1 + Y_0 + TV2_0 + TV1_0 + M_0          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_0 + C3*M_0 + C4*M_0 + A*M_1 + C3*M_1 + C4*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~                       Y_1 + TV2_1 + TV1_1 + poly(M_1, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1,                           data = data_boot, family = gaussian("identity"))
  mTV1_2 <- glm(TV1_2 ~                 M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mTV2_2 <- glm(TV2_2 ~         TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mY_2   <- glm(Y_2   ~ TV2_2 + TV1_2 + M_2 + Y_1 + TV2_1 + TV1_1 + M_1          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_1 + C3*M_1 + C4*M_1 + A*M_2 + C3*M_2 + C4*M_2, data = data_boot, family = binomial("logit"))
  mM_3   <- glm(M_3   ~                       Y_2 + TV2_2 + TV1_2 + poly(M_2, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2,                           data = data_boot, family = gaussian("identity"))
  mTV1_3 <- glm(TV1_3 ~                 M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mTV2_3 <- glm(TV2_3 ~         TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mY_3   <- glm(Y_3   ~ TV2_3 + TV1_3 + M_3 + Y_2 + TV2_2 + TV1_2 + M_2          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_2 + C3*M_2 + C4*M_2 + A*M_3 + C3*M_3 + C4*M_3, data = data_boot, family = binomial("logit"))
  mM_4   <- glm(M_4   ~                       Y_3 + TV2_3 + TV1_3 + poly(M_3, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3,                           data = data_boot, family = gaussian("identity"))
  mTV1_4 <- glm(TV1_4 ~                 M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mTV2_4 <- glm(TV2_4 ~         TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mY_4   <- glm(Y_4   ~ TV2_4 + TV1_4 + M_4 + Y_3 + TV2_3 + TV1_3 + M_3          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_3 + C3*M_3 + C4*M_3 + A*M_4 + C3*M_4 + C4*M_4, data = data_boot, family = binomial("logit"))
  mM_5   <- glm(M_5   ~                       Y_4 + TV2_4 + TV1_4 + poly(M_4, 2) + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4,                           data = data_boot, family = gaussian("identity"))
  mTV1_5 <- glm(TV1_5 ~                 M_5 + Y_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  mTV2_5 <- glm(TV2_5 ~         TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #1b - parametric model for final Y given measured past
  mY_5   <- glm(Y_5   ~ TV2_5 + TV1_5 + M_5 + Y_4 + TV2_4 + TV1_4 + M_4          + A + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C3*C4 + A*C3 + A*C4 + A*M_4 + C3*M_4 + C4*M_4 + A*M_5 + C3*M_5 + C4*M_5, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - only needed for 'M' models
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  eM_3 <- sqrt(summary(mM_3)$dispersion)
  eM_4 <- sqrt(summary(mM_4)$dispersion)
  eM_5 <- sqrt(summary(mM_5)$dispersion)
  
  #Creating a dataset to store the estimated outcomes for each potential outcome
  outcomes <- data.frame(Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Monte Carlo resampling step based on Naimi g-computation notes
  index <- sample(1:NROW(data_boot), size = 50000, replace = T)
  MC <- data_boot[index, ]
  
  #Natural Course - Simulating causal model and saving outcome to compare against observed reality
  pA_nat      <- predict(mA, newdata = data.frame(                                                                                                                                           MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pA_nat      <- noisy_est("binary", pA_nat)
  pM_0nat     <- predict(mM_0, newdata = data.frame(                                                                                                                             A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0nat     <- noisy_est("continuous", pM_0nat, eM_0)
  pTV1_0nat   <- predict(mTV1_0, newdata = data.frame(                                                                                                            M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0nat   <- noisy_est("binary", pTV1_0nat)
  pTV2_0nat   <- predict(mTV2_0, newdata = data.frame(                                                                                         TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0nat   <- noisy_est("binary", pTV2_0nat)
  pY_0nat     <- predict(mY_0, newdata = data.frame(                                                                        TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0nat     <- noisy_est("binary", pY_0nat)
  
  pM_1nat     <- predict(mM_1, newdata = data.frame(                                                         Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1nat     <- noisy_est("continuous", pM_1nat, eM_1)
  pTV1_1nat   <- predict(mTV1_1, newdata = data.frame(                                        M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1nat   <- noisy_est("binary", pTV1_1nat)
  pTV2_1nat   <- predict(mTV2_1, newdata = data.frame(                     TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1nat   <- noisy_est("binary", pTV2_1nat)
  pY_1nat     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, Y_0 = pY_0nat, TV2_0 = pTV2_0nat, TV1_0 = pTV1_0nat, M_0 = pM_0nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1nat     <- noisy_est("binary", pY_1nat)
  
  pM_2nat     <- predict(mM_2, newdata = data.frame(                                                         Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2nat     <- noisy_est("continuous", pM_2nat, eM_2)
  pTV1_2nat   <- predict(mTV1_2, newdata = data.frame(                                        M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2nat   <- noisy_est("binary", pTV1_2nat)
  pTV2_2nat   <- predict(mTV2_2, newdata = data.frame(                     TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2nat   <- noisy_est("binary", pTV2_2nat)
  pY_2nat     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, Y_1 = pY_1nat, TV2_1 = pTV2_1nat, TV1_1 = pTV1_1nat, M_1 = pM_1nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2nat     <- noisy_est("binary", pY_2nat)
  
  pM_3nat     <- predict(mM_3, newdata = data.frame(                                                         Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3nat     <- noisy_est("continuous", pM_3nat, eM_3)
  pTV1_3nat   <- predict(mTV1_3, newdata = data.frame(                                        M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3nat   <- noisy_est("binary", pTV1_3nat)
  pTV2_3nat   <- predict(mTV2_3, newdata = data.frame(                     TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3nat   <- noisy_est("binary", pTV2_3nat)
  pY_3nat     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, Y_2 = pY_2nat, TV2_2 = pTV2_2nat, TV1_2 = pTV1_2nat, M_2 = pM_2nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3nat     <- noisy_est("binary", pY_3nat)
  
  pM_4nat     <- predict(mM_4, newdata = data.frame(                                                         Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4nat     <- noisy_est("continuous", pM_4nat, eM_4)
  pTV1_4nat   <- predict(mTV1_4, newdata = data.frame(                                        M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4nat   <- noisy_est("binary", pTV1_4nat)
  pTV2_4nat   <- predict(mTV2_4, newdata = data.frame(                     TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4nat   <- noisy_est("binary", pTV2_4nat)
  pY_4nat     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, Y_3 = pY_3nat, TV2_3 = pTV2_3nat, TV1_3 = pTV1_3nat, M_3 = pM_3nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4nat     <- noisy_est("binary", pY_4nat)
  
  pM_5nat     <- predict(mM_5, newdata = data.frame(                                                         Y_4 = pY_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5nat     <- noisy_est("continuous", pM_5nat, eM_5)
  pTV1_5nat   <- predict(mTV1_5, newdata = data.frame(                                        M_5 = pM_5nat, Y_4 = pY_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5nat   <- noisy_est("binary", pTV1_5nat)
  pTV2_5nat   <- predict(mTV2_5, newdata = data.frame(                     TV1_5 = pTV1_5nat, M_5 = pM_5nat, Y_4 = pY_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5nat   <- noisy_est("binary", pTV2_5nat)
  pY_5nat     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5nat, TV1_5 = pTV1_5nat, M_5 = pM_5nat, Y_4 = pY_4nat, TV2_4 = pTV2_4nat, TV1_4 = pTV1_4nat, M_4 = pM_4nat, A = pA_nat, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Natural Course - Simulating causal model post-exposure and saving outcome to compare against observed reality
  pM_0natA     <- predict(mM_0, newdata = data.frame(                                                                                                                                    MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0natA     <- noisy_est("continuous", pM_0natA, eM_0)
  pTV1_0natA   <- predict(mTV1_0, newdata = data.frame(                                                                                                                  M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0natA   <- noisy_est("binary", pTV1_0natA)
  pTV2_0natA   <- predict(mTV2_0, newdata = data.frame(                                                                                              TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0natA   <- noisy_est("binary", pTV2_0natA)
  pY_0natA     <- predict(mY_0, newdata = data.frame(                                                                            TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0natA     <- noisy_est("binary", pY_0natA)
  
  pM_1natA     <- predict(mM_1, newdata = data.frame(                                                            Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1natA     <- noisy_est("continuous", pM_1natA, eM_1)
  pTV1_1natA   <- predict(mTV1_1, newdata = data.frame(                                          M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1natA   <- noisy_est("binary", pTV1_1natA)
  pTV2_1natA   <- predict(mTV2_1, newdata = data.frame(                      TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1natA   <- noisy_est("binary", pTV2_1natA)
  pY_1natA     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, Y_0 = pY_0natA, TV2_0 = pTV2_0natA, TV1_0 = pTV1_0natA, M_0 = pM_0natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1natA     <- noisy_est("binary", pY_1natA)
  
  pM_2natA     <- predict(mM_2, newdata = data.frame(                                                            Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2natA     <- noisy_est("continuous", pM_2natA, eM_2)
  pTV1_2natA   <- predict(mTV1_2, newdata = data.frame(                                          M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2natA   <- noisy_est("binary", pTV1_2natA)
  pTV2_2natA   <- predict(mTV2_2, newdata = data.frame(                      TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2natA   <- noisy_est("binary", pTV2_2natA)
  pY_2natA     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, Y_1 = pY_1natA, TV2_1 = pTV2_1natA, TV1_1 = pTV1_1natA, M_1 = pM_1natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2natA     <- noisy_est("binary", pY_2natA)
  
  pM_3natA     <- predict(mM_3, newdata = data.frame(                                                            Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3natA     <- noisy_est("continuous", pM_3natA, eM_3)
  pTV1_3natA   <- predict(mTV1_3, newdata = data.frame(                                          M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3natA   <- noisy_est("binary", pTV1_3natA)
  pTV2_3natA   <- predict(mTV2_3, newdata = data.frame(                      TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3natA   <- noisy_est("binary", pTV2_3natA)
  pY_3natA     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, Y_2 = pY_2natA, TV2_2 = pTV2_2natA, TV1_2 = pTV1_2natA, M_2 = pM_2natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3natA     <- noisy_est("binary", pY_3natA)
  
  pM_4natA     <- predict(mM_4, newdata = data.frame(                                                            Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4natA     <- noisy_est("continuous", pM_4natA, eM_4)
  pTV1_4natA   <- predict(mTV1_4, newdata = data.frame(                                          M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4natA   <- noisy_est("binary", pTV1_4natA)
  pTV2_4natA   <- predict(mTV2_4, newdata = data.frame(                      TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4natA   <- noisy_est("binary", pTV2_4natA)
  pY_4natA     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, Y_3 = pY_3natA, TV2_3 = pTV2_3natA, TV1_3 = pTV1_3natA, M_3 = pM_3natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4natA     <- noisy_est("binary", pY_4natA)
  
  pM_5natA     <- predict(mM_5, newdata = data.frame(                                                            Y_4 = pY_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5natA     <- noisy_est("continuous", pM_5natA, eM_5)
  pTV1_5natA   <- predict(mTV1_5, newdata = data.frame(                                          M_5 = pM_5natA, Y_4 = pY_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_5natA   <- noisy_est("binary", pTV1_5natA)
  pTV2_5natA   <- predict(mTV2_5, newdata = data.frame(                      TV1_5 = pTV1_5natA, M_5 = pM_5natA, Y_4 = pY_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_5natA   <- noisy_est("binary", pTV2_5natA)
  pY_5natA     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5natA, TV1_5 = pTV1_5natA, M_5 = pM_5natA, Y_4 = pY_4natA, TV2_4 = pTV2_4natA, TV1_4 = pTV1_4natA, M_4 = pM_4natA, MC[,c("A","C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  
  #Saving pY_4nat overall, by exposure
  MC$Y_nat <- pY_5nat
  MC$YA_nat <- pY_5natA
  
  outcomes$Y_nat <- mean(MC$Y_nat)
  outcomes$Ya0_nat <- mean(MC$Y_nat[MC$A == 0])
  outcomes$Ya1_nat <- mean(MC$Y_nat[MC$A == 1])
  
  outcomes$YA_nat <- mean(MC$YA_nat)
  outcomes$YA0_nat <- mean(MC$YA_nat[MC$A == 0])
  outcomes$YA1_nat <- mean(MC$YA_nat[MC$A == 1])
  
  #Step 2 - Estimate joint distribution of tv-M given A, A*
  #2a - Set baseline C to observed for each i at each time
  #Fixing A=1
  pM_0a1     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a1     <- noisy_est("continuous", pM_0a1, eM_0)
  pTV1_0a1   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a1   <- noisy_est("binary", pTV1_0a1)
  pTV2_0a1   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a1   <- noisy_est("binary", pTV2_0a1)
  pY_0a1     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a1     <- noisy_est("binary", pY_0a1)
  
  pM_1a1     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a1     <- noisy_est("continuous", pM_1a1, eM_1)
  pTV1_1a1   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a1   <- noisy_est("binary", pTV1_1a1)
  pTV2_1a1   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a1   <- noisy_est("binary", pTV2_1a1)
  pY_1a1     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, Y_0 = pY_0a1, TV2_0 = pTV2_0a1, TV1_0 = pTV1_0a1, M_0 = pM_0a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a1     <- noisy_est("binary", pY_1a1)
  
  pM_2a1     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a1     <- noisy_est("continuous", pM_2a1, eM_2)
  pTV1_2a1   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a1   <- noisy_est("binary", pTV1_2a1)
  pTV2_2a1   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a1   <- noisy_est("binary", pTV2_2a1)
  pY_2a1     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, Y_1 = pY_1a1, TV2_1 = pTV2_1a1, TV1_1 = pTV1_1a1, M_1 = pM_1a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a1     <- noisy_est("binary", pY_2a1)
  
  pM_3a1     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a1     <- noisy_est("continuous", pM_3a1, eM_3)
  pTV1_3a1   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a1   <- noisy_est("binary", pTV1_3a1)
  pTV2_3a1   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a1   <- noisy_est("binary", pTV2_3a1)
  pY_3a1     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, Y_2 = pY_2a1, TV2_2 = pTV2_2a1, TV1_2 = pTV1_2a1, M_2 = pM_2a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a1     <- noisy_est("binary", pY_3a1)
  
  pM_4a1     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a1     <- noisy_est("continuous", pM_4a1, eM_4)
  pTV1_4a1   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a1   <- noisy_est("binary", pTV1_4a1)
  pTV2_4a1   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a1   <- noisy_est("binary", pTV2_4a1)
  pY_4a1     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, Y_3 = pY_3a1, TV2_3 = pTV2_3a1, TV1_3 = pTV1_3a1, M_3 = pM_3a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a1     <- noisy_est("binary", pY_4a1)
  
  pM_5a1     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a1, TV2_4 = pTV2_4a1, TV1_4 = pTV1_4a1, M_4 = pM_4a1, A = 1, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5a1     <- noisy_est("continuous", pM_5a1, eM_5)
  
  #2b - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a1[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_5 = pM_5a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #2c - Repeat 2a assigning A=0, with M at each t predicted based on A=0
  #Fixing A=0
  pM_0a0     <- predict(mM_0, newdata = data.frame(                                                                                                                      A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_0a0     <- noisy_est("continuous", pM_0a0, eM_0)
  pTV1_0a0   <- predict(mTV1_0, newdata = data.frame(                                                                                                      M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_0a0   <- noisy_est("binary", pTV1_0a0)
  pTV2_0a0   <- predict(mTV2_0, newdata = data.frame(                                                                                    TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_0a0   <- noisy_est("binary", pTV2_0a0)
  pY_0a0     <- predict(mY_0, newdata = data.frame(                                                                    TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_0a0     <- noisy_est("binary", pY_0a0)
  
  pM_1a0     <- predict(mM_1, newdata = data.frame(                                                      Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_1a0     <- noisy_est("continuous", pM_1a0, eM_1)
  pTV1_1a0   <- predict(mTV1_1, newdata = data.frame(                                      M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_1a0   <- noisy_est("binary", pTV1_1a0)
  pTV2_1a0   <- predict(mTV2_1, newdata = data.frame(                    TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_1a0   <- noisy_est("binary", pTV2_1a0)
  pY_1a0     <- predict(mY_1, newdata = data.frame(    TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, Y_0 = pY_0a0, TV2_0 = pTV2_0a0, TV1_0 = pTV1_0a0, M_0 = pM_0a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_1a0     <- noisy_est("binary", pY_1a0)
  
  pM_2a0     <- predict(mM_2, newdata = data.frame(                                                      Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_2a0     <- noisy_est("continuous", pM_2a0, eM_2)
  pTV1_2a0   <- predict(mTV1_2, newdata = data.frame(                                      M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_2a0   <- noisy_est("binary", pTV1_2a0)
  pTV2_2a0   <- predict(mTV2_2, newdata = data.frame(                    TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_2a0   <- noisy_est("binary", pTV2_2a0)
  pY_2a0     <- predict(mY_2, newdata = data.frame(    TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, Y_1 = pY_1a0, TV2_1 = pTV2_1a0, TV1_1 = pTV1_1a0, M_1 = pM_1a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_2a0     <- noisy_est("binary", pY_2a0)
  
  pM_3a0     <- predict(mM_3, newdata = data.frame(                                                      Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_3a0     <- noisy_est("continuous", pM_3a0, eM_3)
  pTV1_3a0   <- predict(mTV1_3, newdata = data.frame(                                      M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_3a0   <- noisy_est("binary", pTV1_3a0)
  pTV2_3a0   <- predict(mTV2_3, newdata = data.frame(                    TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_3a0   <- noisy_est("binary", pTV2_3a0)
  pY_3a0     <- predict(mY_3, newdata = data.frame(    TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, Y_2 = pY_2a0, TV2_2 = pTV2_2a0, TV1_2 = pTV1_2a0, M_2 = pM_2a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_3a0     <- noisy_est("binary", pY_3a0)
  
  pM_4a0     <- predict(mM_4, newdata = data.frame(                                                      Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_4a0     <- noisy_est("continuous", pM_4a0, eM_4)
  pTV1_4a0   <- predict(mTV1_4, newdata = data.frame(                                      M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV1_4a0   <- noisy_est("binary", pTV1_4a0)
  pTV2_4a0   <- predict(mTV2_4, newdata = data.frame(                    TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pTV2_4a0   <- noisy_est("binary", pTV2_4a0)
  pY_4a0     <- predict(mY_4, newdata = data.frame(    TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, Y_3 = pY_3a0, TV2_3 = pTV2_3a0, TV1_3 = pTV1_3a0, M_3 = pM_3a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pY_4a0     <- noisy_est("binary", pY_4a0)
  
  pM_5a0     <- predict(mM_5, newdata = data.frame(                                                      Y_4 = pY_4a0, TV2_4 = pTV2_4a0, TV1_4 = pTV1_4a0, M_4 = pM_4a0, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
  pM_5a0     <- noisy_est("continuous", pM_5a0, eM_5)
  
  #2d - Randomly permute the n values of joint mediators under A=1 as in 2a (above), save for below
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_3 = pM_3a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_4 = pM_4a0[sample(1:NROW(MC), NROW(MC), replace = F)],
                           M_5 = pM_5a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 3 - Estimate Q(a,a*), Q(a*,a), Q(a,a), Q(a*,a*) using different 2b/d M for Ga/ Ga*
  #Note - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #       As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  for (a in 0:1) {
    #G changes permuted M
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
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
    pY_5a0g     <- predict(mY_5, newdata = data.frame(    TV2_5 = pTV2_5a0g, TV1_5 = pTV1_5a0g, M_5 = pM_5a0g, Y_4 = pY_4a0g, TV2_4 = pTV2_4a0g, TV1_4 = pTV1_4a0g, M_4 = pM_4a0g, A = 0, MC[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9")]), type = "response")
    
    MC$pY_5a1g <- pY_5a1g
    MC$pY_5a0g <- pY_5a0g
    
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_5a1g)
      outcomes$Ya0g1 <- mean(MC$pY_5a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_5a1g)
      outcomes$Ya0g0 <- mean(MC$pY_5a0g)
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
  results <- data.frame(rep        = boot,
                        
                        Y_nat      = mean(outcomes$Y_nat),          YA_nat     = mean(outcomes$YA_nat),
                        Ya0_nat    = mean(outcomes$Ya0_nat),        Ya1_nat    = mean(outcomes$Ya1_nat),
                        YA0_nat    = mean(outcomes$YA0_nat),        YA1_nat    = mean(outcomes$YA1_nat),
                        
                        Ya0g0      = mean(outcomes$Ya0g0),          Ya0g1      = mean(outcomes$Ya0g1),
                        Ya1g0      = mean(outcomes$Ya1g0),          Ya1g1      = mean(outcomes$Ya1g1),
                        
                        TE         = mean(outcomes$TE),             TE_m       = mean(outcomes$TE_m),
                        pNDE       = mean(outcomes$pNDE),           tNDE       = mean(outcomes$tNDE),
                        pNIE       = mean(outcomes$pNIE),           tNIE       = mean(outcomes$tNIE),
                        pNDE_m     = mean(outcomes$pNDE_m),         tNDE_m     = mean(outcomes$tNDE_m),
                        pNIE_m     = mean(outcomes$pNIE_m),         tNIE_m     = mean(outcomes$tNIE_m))
  return(results)
}


