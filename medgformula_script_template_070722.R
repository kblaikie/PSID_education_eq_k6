#### Mediational g-formula template script #####
#Author - Kieran Blaikie
#Date - 07 July 2022

#Overview - This script provides a simplified template for the mediational g-formula if ran using wide-format data.
#           Running models in wide-format avoids dealing with within-person dependence between observations
#         - Each step in this script follows the procedure laid out in the 2017 paper by Lin et al below:
#           (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/) 
#         - For simplicity, we use:
#            - Time-varying Exposure, A (binary)
#            - Time-varying Mediator, M (continuous)
#            - Time-varying Outcome, Y (binary)
#            - 2 time-fixed baseline confounders, C1, C2
#            - 1 time-varying confounder, TV1 (binary)
#         - We assume the temporal order: C > A(t) > TV1(t) > M(t) > Y(t) > ... > Y(T)
#         - We have three time-points in this analysis: Baseline (0), Intermediate (1), Final (2)
#            - Each time-varying variable will have the suffix _N indicating the time-point it corresponds to

#Changes - Compared to the medgformula_script_template_07122021.R file, this script fixes prior errors in
#          1) correctly permuting mediator values over time, 2) retaining a g-computation run when using the 
#          original data (where 'boot' = 0 in script line 306, 3) removing the Step 3d K-repetition loop given
#          it isn't necessary to account for simulation error when also using an MC sample, 4) saving continuous
#          variable prediction error terms before looping for efficiency in lines 83-86, and 5) adding Natural Course
#          simulations to assess for gross model misspecification (indicated where Y_obs and Y_nat are dissimilar)

#### Loading necessary packages ####
library(parallel); library(doParallel)

#### Saving mediational g-formula process as function ####
medgf <- function(dataset) {
  #### Creating a function to add noise to predictions for each i and return 0/1 values from predicted probabilities for binary variables ####
  #   - In Steps 2 and 3, the predict() function is used extensively to simulate person-specific variable values for time-varying exposure, mediator, and outcome
  #   - Via the predict() function, each i within covariate strata are assigned the same predicted mean (or response probability)
  #   - To reflect prediction uncertainty for these individuals, we add noise based on 1) the overall prediction uncertainty for continuous variables,
  #     or 2) the strata-specific prediction uncertainty for binary variables 
  #   - Note this function doesn't currently accommodate categorical variables or non-Normal continuous variables, though this could be in theory be accounted for
  noisy_est <- function(type, prediction, error) {
    if (type == "binary") {
      noisy <- rbinom(10000,1,prediction)
    } else {
      noisy <- rnorm(10000, prediction, error)
    }
    return(noisy)
  }
  
  #Creating bootstrap resample of original dataset
  #Creating bootstrap resample of original dataset, keeping original where boot = 0)
  if (boot == 0) {
    index_boot <- 1:NROW(dataset)
  } else {
    index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
  }
  data_boot <- dataset[index_boot, ]
  
  #### Step 1 - Fit parametric models for the observed data ####
  #           - Each forward prediction in Steps 2 and 3 are done using these fitted models. This ensures each prediction maintains the
  #             variance of our observed data, even if we make preductions using expanded monte carlo resamples for simulation stability
  #           - Each model should adjusting for all covariates occuring prior to it in the assumed causal DAG
  #           - For computational simplicity, we make a first-order Markov assumption, assuming each factor is independent 
  #             of all covariates at times <t-1 except baseline confounders, conditional on all covariates at times t and t-1 
  #           - For each model, we include an interaction term between the most recent A and M variables prior to the model outcome at time t 
  
  #Step 1a - For times t>=0, fit parametric models for the joint density of C, A, and M at t given measured past
  mA_0   <- glm(A_0   ~                                               C1 + C2,           data = data_boot, family = binomial("logit"))
  mTV1_0 <- glm(TV1_0 ~                                         A_0 + C1 + C2,           data = data_boot, family = binomial("logit"))
  mM_0   <- glm(M_0   ~                                 TV1_0 + A_0 + C1 + C2,           data = data_boot, family = gaussian("identity"))
  mY_0   <- glm(Y_0   ~                           M_0 + TV1_0 + A_0 + C1 + C2 + A_0*M_0, data = data_boot, family = binomial("logit"))
  mA_1   <- glm(A_1   ~                     Y_0 + M_0 + TV1_0 + A_0 + C1 + C2 + A_0*M_0, data = data_boot, family = binomial("logit"))
  mTV1_1 <- glm(TV1_1 ~               A_1 + Y_0 + M_0 + TV1_0 + A_0 + C1 + C2 + A_1*M_0, data = data_boot, family = binomial("logit"))
  mM_1   <- glm(M_1   ~       TV1_1 + A_1 + Y_0 + M_0 + TV1_0 + A_0 + C1 + C2 + A_1*M_0, data = data_boot, family = gaussian("identity"))
  mY_1   <- glm(Y_1   ~ M_1 + TV1_1 + A_1 + Y_0 + M_0 + TV1_0 + A_0 + C1 + C2 + A_1*M_1, data = data_boot, family = binomial("logit"))
  mA_2   <- glm(A_2   ~                     Y_1 + M_1 + TV1_1 + A_1 + C1 + C2 + A_1*M_1, data = data_boot, family = binomial("logit"))
  mTV1_2 <- glm(TV1_2 ~               A_2 + Y_1 + M_1 + TV1_1 + A_1 + C1 + C2 + A_2*M_1, data = data_boot, family = binomial("logit"))
  mM_2   <- glm(M_2   ~       TV1_2 + A_2 + Y_1 + M_1 + TV1_1 + A_1 + C1 + C2 + A_2*M_1, data = data_boot, family = gaussian("identity"))
  
  #Step 1b - Fit a parametric model for the mean outcome at end of follow-up given the measured past
  mY_2   <- glm(Y_2   ~ M_2 + TV1_2 + A_2 + Y_1 + M_1 + TV1_1 + A_1 + C1 + C2 + A_2*M_2, data = data_boot, family = binomial("logit"))
  
  #Storing model error terms for faster computation - needed for continuous var models (only 'M' in this case)
  eM_0 <- sqrt(summary(mM_0)$dispersion)
  eM_1 <- sqrt(summary(mM_1)$dispersion)
  eM_2 <- sqrt(summary(mM_2)$dispersion)
  
  #For Steps 2-3, first create a dataset to store our potential outcomes from Step 3c as well as Natural Course predictions
  outcomes <- data.frame(Y_obs = NA, Y_nat = NA, Ya0g0 = NA, Ya0g1 = NA, Ya1g0 = NA, Ya1g1 = NA)
  
  #Next resample with replacement from our original data to ensure we obtain an accurate simulation/well-behaved monte carlo estimator
  # - Note that by using an expanded Monte Carlo sample, Step 3d (the 'k-loop') is no longer required under the assumption
  #   that our MC sample is sufficiently large to make simulation error minimal
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #### Step 2 - Estimate the joint distribution of time-varying mediators under time-varying exposure interventions a, a* ####
  #           - Step 2 is solely to obtain a set of mediator values we would expect to see under each exposure intervention,
  #             that we can then permute and use in Step 3 for each joint A-M intervention
  
  #Step 2a - Set baseline covariates to the observed values for subject i, then recursively for each time t from 0-T-1:
  #          i  - Generate time t TV, A, M based on model coefficients generated in Step 1
  #          ii - Assign time t exposure under the intervention A=a
  #        - Note no need to predict A_N except for the Natural Course, given we are assigning A = a for each person  
  
  #Natural Course
  pA_0nat      <- predict(mA_0, newdata = data.frame(                                                                                                                 MC[,c("C1","C2")]), type = "response")
  pA_0nat      <- noisy_est("binary", pA_0nat)
  pTV1_0nat    <- predict(mTV1_0, newdata = data.frame(                                                                                                A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pTV1_0nat    <- noisy_est("binary", pTV1_0nat)
  pM_0nat      <- predict(mM_0, newdata = data.frame(                                                                               TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pM_0nat      <- noisy_est("continuous", pM_0nat, eM_0)
  pY_0nat      <- predict(mY_0, newdata = data.frame(                                                                M_0 = pM_0nat, TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pY_0nat      <- noisy_est("binary", pY_0nat)
  
  pA_1nat      <- predict(mA_1, newdata = data.frame(                                                 Y_0 = pY_0nat, M_0 = pM_0nat, TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pA_1nat      <- noisy_est("binary", pA_1nat)
  pTV1_1nat    <- predict(mTV1_1, newdata = data.frame(                                A_1 = pA_1nat, Y_0 = pY_0nat, M_0 = pM_0nat, TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pTV1_1nat    <- noisy_est("binary", pTV1_1nat)
  pM_1nat      <- predict(mM_1, newdata = data.frame(               TV1_1 = pTV1_1nat, A_1 = pA_1nat, Y_0 = pY_0nat, M_0 = pM_0nat, TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pM_1nat      <- noisy_est("continuous", pM_1nat, eM_1)
  pY_1nat      <- predict(mY_1, newdata = data.frame(M_1 = pM_1nat, TV1_1 = pTV1_1nat, A_1 = pA_1nat, Y_0 = pY_0nat, M_0 = pM_0nat, TV1_0 = pTV1_0nat, A_0 = pA_0nat, MC[,c("C1","C2")]), type = "response")
  pY_1nat      <- noisy_est("binary", pY_1nat)
  
  pA_2nat      <- predict(mA_2, newdata = data.frame(                                                 Y_1 = pY_1nat, M_1 = pM_1nat, TV1_1 = pTV1_0nat, A_1 = pA_1nat, MC[,c("C1","C2")]), type = "response")
  pA_2nat      <- noisy_est("binary", pA_2nat)
  pTV1_2nat    <- predict(mTV1_2, newdata = data.frame(                                A_2 = pA_2nat, Y_1 = pY_1nat, M_1 = pM_1nat, TV1_1 = pTV1_0nat, A_1 = pA_1nat, MC[,c("C1","C2")]), type = "response")
  pTV1_2nat    <- noisy_est("binary", pTV1_2nat)
  pM_2nat      <- predict(mM_2, newdata = data.frame(               TV1_2 = pTV1_2nat, A_2 = pA_2nat, Y_1 = pY_1nat, M_1 = pM_1nat, TV1_1 = pTV1_1nat, A_1 = pA_1nat, MC[,c("C1","C2")]), type = "response")
  pM_2nat      <- noisy_est("continuous", pM_2nat, eM_2)
  pY_2nat      <- predict(mY_2, newdata = data.frame(M_2 = pM_2nat, TV1_2 = pTV1_2nat, A_2 = pA_2nat, Y_1 = pY_1nat, M_1 = pM_1nat, TV1_1 = pTV1_1nat, A_1 = pA_1nat, MC[,c("C1","C2")]), type = "response")
  
  #Storing the Natural Course prediction from Step 2a (pY_2nat) as well as the individual-specific observed Y_2 for comparison
  MC$pY_2nat <- pY_2nat
  outcomes$Y_nat <- mean(MC$pY_2nat)
  outcomes$Y_obs <- mean(MC$Y_2)
  
  #Fixing a=1
  pTV1_0a1     <- predict(mTV1_0, newdata = data.frame(                                                                                     A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_0a1     <- noisy_est("binary", pTV1_0a1)
  pM_0a1       <- predict(mM_0, newdata = data.frame(                                                                     TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pM_0a1       <- noisy_est("continuous", pM_0a1, eM_0)
  pY_0a1       <- predict(mY_0, newdata = data.frame(                                                       M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pY_0a1       <- noisy_est("binary", pY_0a1)
  
  pTV1_1a1     <- predict(mTV1_1, newdata = data.frame(                              A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_1a1     <- noisy_est("binary", pTV1_1a1)
  pM_1a1       <- predict(mM_1, newdata = data.frame(              TV1_1 = pTV1_1a1, A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pM_1a1       <- noisy_est("continuous", pM_1a1, eM_1)
  pY_1a1       <- predict(mY_1, newdata = data.frame(M_1 = pM_1a1, TV1_1 = pTV1_1a1, A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pY_1a1       <- noisy_est("binary", pY_1a1)
  
  pTV1_2a1     <- predict(mTV1_2, newdata = data.frame(                              A_2 = 1, Y_1 = pY_1a1, M_1 = pM_1a1, TV1_1 = pTV1_0a1, A_1 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_2a1     <- noisy_est("binary", pTV1_2a1)
  pM_2a1       <- predict(mM_2, newdata = data.frame(              TV1_2 = pTV1_2a1, A_2 = 1, Y_1 = pY_1a1, M_1 = pM_1a1, TV1_1 = pTV1_1a1, A_1 = 1, MC[,c("C1","C2")]), type = "response")
  pM_2a1       <- noisy_est("continuous", pM_2a1, eM_2)
  
  #Step 2b - Randomly permute the n values of joint mediators under A=1, provided in 2a above. Save for Step 3
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a1[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a1[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #Step 2c - Repeat Step 2a assigning A=0, with M at each time t predicted based on A=0
  
  #Fixing a=0
  pTV1_0a0     <- predict(mTV1_0, newdata = data.frame(                                                                                     A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_0a0     <- noisy_est("binary", pTV1_0a0)
  pM_0a0       <- predict(mM_0, newdata = data.frame(                                                                     TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pM_0a0       <- noisy_est("continuous", pM_0a0, eM_0)
  pY_0a0       <- predict(mY_0, newdata = data.frame(                                                       M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pY_0a0       <- noisy_est("binary", pY_0a0)
  
  pTV1_1a0     <- predict(mTV1_1, newdata = data.frame(                              A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_1a0     <- noisy_est("binary", pTV1_1a0)
  pM_1a0       <- predict(mM_1, newdata = data.frame(              TV1_1 = pTV1_1a0, A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pM_1a0       <- noisy_est("continuous", pM_1a0, eM_1)
  pY_1a0       <- predict(mY_1, newdata = data.frame(M_1 = pM_1a0, TV1_1 = pTV1_1a0, A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pY_1a0       <- noisy_est("binary", pY_1a0)
  
  pTV1_2a0     <- predict(mTV1_2, newdata = data.frame(                              A_2 = 0, Y_1 = pY_1a0, M_1 = pM_1a0, TV1_1 = pTV1_0a0, A_1 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_2a0     <- noisy_est("binary", pTV1_2a0)
  pM_2a0       <- predict(mM_2, newdata = data.frame(              TV1_2 = pTV1_2a0, A_2 = 0, Y_1 = pY_1a0, M_1 = pM_1a0, TV1_1 = pTV1_1a0, A_1 = 0, MC[,c("C1","C2")]), type = "response")
  pM_2a0       <- noisy_est("continuous", pM_2a0, eM_2)
  
  #Step 2d - Repeat Step 2b storing the randomly permuted predicted mediator for each t generated in Step 2c
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_1 = pM_1a0[sample(1:NROW(MC), NROW(MC), replace = F)], 
                           M_2 = pM_2a0[sample(1:NROW(MC), NROW(MC), replace = F)])
  
  #### Step 3 - Estimate Q(a,a), Q(a,a*), Q(a*,a), Q(a*,a*) by repeating the below steps for each exposure and mediator assignment ####
  #           - In this script, this is accomplished by looping through A- 0:1, where for each A, the mediator 'G' is assigned as under that exposure intervention
  #           - The exposure is separately assigned within each loop in order to obtain each A,M, A,M*, A*,M, and A*,M* expected mean outcome Y at end of follow-up
  #           - This step originally involved looping through 30 times (Step 3d) in order to account for simulation error
  #             As a large MC sample is used, this step of repermuting the mediator isn't necessary 
  #Step 3d - Repeat Step 3a-c K (30) times using different random permutations of joint mediators M stored in Steps 2b and 2d
  #        - In this script, Step 3d is performed through creating a loop for Steps 3a-c, hence why it appears beforehand
  
  #Step 3a - Recursively for each time t from 0 -to T-1:
  #          i   - Repeat Step 2ai replacing the time-varying exposure intervention a with the joint time-varying exposure and mediator intervention a-m 
  #          ii  - Assign the time t mediator as the ith component of the permuted vector for time t from 2b or 2d
  #          iii - Assign time t exposure under the intervention a or a*
  #        - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
  
  for (a in 0:1) {
    #Step 2b/2d random permutation of joint mediator distribution K times, stored as intervention 'G' 
    if (a == 1) {
      G <- perm_Mt_a1
    }
    if (a == 0) {
      G <- perm_Mt_a0
    }
    
    #Q(a,a) or Q(a,a*)
    pTV1_0a1g    <- predict(mTV1_0, newdata = data.frame(                                                                                          A_0 = 1, MC[,c("C1","C2")]), type = "response")
    pTV1_0a1g    <- noisy_est("binary", pTV1_0a1g)
    pM_0a1g      <- G$M_0
    pY_0a1g      <- predict(mY_0, newdata = data.frame(                                                          M_0 = pM_0a1g, TV1_0 = pTV1_0a1g, A_0 = 1, MC[,c("C1","C2")]), type = "response")
    pY_0a1g      <- noisy_est("binary", pY_0a1g)
    
    pTV1_1a1g    <- predict(mTV1_1, newdata = data.frame(                                A_1 = 1, Y_0 = pY_0a1g, M_0 = pM_0a1g, TV1_0 = pTV1_0a1g, A_0 = 1, MC[,c("C1","C2")]), type = "response")
    pTV1_1a1g    <- noisy_est("binary", pTV1_1a1g)
    pM_1a1g      <- G$M_1
    pY_1a1g      <- predict(mY_1, newdata = data.frame(M_1 = pM_1a1g, TV1_1 = pTV1_1a1g, A_1 = 1, Y_0 = pY_0a1g, M_0 = pM_0a1g, TV1_0 = pTV1_0a1g, A_0 = 1, MC[,c("C1","C2")]), type = "response")
    pY_1a1g      <- noisy_est("binary", pY_1a1g)
    
    pTV1_2a1g    <- predict(mTV1_2, newdata = data.frame(                                A_2 = 1, Y_1 = pY_1a1g, M_1 = pM_1a1g, TV1_1 = pTV1_0a1g, A_1 = 1, MC[,c("C1","C2")]), type = "response")
    pTV1_2a1g    <- noisy_est("binary", pTV1_2a1g)
    pM_2a1g      <- G$M_2
    
    #Q(a*,a) or Q(a*,a*)
    pTV1_0a0g    <- predict(mTV1_0, newdata = data.frame(                                                                                          A_0 = 0, MC[,c("C1","C2")]), type = "response")
    pTV1_0a0g    <- noisy_est("binary", pTV1_0a0g)
    pM_0a0g      <- G$M_0
    pY_0a0g      <- predict(mY_0, newdata = data.frame(                                                          M_0 = pM_0a0g, TV1_0 = pTV1_0a0g, A_0 = 0, MC[,c("C1","C2")]), type = "response")
    pY_0a0g      <- noisy_est("binary", pY_0a0g)
    
    pTV1_1a0g    <- predict(mTV1_1, newdata = data.frame(                                A_1 = 0, Y_0 = pY_0a0g, M_0 = pM_0a0g, TV1_0 = pTV1_0a0g, A_0 = 0, MC[,c("C1","C2")]), type = "response")
    pTV1_1a0g    <- noisy_est("binary", pTV1_1a0g)
    pM_1a0g      <- G$M_1
    pY_1a0g      <- predict(mY_1, newdata = data.frame(M_1 = pM_1a0g, TV1_1 = pTV1_1a0g, A_1 = 0, Y_0 = pY_0a0g, M_0 = pM_0a0g, TV1_0 = pTV1_0a0g, A_0 = 0, MC[,c("C1","C2")]), type = "response")
    pY_1a0g      <- noisy_est("binary", pY_1a0g)
    
    pTV1_2a0g    <- predict(mTV1_2, newdata = data.frame(                                A_2 = 0, Y_1 = pY_1a0g, M_1 = pM_1a0g, TV1_1 = pTV1_0a0g, A_1 = 0, MC[,c("C1","C2")]), type = "response")
    pTV1_2a0g    <- noisy_est("binary", pTV1_2a0g)
    pM_2a0g      <- G$M_2
    
    #Step 3b - Simulate the outcome given each of the person-specific i histories based on the estimated model coefficients of Step 1b and histories generated in Step 3a
    pY_2a1g      <- predict(mY_2, newdata = data.frame(M_2 = pM_2a1g, TV1_2 = pTV1_2a1g, A_2 = 1, Y_1 = pY_1a1g, M_1 = pM_1a1g, TV1_1 = pTV1_1a1g, A_1 = 1, MC[,c("C1","C2")]), type = "response")
    pY_2a0g      <- predict(mY_2, newdata = data.frame(M_2 = pM_2a0g, TV1_2 = pTV1_2a0g, A_2 = 0, Y_1 = pY_1a0g, M_1 = pM_1a0g, TV1_1 = pTV1_1a0g, A_1 = 0, MC[,c("C1","C2")]), type = "response")
    
    MC$pY_2a1g <- pY_2a1g
    MC$pY_2a0g <- pY_2a0g
    
    #Step 3c - Estimate and store the mean outcome under each joint exposure-mediator intervention within each Step 3d K repetition
    if (a == 1) {
      outcomes$Ya1g1 <- mean(MC$pY_2a1g)
      outcomes$Ya0g1 <- mean(MC$pY_2a0g)
    }
    if (a == 0) {
      outcomes$Ya1g0 <- mean(MC$pY_2a1g)
      outcomes$Ya0g0 <- mean(MC$pY_2a0g)
    }
  }
  
  #Step 3c (continued) - Calculate each causal effect of interest including the TE, rPNDE, rTNDE, rPNIE, and rTNIE on the absolute and multiplicative scales
  outcomes$TE     <- outcomes$Ya1g1 - outcomes$Ya0g0
  outcomes$pNDE   <- outcomes$Ya1g0 - outcomes$Ya0g0
  outcomes$tNDE   <- outcomes$Ya1g1 - outcomes$Ya0g1
  outcomes$pNIE   <- outcomes$Ya0g1 - outcomes$Ya0g0
  outcomes$tNIE   <- outcomes$Ya1g1 - outcomes$Ya1g0
  outcomes$TE_m   <- outcomes$Ya1g1 / outcomes$Ya0g0
  outcomes$pNDE_m <- outcomes$Ya1g0 / outcomes$Ya0g0
  outcomes$tNDE_m <- outcomes$Ya1g1 / outcomes$Ya0g1
  outcomes$pNIE_m <- outcomes$Ya0g1 / outcomes$Ya0g0
  outcomes$tNIE_m <- outcomes$Ya1g1 / outcomes$Ya1g0
  
  #Step 3e - Estimate each joint intervention mean outcome as the mean of each Step 3d K repetition, and store it within the 'results' dataset for estimating the 95% CI
  #Adding above summary measures to 'results' dataset for 95% CI estimation  
  results <- data.frame(rep    = boot,
                        Y_obs  = mean(outcomes$Y_obs),   Y_nat  = mean(outcomes$Y_nat),
                        Ya0g0  = mean(outcomes$Ya0g0),   Ya0g1  = mean(outcomes$Ya0g1),
                        Ya1g0  = mean(outcomes$Ya1g0),   Ya1g1  = mean(outcomes$Ya1g1),
                        TE     = mean(outcomes$TE),
                        pNDE   = mean(outcomes$pNDE),    tNDE   = mean(outcomes$tNDE),
                        pNIE   = mean(outcomes$pNIE),    tNIE   = mean(outcomes$tNIE),
                        TE_m   = mean(outcomes$TE_m),
                        pNDE_m = mean(outcomes$pNDE_m),  tNDE_m = mean(outcomes$tNDE_m),
                        pNIE_m = mean(outcomes$pNIE_m),  tNIE_m = mean(outcomes$tNIE_m))
  return(results)
}

#Set seed for replicability
set.seed(2021)

#### Performing mediational g-formula in parallel ####
start <- Sys.time()
print(paste0("Mediational g-formula started at: ", start))

#### Load template data ####
dataset <- read.csv("template_data.csv")

#Setting number of cores to use as N-1 of available CPU cores
clusters <- parallel::makeCluster((detectCores() - 1))
doParallel::registerDoParallel(clusters)

#Performing mediational g-formula and saving output from each bootstrap in 'boot_output'
# - The number of bootstraps is set by boot (here as 200, +0 for the original dataset)
assign("boot_output", foreach(boot = 0:200, .combine = rbind) %dopar% medgf(dataset))
parallel::stopCluster(clusters)
print(paste0("Mediational g-formula finished at: ", Sys.time()))
print(paste0("Time elapsed: ", round((Sys.time() - start),3), " secs/mins"))
rm(start, clusters)

#Saving bootstrap findings
write.csv(boot_output, "boot_results_date.csv", row.names = FALSE, na = "")

#### Combining bootstrap findings to give summary results ####
boot_output_est <- boot_output[1,] #Take point-estimates from original data
boot_output_ci <- boot_output[2:201,] #Base SE, CI on bootstrap resampled data
medgf_output <- data.frame(metric = names(boot_output)[2:17], est = rep(NA, 16), SE = rep(NA, 16), LCI = rep(NA, 16), UCI = rep(NA, 16))

for (metric in medgf_output$metric) {
  medgf_output$est[medgf_output$metric == metric]    <- boot_output_est[[metric]] 
  medgf_output$SE[medgf_output$metric == metric]     <- sd(boot_output_ci[[metric]])
  medgf_output$LCI[medgf_output$metric == metric]    <- mean(boot_output_ci[[metric]]) - 1.96*sd(boot_output_ci[[metric]])
  medgf_output$UCI[medgf_output$metric == metric]    <- mean(boot_output_ci[[metric]]) + 1.96*sd(boot_output_ci[[metric]])
}
rm(metric)

#Summary results preview
cbind(metric = medgf_output[,1], round(medgf_output[, 2:5], 3))
