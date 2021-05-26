#### Mediational g-formula template script #####
#Author - Kieran Blaikie
#Date - 26 May 2021

#Overview - This script provides a simplified template for the mediational g-formula if ran using wide-format data.
#           Running models in wide-format avoids having to deal with within-person dependence between observations
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

#### Creating a function to add noise to predictions for each i and return 0/1 values from predicted probabilities for binary variables ####
#   - In Steps 2 and 3, the predict() function is used extensively to simulate person-specific variable values for time-varying exposure, mediator, and outcome
#   - Via the predict() function, each i within covariate strata are assigned the same predicted mean (or response response probability)
#   - To reflect prediction uncertainty for these individuals, we add noise based on 1) the overall prediction uncertainty for continuous variables,
#     or 2) the strata-specific prediction uncertainty for binary variables 
#   - Note this function doesn't currently accommodate categorical variables or non-Normal continuous variables, though this could be in theory be accounted for
noisy_est <- function(type, prediction) {
  if (type == "binary") {
    noisy <- sapply(1:10000, function(x) {rbinom(1,1,prediction[x])})
  } else {
    noisy <- sapply(1:10000, function(x) {rnorm(1,prediction$fit[x], prediction$residual.scale)})
  }
  return(noisy)
}

#### Bootstraping loop for 95% CI ####
#Create a dataset to store each bootstrap output for computing 95% CI
#  - 200 repetitions seems sufficient for moment-based 95% CI, though more (e.g. 500) may be needed for quantile-based 95% CI
results <- data.frame(bootstrap  = rep(NA, 200),
                      Ya0g0  = rep(NA, 200),  Ya0g1  = rep(NA, 200),  Ya1g0  = rep(NA, 200), Ya1g1  = rep(NA, 200),  
                      pNDE   = rep(NA, 200),  tNDE   = rep(NA, 200),  pNDE_m = rep(NA, 200), tNDE_m = rep(NA, 200), 
                      pNIE   = rep(NA, 200),  tNIE   = rep(NA, 200),  pNIE_m = rep(NA, 200), tNIE_m = rep(NA, 200),
                      TE     = rep(NA, 200),  TE_m   = rep(NA, 200))

#Set seed for replicability
set.seed(2021)

#The bootstrapping process can run long periods of time. Creating a progress bar protects against insanity. 
pb <- txtProgressBar(min = 0, max = 200, initial = 0, style = 3)

for (boot in 1:200) {
  #Print current Bootstrap for guaging progress
  print(boot)
  
  #Store the current bootstrap number
  results$bootstrap[boot] <- boot
  
  #Create an equal size bootstrap resample of the original dataset to perform Steps 1-3 on
  index_boot <- sample(1:NROW(dataset), NROW(dataset), replace = T)
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
  
  #For Steps 2-3, first create a dataset to store our potential outcomes from each Step 3d K repetition 
  outcomes <- data.frame(K = rep(1:30), Ya0g0 = rep(NA, 30), Ya0g1 = rep(NA, 30), Ya1g0 = rep(NA, 30), Ya1g1 = rep(NA, 30))
  
  #Next resample with replacement from our original data to ensure we obtain an accurate simulation/well-behaved monte carlo estimator
  index <- sample(1:NROW(data_boot), size = 10000, replace = T)
  MC <- data_boot[index, ]
  
  #### Step 2 - Estimate the joint distribution of time-varying mediators under time-varying exposure interventions A, A* ####
  #           - Step 2 is solely to obtain a set of mediator values we would expect to see under each exposure intervention,
  #             that we can then permute and use in Step 3 for each joint A-M intervention
  
  #Step 2a - Set baseline covariates to the observed values for subject i, then recursively for each time t from 0-T-1:
  #          i  - Generate time t TV, A, M based on model coefficients generated in Step 1
  #          ii - Assign time t exposure under the intervention A=a
  #        - Note no need to predict A_N, given we are assigning A = a for each person  
  
  #Fixing A=1
  pTV1_0a1     <- predict(mTV1_0, newdata = data.frame(                                                                                     A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_0a1     <- noisy_est("binary", pTV1_0a1)
  pM_0a1       <- predict(mM_0, newdata = data.frame(                                                                     TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_0a1       <- noisy_est("continuous", pM_0a1)
  pY_0a1       <- predict(mY_0, newdata = data.frame(                                                       M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pY_0a1       <- noisy_est("binary", pY_0a1)
  
  pTV1_1a1     <- predict(mTV1_1, newdata = data.frame(                              A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_1a1     <- noisy_est("binary", pTV1_1a1)
  pM_1a1       <- predict(mM_1, newdata = data.frame(              TV1_1 = pTV1_1a1, A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_1a1       <- noisy_est("continuous", pM_1a1)
  pY_1a1       <- predict(mY_1, newdata = data.frame(M_1 = pM_1a1, TV1_1 = pTV1_1a1, A_1 = 1, Y_0 = pY_0a1, M_0 = pM_0a1, TV1_0 = pTV1_0a1, A_0 = 1, MC[,c("C1","C2")]), type = "response")
  pY_1a1       <- noisy_est("binary", pY_1a1)
  
  pTV1_2a1     <- predict(mTV1_2, newdata = data.frame(                              A_2 = 1, Y_1 = pY_1a1, M_1 = pM_1a1, TV1_1 = pTV1_0a1, A_1 = 1, MC[,c("C1","C2")]), type = "response")
  pTV1_2a1     <- noisy_est("binary", pTV1_2a1)
  pM_2a1       <- predict(mM_2, newdata = data.frame(              TV1_2 = pTV1_2a1, A_2 = 1, Y_1 = pY_1a1, M_1 = pM_1a1, TV1_1 = pTV1_1a1, A_1 = 1, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_2a1       <- noisy_est("continuous", pM_2a1)
  pY_2a1       <- predict(mY_2, newdata = data.frame(M_2 = pM_2a1, TV1_2 = pTV1_2a1, A_2 = 1, Y_1 = pY_1a1, M_1 = pM_1a1, TV1_1 = pTV1_1a1, A_1 = 1, MC[,c("C1","C2")]), type = "response")
  pY_2a1       <- noisy_est("binary", pY_2a1)
  
  #Step 2b - Below stores the predicted mediator for each t generated in Step 2a for use in Step 3. Random permutation is performed prior to Step 3 later.
  #   - Below is the predicted mediator for each i under A=1 at T=0 - T=4 (M_0 to M_4)
  perm_Mt_a1 <- data.frame(M_0 = pM_0a1,
                           M_1 = pM_1a1,
                           M_2 = pM_2a1)
  
  #Step 2c - Repeat Step 2a assigning A=0, with M at each time t predicted based on A=0
  
  #Fixing A=0
  pTV1_0a0     <- predict(mTV1_0, newdata = data.frame(                                                                                     A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_0a0     <- noisy_est("binary", pTV1_0a0)
  pM_0a0       <- predict(mM_0, newdata = data.frame(                                                                     TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_0a0       <- noisy_est("continuous", pM_0a0)
  pY_0a0       <- predict(mY_0, newdata = data.frame(                                                       M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pY_0a0       <- noisy_est("binary", pY_0a0)
  
  pTV1_1a0     <- predict(mTV1_1, newdata = data.frame(                              A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_1a0     <- noisy_est("binary", pTV1_1a0)
  pM_1a0       <- predict(mM_1, newdata = data.frame(              TV1_1 = pTV1_1a0, A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_1a0       <- noisy_est("continuous", pM_1a0)
  pY_1a0       <- predict(mY_1, newdata = data.frame(M_1 = pM_1a0, TV1_1 = pTV1_1a0, A_1 = 0, Y_0 = pY_0a0, M_0 = pM_0a0, TV1_0 = pTV1_0a0, A_0 = 0, MC[,c("C1","C2")]), type = "response")
  pY_1a0       <- noisy_est("binary", pY_1a0)
  
  pTV1_2a0     <- predict(mTV1_2, newdata = data.frame(                              A_2 = 0, Y_1 = pY_1a0, M_1 = pM_1a0, TV1_1 = pTV1_0a0, A_1 = 0, MC[,c("C1","C2")]), type = "response")
  pTV1_2a0     <- noisy_est("binary", pTV1_2a0)
  pM_2a0       <- predict(mM_2, newdata = data.frame(              TV1_2 = pTV1_2a0, A_2 = 0, Y_1 = pY_1a0, M_1 = pM_1a0, TV1_1 = pTV1_1a0, A_1 = 0, MC[,c("C1","C2")]), se.fit = T, interval = "prediction")
  pM_2a0       <- noisy_est("continuous", pM_2a0)
  pY_2a0       <- predict(mY_2, newdata = data.frame(M_2 = pM_2a0, TV1_2 = pTV1_2a0, A_2 = 0, Y_1 = pY_1a0, M_1 = pM_1a0, TV1_1 = pTV1_1a0, A_1 = 0, MC[,c("C1","C2")]), type = "response")
  pY_2a0       <- noisy_est("binary", pY_2a0)
  
  #Step 2d - Repeat Step 2b storing the predicted mediator for each t generated in Step 2c. Again, random permutation is performed prior to Step 3 later.
  perm_Mt_a0 <- data.frame(M_0 = pM_0a0,
                           M_1 = pM_1a0,
                           M_2 = pM_2a0)
  
  #### Step 3 - Estimate Q(a,a), Q(a,a*), Q(a*,a), Q(a*,a*) by repeating the below steps for each exposure and mediator assignment ####
  #           - In this script, this is accomplished by looping through A- 0:1, where for each A, the mediator 'G' is assigned as under that exposure intervention
  #           - The exposure is separately assigned within each loop in order to obtain each A,M, A,M*, A*,M, and A*,M* expected mean outcome Y at end of follow-up
  
  #Step 3d - Repeat Step 3a-c K (30) times using different random permutations of joint mediators M stored in Steps 2b and 2d
  #        - In this script, Step 3d is performed through creating a loop for Steps 3a-c, hence why it appears beforehand
  
  #Create another progress bar to keep track of each K repetition
  pk <- txtProgressBar(min = 0, max = 30, initial = 0, style = 3)
  
  for (k in 1:30) {
    
    #Step 3a - Recursively for each time t from 0 -to T-1:
    #          i   - Repeat Step 2ai replacing the time-varying exposure intervention a with the joint time-varying exposure and mediator intervention a-m 
    #          ii  - Assign the time t mediator as the ith component of the permuted vector for time t from 2b or 2d
    #          iii - Assign time t exposure under the intervention a or a*
    #        - The below simulates Y at T=3 K times, saving the mean Y for each K in the 'outcomes' dataset where A=0/1 and M=Ga
    #        - Where a=1 below, M=Ga based on M|A=1, where a=0 below, M=Ga* based on M|A=0
    
    for (a in 0:1) {
      #Step 2b/2d random permutation of joint mediator distribution K times, stored as intervention 'G' 
      if (a == 1) {
        G <- perm_Mt_a1[sample(NROW(perm_Mt_a1)), ]
      }
      if (a == 0) {
        G <- perm_Mt_a0[sample(NROW(perm_Mt_a0)), ]
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
      pY_2a1g      <- noisy_est("binary", pY_2a1g)
      pY_2a0g      <- predict(mY_2, newdata = data.frame(M_2 = pM_2a0g, TV1_2 = pTV1_2a0g, A_2 = 0, Y_1 = pY_1a0g, M_1 = pM_1a0g, TV1_1 = pTV1_1a0g, A_1 = 0, MC[,c("C1","C2")]), type = "response")
      pY_2a0g      <- noisy_est("binary", pY_2a0g)
      
      MC$pY_2a1g <- pY_2a1g
      MC$pY_2a0g <- pY_2a0g
      
      #Step 3c - Estimate and store the mean outcome under each joint exposure-mediator intervention within each Step 3d K repetition
      if (a == 1) {
        outcomes$Ya1g1[k] <- mean(MC$pY_2a1g)
        outcomes$Ya0g1[k] <- mean(MC$pY_2a0g)
      }
      if (a == 0) {
        outcomes$Ya1g0[k] <- mean(MC$pY_2a1g)
        outcomes$Ya0g0[k] <- mean(MC$pY_2a0g)
      }
    }
    
    #Update the K-repetition progress bar
    setTxtProgressBar(pk,k)
    
    #Clearing environment of all objects created in loop but not needed to improve speed
    rm(a, G, pM_0a0g, pM_0a1g, pM_1a0g, pM_1a1g, pM_2a0g, pM_2a1g, pTV1_0a0g, pTV1_0a1g, pTV1_1a0g, pTV1_1a1g, pTV1_2a0g, pTV1_2a1g, pY_0a0g, pY_0a1g, pY_1a0g, pY_1a1g, pY_2a0g, pY_2a1g)
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
  results$Ya0g0[boot]  <- mean(outcomes$Ya0g0)
  results$Ya0g1[boot]  <- mean(outcomes$Ya0g1)
  results$Ya1g0[boot]  <- mean(outcomes$Ya1g0)
  results$Ya1g1[boot]  <- mean(outcomes$Ya1g1)
  results$TE[boot]     <- mean(outcomes$TE)
  results$pNDE[boot]   <- mean(outcomes$pNDE)
  results$tNDE[boot]   <- mean(outcomes$tNDE)
  results$pNIE[boot]   <- mean(outcomes$pNIE)
  results$tNIE[boot]   <- mean(outcomes$tNIE)
  results$TE_m[boot]    <- mean(outcomes$TE_m)
  results$pNDE_m[boot]  <- mean(outcomes$pNDE_m)
  results$tNDE_m[boot]  <- mean(outcomes$tNDE_m)
  results$pNIE_m[boot]  <- mean(outcomes$pNIE_m)
  results$tNIE_m[boot]  <- mean(outcomes$tNIE_m)
  
  #Updating progress bar
  setTxtProgressBar(pb,boot)
  
  #Within each bootstrap, save the results dataset repeatedly in case of crashing
  write.csv(results, paste0(directory, "boot_results_date.csv"), row.names = FALSE, na = "")
  
  #Clearing environment of all objects created but not needed to improve speed
  rm(k, index, index_boot, pM_0a0, pM_0a1, pM_1a0, pM_1a1, pM_2a0, pM_2a1, pTV1_0a0, pTV1_0a1, pTV1_1a0, pTV1_1a1, pTV1_2a0, pTV1_2a1, pY_0a0, pY_0a1, pY_1a0, pY_1a1, pY_2a0, pY_2a1,
     data_boot, mA_0, mA_1, mA_2, mM_0, mM_1, mM_2, mTV1_0, mTV1_1, mTV1_2, mY_0, mY_1, mY_2, outcomes, MC, perm_Mt_a0, perm_Mt_a1)
} 
