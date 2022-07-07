#### Mediational g-formula script - marginal models - same timepoint Y #####
#Author - Kieran Blaikie
#Date - 11 Oct 2021
#Total runtime - ~20 hours

#Overview - This script creates 3-6 time-point marginal model with 95% CI
#           for multiple imputed datasets and complete-case samples
#         - EQ is treated as continuous, both PCA-based or linear
#         - K-6 is treated as binary, with the following assumed DAG
#           Base C > HS > rep[EQ > Unemp > SRH > MAR > K-6 ] > ... > T
#         - Baseline C 1-9: 1) year, 2) age, 3) sex, 4) race, 5) nativity, 
#                           6) disability, 7) parental wealth, 8) region, 9) occupation

#Changes - Compared to the medgformula_script_091221_parallel.R script:
#           1 - Only 40 multiply imputed datasets are used instead of 50
#           2 - Following model checks, polynomial terms are included in the 
#               making_medgformula_functions scripts
#           3 - The imputed datasets used differ as described in dataset_build_100121.R

#Loading libraries
library(doParallel)
library(tidyverse)

#Setting directory
directory <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"

#### Loading mediational g-formula functions for multiple imputation analysis ####
source("R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Scripts/making_medgformula_functions_101121.R")

#### Setting seed for replicability ####
set.seed(2021)

#### Performing mediational g-formula using PCA-based EQ with MI datasets ####
for (int in 3:6) {
  for (mi in 1:40) {
    start <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " Start time: ", start))
    
    #Loading multiply imputed dataset
    medgf_wide <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_", mi, ".csv"))
    
    #Subsetting to necessary variables for analysis
    medgf_wide %>% select(unique_id, base_ed, base_year, base_age, female, poc_hisp, 
                          nativity_not_US, base_disability, parents_poor, base_region, base_occupation, 
                          starts_with("unemp_or_nilf"), starts_with("srh_vgood_exc"), 
                          starts_with("married"), starts_with("eq"),
                          starts_with("k6_bin")) -> dataset
    rm(medgf_wide)
    
    #Simplifying naming in working dataset
    if (int == 3) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", 
                          "TV2_0", "TV2_1", "TV2_2", 
                          "TV3_0", "TV3_1", "TV3_2", 
                          "M_0", "M_1", "M_2",  
                          "Y_0", "Y_1", "Y_2")
    }
    if (int == 4) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3",
                          "M_0", "M_1", "M_2", "M_3", 
                          "Y_0", "Y_1", "Y_2", "Y_3")
    }
    if (int == 5) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3", "TV3_4",
                          "M_0", "M_1", "M_2", "M_3", "M_4", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4")
    }
    if (int == 6) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", "TV1_5", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4", "TV2_5",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3", "TV3_4", "TV3_5",
                          "M_0", "M_1", "M_2", "M_3", "M_4", "M_5", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5")
    }
    
    #Performing mediational g-formula in parallel and exporting output
    clusters <- parallel::makeCluster((detectCores()-1))
    doParallel::registerDoParallel(clusters)
    if (int == 3) {
      assign(paste0("int_3_mi_", mi), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_3(dataset))
    }
    if (int == 4) {
      assign(paste0("int_4_mi_", mi), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_4(dataset))
    }
    if (int == 5) {
      assign(paste0("int_5_mi_", mi), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_5(dataset))
    }
    if (int == 6) {
      assign(paste0("int_6_mi_", mi), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_6(dataset))
    }
    parallel::stopCluster(clusters)
    write.csv(eval(parse(text = paste0("int_", int, "_mi_", mi))), 
              paste0(directory, "MI_data/boot_results_101121_", int, "_mi_", mi, ".csv"), row.names = FALSE, na = "")
    end <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " completed. Time elapsed: ", (end - start)))
    print("")
    rm(dataset, start, end, clusters)
  }
}
rm(int, mi)

#### Extracting estimates from each PCA MI dataset where using the observed data ####
for (int in 3:6) {
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_", mi)))
    assign(paste0("int_", int, "_mi_observed_", mi), data[1, ])
    assign(paste0("int_", int, "_mi_", mi), data[2:201, ]) 
  }
  assign(paste0("int_", int, "_mi_mu_est"), eval(parse(text = paste0("int_", int, "_mi_1[0, ]"))))
  obs_mu <- eval(parse(text = paste0("int_", int, "_mi_mu_est")))
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_observed_", mi)))
    obs_mu <- rbind(obs_mu, data)
  }
  obs_mu %>% 
    select(Y_nat, YA_nat, Ya0_nat, Ya1_nat, YA0_nat, YA1_nat,
           Ya0g0, Ya0g1, Ya1g0, Ya1g1, TE, TE_m, pNDE, tNDE, 
           pNIE, tNIE, pNDE_m, tNDE_m, pNIE_m, tNIE_m) %>%
    mutate(Y_nat      = mean(Y_nat),          YA_nat = mean(YA_nat),
           Ya0_nat    = mean(Ya0_nat),        Ya1_nat    = mean(Ya1_nat),
           YA0_nat    = mean(YA0_nat),        YA1_nat    = mean(YA1_nat),
           Ya0g0      = mean(Ya0g0),          Ya0g1      = mean(Ya0g1),
           Ya1g0      = mean(Ya1g0),          Ya1g1      = mean(Ya1g1),
           TE         = mean(TE),             TE_m       = mean(TE_m),
           pNDE       = mean(pNDE),           tNDE       = mean(tNDE),
           pNIE       = mean(pNIE),           tNIE       = mean(tNIE),
           pNDE_m     = mean(pNDE_m),         tNDE_m     = mean(tNDE_m),
           pNIE_m     = mean(pNIE_m),         tNIE_m     = mean(tNIE_m)) %>% 
    distinct() -> obs_mu
  assign(paste0("int_", int, "_mi_mu_est"), obs_mu)
  rm(list = grep("_observed_", ls(), value = T))
  rm(obs_mu, data)
}

#### Combining findings to give summary findings for PCA MI analyses above ####
#This step provides:
# - Point estimates (calculated as the mean of MI-specific estimates)
# - Pooled SE estimates (calculated applying Rubin's Rules to MI-specific bootstrap SE)
# - 95% CI (calculated using pooled SE)
# - p-values (calculated using pooled SE, estimating df via Barnard & Rubin 1999 formula)
summary_all <- data.frame(Period = NA, Measure = NA, Point = NA, SE = NA, LCI = NA, UCI = NA, p_val = NA, Wald_T = NA, df_old = NA, df_obs = NA, df_adj = NA)
summary_all <- summary_all[0, ]
for (int in 3:6) {
  int_results <- data.frame(mi_set = NA,
                            eval(parse(text = paste0("int_", int, "_mi_", 1, "[1, c(2:21)]"))))
  int_results <- int_results[0, ]
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_", mi)))
    results <- data.frame(mi_set = mi,
                          Y_nat      = mean(data$Y_nat),          YA_nat     = mean(data$YA_nat),    
                          Ya0_nat    = mean(data$Ya0_nat),        Ya1_nat    = mean(data$Ya1_nat),
                          YA0_nat    = mean(data$YA0_nat),        YA1_nat    = mean(data$YA1_nat),
                          Ya0g0      = mean(data$Ya0g0),          Ya0g1      = mean(data$Ya0g1),
                          Ya1g0      = mean(data$Ya1g0),          Ya1g1      = mean(data$Ya1g1),
                          TE         = mean(data$TE),             TE_m       = mean(data$TE_m),
                          pNDE       = mean(data$pNDE),           tNDE       = mean(data$tNDE),
                          pNIE       = mean(data$pNIE),           tNIE       = mean(data$tNIE),
                          pNDE_m     = mean(data$pNDE_m),         tNDE_m     = mean(data$tNDE_m),
                          pNIE_m     = mean(data$pNIE_m),         tNIE_m     = mean(data$tNIE_m))
    
    results_se <- data.frame(mi_set = mi,
                             Y_nat      = sd(data$Y_nat),          YA_nat     = sd(data$YA_nat),    
                             Ya0_nat    = sd(data$Ya0_nat),        Ya1_nat    = sd(data$Ya1_nat),
                             YA0_nat    = sd(data$YA0_nat),        YA1_nat    = sd(data$YA1_nat),
                             Ya0g0      = sd(data$Ya0g0),          Ya0g1      = sd(data$Ya0g1),
                             Ya1g0      = sd(data$Ya1g0),          Ya1g1      = sd(data$Ya1g1),
                             TE         = sd(data$TE),             TE_m       = sd(data$TE_m),
                             pNDE       = sd(data$pNDE),           tNDE       = sd(data$tNDE),
                             pNIE       = sd(data$pNIE),           tNIE       = sd(data$tNIE),
                             pNDE_m     = sd(data$pNDE_m),         tNDE_m     = sd(data$tNDE_m),
                             pNIE_m     = sd(data$pNIE_m),         tNIE_m     = sd(data$tNIE_m))
    
    results <- rbind(results, results_se)
    int_results <- rbind(int_results, results)
    rm(results, results_se)
  }
  assign(paste0("int_", int, "_results"), int_results)
  int_point <- int_results[c(seq(from = 1, to = 79, by = 2)), c(2:21)]
  int_se    <- int_results[c(seq(from = 2, to = 80, by = 2)), c(2:21)]
  int_point_summary <- data.frame(Y_nat          = mean(int_point$Y_nat),          YA_nat     = mean(int_point$YA_nat),     
                                  Ya0_nat        = mean(int_point$Ya0_nat),        Ya1_nat    = mean(int_point$Ya1_nat),
                                  YA0_nat        = mean(int_point$YA0_nat),        YA1_nat    = mean(int_point$YA1_nat),
                                  Ya0g0          = mean(int_point$Ya0g0),          Ya0g1      = mean(int_point$Ya0g1),
                                  Ya1g0          = mean(int_point$Ya1g0),          Ya1g1      = mean(int_point$Ya1g1),
                                  TE             = mean(int_point$TE),             TE_m       = mean(int_point$TE_m),
                                  pNDE           = mean(int_point$pNDE),           tNDE       = mean(int_point$tNDE),
                                  pNIE           = mean(int_point$pNIE),           tNIE       = mean(int_point$tNIE),
                                  pNDE_m         = mean(int_point$pNDE_m),         tNDE_m     = mean(int_point$tNDE_m),
                                  pNIE_m         = mean(int_point$pNIE_m),         tNIE_m     = mean(int_point$tNIE_m))
  
  int_within_var_summary <- data.frame(Y_nat      = mean((int_se$Y_nat^2)),       YA_nat     = mean((int_se$YA_nat^2)),     
                                       Ya0_nat    = mean((int_se$Ya0_nat^2)),     Ya1_nat    = mean((int_se$Ya1_nat^2)),
                                       YA0_nat    = mean((int_se$YA0_nat^2)),     YA1_nat    = mean((int_se$YA1_nat^2)),
                                       Ya0g0      = mean((int_se$Ya0g0^2)),       Ya0g1  = mean((int_se$Ya0g1^2)),
                                       Ya1g0      = mean((int_se$Ya1g0^2)),       Ya1g1  = mean((int_se$Ya1g1^2)),
                                       TE         = mean((int_se$TE^2)),          TE_m   = mean((int_se$TE_m^2)),
                                       pNDE       = mean((int_se$pNDE^2)),        tNDE   = mean((int_se$tNDE^2)),
                                       pNIE       = mean((int_se$pNIE^2)),        tNIE   = mean((int_se$tNIE^2)),
                                       pNDE_m     = mean((int_se$pNDE_m^2)),      tNDE_m = mean((int_se$tNDE_m^2)),
                                       pNIE_m     = mean((int_se$pNIE_m^2)),      tNIE_m = mean((int_se$tNIE_m^2)))
  
  int_between_var_summary <- data.frame(Y_nat = sum((int_point$Y_nat - int_point_summary$Y_nat)^2)/39,
                                        YA_nat = sum((int_point$YA_nat - int_point_summary$YA_nat)^2)/39,
                                        Ya0_nat = sum((int_point$Ya0_nat - int_point_summary$Ya0_nat)^2)/39,
                                        Ya1_nat = sum((int_point$Ya1_nat - int_point_summary$Ya1_nat)^2)/39,
                                        YA0_nat = sum((int_point$YA0_nat - int_point_summary$YA0_nat)^2)/39,
                                        YA1_nat = sum((int_point$YA1_nat - int_point_summary$YA1_nat)^2)/39,
                                        Ya0g0 = sum((int_point$Ya0g0 - int_point_summary$Ya0g0)^2)/39,
                                        Ya0g1 = sum((int_point$Ya0g1 - int_point_summary$Ya0g1)^2)/39,
                                        Ya1g0 = sum((int_point$Ya1g0 - int_point_summary$Ya1g0)^2)/39,
                                        Ya1g1 = sum((int_point$Ya1g1 - int_point_summary$Ya1g1)^2)/39,
                                        TE   = sum((int_point$TE - int_point_summary$TE)^2)/39,
                                        TE_m = sum((int_point$TE_m - int_point_summary$TE_m)^2)/39,
                                        pNDE = sum((int_point$pNDE - int_point_summary$pNDE)^2)/39,
                                        tNDE = sum((int_point$tNDE - int_point_summary$tNDE)^2)/39,
                                        pNIE = sum((int_point$pNIE - int_point_summary$pNIE)^2)/39,
                                        tNIE = sum((int_point$tNIE - int_point_summary$tNIE)^2)/39,
                                        pNDE_m = sum((int_point$pNDE_m - int_point_summary$pNDE_m)^2)/39,
                                        tNDE_m = sum((int_point$tNDE_m - int_point_summary$tNDE_m)^2)/39,
                                        pNIE_m = sum((int_point$pNIE_m - int_point_summary$pNIE_m)^2)/39,
                                        tNIE_m = sum((int_point$tNIE_m - int_point_summary$tNIE_m)^2)/39)
  int_se_summary <- data.frame(Y_nat       = sqrt(sum(int_within_var_summary$Y_nat + int_between_var_summary$Y_nat + (int_between_var_summary$Y_nat/40))),
                               YA_nat      = sqrt(sum(int_within_var_summary$YA_nat + int_between_var_summary$YA_nat + (int_between_var_summary$YA_nat/40))),
                               Ya0_nat     = sqrt(sum(int_within_var_summary$Ya0_nat + int_between_var_summary$Ya0_nat + (int_between_var_summary$Ya0_nat/40))),
                               Ya1_nat     = sqrt(sum(int_within_var_summary$Ya1_nat + int_between_var_summary$Ya1_nat + (int_between_var_summary$Ya1_nat/40))),
                               YA0_nat     = sqrt(sum(int_within_var_summary$YA0_nat + int_between_var_summary$YA0_nat + (int_between_var_summary$YA0_nat/40))),
                               YA1_nat     = sqrt(sum(int_within_var_summary$YA1_nat + int_between_var_summary$YA1_nat + (int_between_var_summary$YA1_nat/40))),
                               Ya0g0  = sqrt(sum(int_within_var_summary$Ya0g0 + int_between_var_summary$Ya0g0 + (int_between_var_summary$Ya0g0/40))),
                               Ya0g1  = sqrt(sum(int_within_var_summary$Ya0g1 + int_between_var_summary$Ya0g1 + (int_between_var_summary$Ya0g1/40))),
                               Ya1g0  = sqrt(sum(int_within_var_summary$Ya1g0 + int_between_var_summary$Ya1g0 + (int_between_var_summary$Ya1g0/40))),
                               Ya1g1  = sqrt(sum(int_within_var_summary$Ya1g1 + int_between_var_summary$Ya1g1 + (int_between_var_summary$Ya1g1/40))),
                               TE     = sqrt(sum(int_within_var_summary$TE + int_between_var_summary$TE + (int_between_var_summary$TE/40))),
                               TE_m   = sqrt(sum(int_within_var_summary$TE_m + int_between_var_summary$TE_m + (int_between_var_summary$TE_m/40))),
                               pNDE   = sqrt(sum(int_within_var_summary$pNDE + int_between_var_summary$pNDE + (int_between_var_summary$pNDE/40))),
                               tNDE   = sqrt(sum(int_within_var_summary$tNDE + int_between_var_summary$tNDE + (int_between_var_summary$tNDE/40))),
                               pNIE   = sqrt(sum(int_within_var_summary$pNIE + int_between_var_summary$pNIE + (int_between_var_summary$pNIE/40))),
                               tNIE   = sqrt(sum(int_within_var_summary$tNIE + int_between_var_summary$tNIE + (int_between_var_summary$tNIE/40))),
                               pNDE_m = sqrt(sum(int_within_var_summary$pNDE_m + int_between_var_summary$pNDE_m + (int_between_var_summary$pNDE_m/40))),
                               tNDE_m = sqrt(sum(int_within_var_summary$tNDE_m + int_between_var_summary$tNDE_m + (int_between_var_summary$tNDE_m/40))),
                               pNIE_m = sqrt(sum(int_within_var_summary$pNIE_m + int_between_var_summary$pNIE_m + (int_between_var_summary$pNIE_m/40))),
                               tNIE_m = sqrt(sum(int_within_var_summary$tNIE_m + int_between_var_summary$tNIE_m + (int_between_var_summary$tNIE_m/40))))
  int_lambda_summary <- data.frame(Y_nat = NA,
                                   YA_nat = NA,
                                   Ya0_nat = NA,
                                   Ya1_nat = NA,
                                   YA0_nat = NA,
                                   YA1_nat = NA,
                                   Ya0g0 = NA,
                                   Ya0g1 = NA,
                                   Ya1g0 = NA,
                                   Ya1g1 = NA,
                                   TE   = (int_between_var_summary$TE + (int_between_var_summary$TE/40))/int_se_summary$TE,
                                   TE_m = (int_between_var_summary$TE_m + (int_between_var_summary$TE_m/40))/int_se_summary$TE_m,
                                   pNDE = (int_between_var_summary$pNDE + (int_between_var_summary$pNDE/40))/int_se_summary$pNDE,
                                   tNDE = (int_between_var_summary$tNDE + (int_between_var_summary$tNDE/40))/int_se_summary$tNDE,
                                   pNIE = (int_between_var_summary$pNIE + (int_between_var_summary$pNIE/40))/int_se_summary$pNIE,
                                   tNIE = (int_between_var_summary$tNIE + (int_between_var_summary$tNIE/40))/int_se_summary$tNIE,
                                   pNDE_m = (int_between_var_summary$pNDE_m + (int_between_var_summary$pNDE_m/40))/int_se_summary$pNDE_m,
                                   tNDE_m = (int_between_var_summary$tNDE_m + (int_between_var_summary$tNDE_m/40))/int_se_summary$tNDE_m,
                                   pNIE_m = (int_between_var_summary$pNIE_m + (int_between_var_summary$pNIE_m/40))/int_se_summary$pNIE_m,
                                   tNIE_m = (int_between_var_summary$tNIE_m + (int_between_var_summary$tNIE_m/40))/int_se_summary$tNIE_m)
  int_lci_summary <- data.frame(Y_nat       = int_point_summary$Y_nat - 1.96*int_se_summary$Y_nat,
                                YA_nat      = int_point_summary$YA_nat - 1.96*int_se_summary$YA_nat,
                                Ya0_nat     = int_point_summary$Ya0_nat - 1.96*int_se_summary$Ya0_nat,
                                Ya1_nat     = int_point_summary$Ya1_nat - 1.96*int_se_summary$Ya1_nat,
                                YA0_nat     = int_point_summary$YA0_nat - 1.96*int_se_summary$YA0_nat,
                                YA1_nat     = int_point_summary$YA1_nat - 1.96*int_se_summary$YA1_nat,
                                Ya0g0  = int_point_summary$Ya0g0 - 1.96*int_se_summary$Ya0g0,
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
  int_uci_summary <- data.frame(Y_nat       = int_point_summary$Y_nat + 1.96*int_se_summary$Y_nat,
                                YA_nat      = int_point_summary$YA_nat + 1.96*int_se_summary$YA_nat,
                                Ya0_nat     = int_point_summary$Ya0_nat + 1.96*int_se_summary$Ya0_nat,
                                Ya1_nat     = int_point_summary$Ya1_nat + 1.96*int_se_summary$Ya1_nat,
                                YA0_nat     = int_point_summary$YA0_nat + 1.96*int_se_summary$YA0_nat,
                                YA1_nat     = int_point_summary$YA1_nat + 1.96*int_se_summary$YA1_nat,
                                Ya0g0  = int_point_summary$Ya0g0 + 1.96*int_se_summary$Ya0g0,
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
  int_summary <- rbind(int_point_summary, int_se_summary, int_lci_summary, int_uci_summary, int_lambda_summary)
  int_summary <- cbind(data.frame(Measure = c("Point", "SE", "LCI", "UCI", "lambda")), int_summary)
  int_summary <- as.data.frame(t(int_summary))
  names(int_summary) <- int_summary[1,]
  int_summary <- int_summary[2:21, ]
  int_summary$Measure <- row.names(int_summary)
  int_summary$Period <- int*2
  int_summary$T_val <- ifelse((grepl("Y", int_summary$Measure)==T), NA, 
                              ifelse(((grepl("_m", int_summary$Measure)==T) & as.numeric(int_summary$Point) >= 1), ((as.numeric(int_summary$Point)-1)/as.numeric(int_summary$SE)), 
                                     ifelse(((grepl("_m", int_summary$Measure)==T) & as.numeric(int_summary$Point) <1), ((1/(as.numeric(int_summary$Point)))-1/as.numeric(int_summary$SE)), 
                                            ((as.numeric(int_summary$Point)/as.numeric(int_summary$SE))))))
  int_summary$df_old <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               (39/as.numeric(int_summary$lambda)^2))
  dataset <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_1.csv")) #Needed for estimating observed df via Barnard and Rubin (1999) method
  sample <- NROW(dataset)
  rm(dataset)
  int_summary$df_obs <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               ((((sample - 10)+1)/((sample - 10)+3))*(sample - 10)*(1 - as.numeric(int_summary$lambda)))) #I.e. [(n-k)+1/(n-k)+3]*(n-k)*(1-lambda)
  int_summary$df_adj <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               ((int_summary$df_old*int_summary$df_obs)/(int_summary$df_old+int_summary$df_obs)))
  int_summary$p_val <- (2*(1 - pt(q = abs(int_summary$T_val), df = int_summary$df_adj)))
  int_summary <- int_summary[,c("Period", "Measure", "Point", "SE", "LCI", "UCI", "p_val", "T_val", "df_old", "df_obs", "df_adj")]
  assign(paste0("int_", int, "_summary"), int_summary)
  write.csv(eval(parse(text = paste0("int_", int, "_summary"))), paste0(directory, "MI_data/pooled_summary_10112021_int_", int, ".csv"), row.names = FALSE, na = "")
  summary_all <- rbind(summary_all, int_summary)
  write.csv(summary_all, paste0(directory, "MI_data/pooled_summary_10112021.csv"), row.names = FALSE, na = "")
  rm(int_results, int_point, int_se, 
     int_point_summary, int_within_var_summary, int_between_var_summary, 
     int_se_summary, int_lci_summary, int_uci_summary, int_lambda_summary)
}
rm(int, mi)





#### Setting seed for replicability ####
set.seed(2021)
#### Performing mediational g-formula using PCA-based EQ with CC datasets ####
for (int in 3:6) {
  start <- Sys.time()
  print(paste0("Int ", int, " Start time: ", start))
  
  #Loading multiply imputed dataset
  medgf_wide <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_cc.csv"))
  
  #Subsetting to necessary variables for analysis
  medgf_wide %>% select(unique_id, base_ed, base_year, base_age, female, poc_hisp, 
                        nativity_not_US, base_disability, parents_poor, base_region, base_occupation, 
                        starts_with("srh_vgood_exc"), starts_with("married"), starts_with("eq"),
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
  
  #Performing mediational g-formula in parallel and exporting output
  clusters <- parallel::makeCluster((detectCores()-1))
  doParallel::registerDoParallel(clusters)
  if (int == 3) {
    assign(paste0("int_3_cc"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_3_cc(dataset))
  }
  if (int == 4) {
    assign(paste0("int_4_cc"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_4_cc(dataset))
  }
  if (int == 5) {
    assign(paste0("int_5_cc"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_5_cc(dataset))
  }
  if (int == 6) {
    assign(paste0("int_6_cc"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_6_cc(dataset))
  }
  parallel::stopCluster(clusters)
  write.csv(eval(parse(text = paste0("int_", int, "_cc"))), 
            paste0(directory, "MI_data/boot_results_101121_", int, "_cc.csv"), row.names = FALSE, na = "")
  end <- Sys.time()
  print(paste0("Int ", int, " completed. Time elapsed: ", (end - start)))
  print("")
  rm(dataset, start, end, clusters)
}
rm(int)



#### Extracting estimates from each PCA CC dataset where using the observed data ####
for (int in 3:6) {
  data <- eval(parse(text = paste0("int_", int, "_cc")))
  assign(paste0("int_", int, "_cc_obs"), data[1, c(2:21)])
  assign(paste0("int_", int, "_cc"), data[2:201, ]) 
}
rm(int)

#### Combining findings to give summary findings for PCA CC analyses above ####
#This step provides:
# - Point estimates (calculated as the mean of bootstrap-specific estimates)
# - SE estimates (calculated as SD of bootstrap estimates)
# - 95% CI (calculated using bootstrap SE)
# - p-value (calculated using dataset size, bootstrap SE for T-value, )
for (int in 3:6) {
  data <- eval(parse(text = paste0("int_", int, "_cc")))
  
  int_results <- data.frame(Y_nat      = mean(data$Y_nat),      YA_nat     = mean(data$YA_nat),
                            Ya0_nat    = mean(data$Ya0_nat),    Ya1_nat    = mean(data$Ya1_nat),
                            YA0_nat    = mean(data$YA0_nat),    YA1_nat    = mean(data$YA1_nat),
                            Ya0g0      = mean(data$Ya0g0),      Ya0g1      = mean(data$Ya0g1),
                            Ya1g0      = mean(data$Ya1g0),      Ya1g1      = mean(data$Ya1g1),
                            TE         = mean(data$TE),         TE_m       = mean(data$TE_m),
                            pNDE       = mean(data$pNDE),       tNDE       = mean(data$tNDE),
                            pNIE       = mean(data$pNIE),       tNIE       = mean(data$tNIE),
                            pNDE_m     = mean(data$pNDE_m),     tNDE_m     = mean(data$tNDE_m),
                            pNIE_m     = mean(data$pNIE_m),     tNIE_m     = mean(data$tNIE_m))
  
  int_stderror <- data.frame(Y_nat      = sd(data$Y_nat),      YA_nat     = sd(data$YA_nat),
                             Ya0_nat    = sd(data$Ya0_nat),    Ya1_nat    = sd(data$Ya1_nat),
                             YA0_nat    = sd(data$YA0_nat),    YA1_nat    = sd(data$YA1_nat),
                             Ya0g0      = sd(data$Ya0g0),      Ya0g1      = sd(data$Ya0g1),
                             Ya1g0      = sd(data$Ya1g0),      Ya1g1      = sd(data$Ya1g1),
                             TE         = sd(data$TE),         TE_m       = sd(data$TE_m),
                             pNDE       = sd(data$pNDE),       tNDE       = sd(data$tNDE),
                             pNIE       = sd(data$pNIE),       tNIE       = sd(data$tNIE),
                             pNDE_m     = sd(data$pNDE_m),     tNDE_m     = sd(data$tNDE_m),
                             pNIE_m     = sd(data$pNIE_m),     tNIE_m     = sd(data$tNIE_m))
  
  dataset <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_cc.csv"))
  sample <- NROW(dataset)
  
  results <- rbind(int_results, int_stderror)
  results <- as.data.frame(t(results))
  results$Measure <- rownames(results)
  rownames(results) <- NULL
  names(results) <- c("Point", "SE", "Measure")
  
  results$T_val <- ifelse((grepl("Y", results$Measure)==T), NA, 
                          ifelse(((grepl("_m", results$Measure)==T) & as.numeric(results$Point) >= 1), ((as.numeric(results$Point)-1)/(as.numeric(results$SE)/sqrt(sample))), 
                                 ifelse(((grepl("_m", results$Measure)==T) & as.numeric(results$Point) <1), ((1/(as.numeric(results$Point)))-1/(as.numeric(results$SE)/sqrt(sample))), 
                                        ((as.numeric(results$Point)/(as.numeric(results$SE)/sqrt(sample)))))))
  results$p_val <- (2*(1 - pt(q = abs(results$T_val), df = (sample-10))))
  results %>% mutate(LCI = Point - 1.96*SE, UCI = Point + 1.96*SE) %>%
    select(Measure, Point, SE, LCI, UCI, T_val, p_val) -> results
  
  assign(paste0("int_", int, "_cc_results"), results)
  write.csv(eval(parse(text = paste0("int_", int, "_cc_results"))), paste0(directory, "boot_summary_10112021_int_", int, "_cc.csv"), row.names = FALSE, na = "")
  rm(int_results, int_stderror, results)
}

results_all <- cbind(data.frame(Time = (rep(3:6, each = 20)*2)), (rbind(int_3_cc_results, int_4_cc_results, int_5_cc_results, int_6_cc_results)))
write.csv(results_all, paste0(directory, "boot_summary_10112021_cc.csv"), row.names = FALSE, na = "")



#### Setting seed for replicability ####
set.seed(2021)

#### Performing mediational g-formula using linear EQ with MI datasets ####
for (int in 3:6) {
  for (mi in 1:40) {
    start <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " Start time: ", start))
    
    #Loading multiply imputed dataset
    medgf_wide <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_", mi, ".csv"))
    
    #Subsetting to necessary variables for analysis
    medgf_wide %>% select(unique_id, base_ed, base_year, base_age, female, poc_hisp, 
                          nativity_not_US, base_disability, parents_poor, base_region, base_occupation, 
                          starts_with("unemp_or_nilf"), starts_with("srh_vgood_exc"), 
                          starts_with("married"), starts_with("linear_eq"),
                          starts_with("k6_bin")) -> dataset
    rm(medgf_wide)
    
    #Simplifying naming in working dataset
    if (int == 3) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", 
                          "TV2_0", "TV2_1", "TV2_2", 
                          "TV3_0", "TV3_1", "TV3_2", 
                          "M_0", "M_1", "M_2",  
                          "Y_0", "Y_1", "Y_2")
    }
    if (int == 4) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3",
                          "M_0", "M_1", "M_2", "M_3", 
                          "Y_0", "Y_1", "Y_2", "Y_3")
    }
    if (int == 5) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3", "TV3_4",
                          "M_0", "M_1", "M_2", "M_3", "M_4", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4")
    }
    if (int == 6) {
      names(dataset) <- c("ID", "A", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","C9",
                          "TV1_0", "TV1_1", "TV1_2", "TV1_3", "TV1_4", "TV1_5", 
                          "TV2_0", "TV2_1", "TV2_2", "TV2_3", "TV2_4", "TV2_5",
                          "TV3_0", "TV3_1", "TV3_2", "TV3_3", "TV3_4", "TV3_5",
                          "M_0", "M_1", "M_2", "M_3", "M_4", "M_5", 
                          "Y_0", "Y_1", "Y_2", "Y_3", "Y_4", "Y_5")
    }
    
    #Performing mediational g-formula in parallel and exporting output
    clusters <- parallel::makeCluster((detectCores()-1))
    doParallel::registerDoParallel(clusters)
    if (int == 3) {
      assign(paste0("int_3_mi_", mi, "_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_3(dataset))
    }
    if (int == 4) {
      assign(paste0("int_4_mi_", mi, "_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_4(dataset))
    }
    if (int == 5) {
      assign(paste0("int_5_mi_", mi, "_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_5(dataset))
    }
    if (int == 6) {
      assign(paste0("int_6_mi_", mi, "_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_6(dataset))
    }
    parallel::stopCluster(clusters)
    write.csv(eval(parse(text = paste0("int_", int, "_mi_", mi, "_lin"))), 
              paste0(directory, "MI_data/boot_results_101121_", int, "_mi_", mi, "_lin.csv"), row.names = FALSE, na = "")
    end <- Sys.time()
    print(paste0("Int ", int, " MI ", mi, " completed. Time elapsed: ", (end - start)))
    print("")
    rm(dataset, start, end, clusters)
  }
}
rm(int, mi)

#### Extracting estimates from each linear MI dataset where using the observed data ####
for (int in 3:6) {
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_", mi, "_lin")))
    assign(paste0("int_", int, "_mi_observed_", mi, "_lin"), data[1, ])
    assign(paste0("int_", int, "_mi_", mi, "_lin"), data[2:201, ]) 
  }
  assign(paste0("int_", int, "_mi_mu_est_lin"), eval(parse(text = paste0("int_", int, "_mi_1_lin[0, ]"))))
  obs_mu <- eval(parse(text = paste0("int_", int, "_mi_mu_est_lin")))
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_observed_", mi, "_lin")))
    obs_mu <- rbind(obs_mu, data)
  }
  obs_mu %>% 
    select(Y_nat, YA_nat, Ya0_nat, Ya1_nat, YA0_nat, YA1_nat,
           Ya0g0, Ya0g1, Ya1g0, Ya1g1, TE, TE_m, pNDE, tNDE, 
           pNIE, tNIE, pNDE_m, tNDE_m, pNIE_m, tNIE_m) %>%
    mutate(Y_nat      = mean(Y_nat),          YA_nat = mean(YA_nat),
           Ya0_nat    = mean(Ya0_nat),        Ya1_nat    = mean(Ya1_nat),
           YA0_nat    = mean(YA0_nat),        YA1_nat    = mean(YA1_nat),
           Ya0g0      = mean(Ya0g0),          Ya0g1      = mean(Ya0g1),
           Ya1g0      = mean(Ya1g0),          Ya1g1      = mean(Ya1g1),
           TE         = mean(TE),             TE_m       = mean(TE_m),
           pNDE       = mean(pNDE),           tNDE       = mean(tNDE),
           pNIE       = mean(pNIE),           tNIE       = mean(tNIE),
           pNDE_m     = mean(pNDE_m),         tNDE_m     = mean(tNDE_m),
           pNIE_m     = mean(pNIE_m),         tNIE_m     = mean(tNIE_m)) %>% 
    distinct() -> obs_mu
  assign(paste0("int_", int, "_mi_mu_est_lin"), obs_mu)
  rm(list = grep("_observed_", ls(), value = T))
  rm(obs_mu, data)
}

#### Combining findings to give summary findings for linear MI analyses above ####
#This step provides:
# - Point estimates (calculated as the mean of MI-specific estimates)
# - Pooled SE estimates (calculated applying Rubin's Rules to MI-specific bootstrap SE)
# - 95% CI (calculated using pooled SE)
# - p-values (calculated using pooled SE, estimating df via Barnard & Rubin 1999 formula)
summary_all <- data.frame(Period = NA, Measure = NA, Point = NA, SE = NA, LCI = NA, UCI = NA, p_val = NA, Wald_T = NA, df_old = NA, df_obs = NA, df_adj = NA)
summary_all <- summary_all[0, ]
for (int in 3:6) {
  int_results <- data.frame(mi_set = NA,
                            eval(parse(text = paste0("int_", int, "_mi_", 1, "_lin[1, c(2:21)]"))))
  int_results <- int_results[0, ]
  for (mi in 1:40) {
    data <- eval(parse(text = paste0("int_", int, "_mi_", mi, "_lin")))
    results <- data.frame(mi_set = mi,
                          Y_nat      = mean(data$Y_nat),          YA_nat     = mean(data$YA_nat),    
                          Ya0_nat    = mean(data$Ya0_nat),        Ya1_nat    = mean(data$Ya1_nat),
                          YA0_nat    = mean(data$YA0_nat),        YA1_nat    = mean(data$YA1_nat),
                          Ya0g0      = mean(data$Ya0g0),          Ya0g1      = mean(data$Ya0g1),
                          Ya1g0      = mean(data$Ya1g0),          Ya1g1      = mean(data$Ya1g1),
                          TE         = mean(data$TE),             TE_m       = mean(data$TE_m),
                          pNDE       = mean(data$pNDE),           tNDE       = mean(data$tNDE),
                          pNIE       = mean(data$pNIE),           tNIE       = mean(data$tNIE),
                          pNDE_m     = mean(data$pNDE_m),         tNDE_m     = mean(data$tNDE_m),
                          pNIE_m     = mean(data$pNIE_m),         tNIE_m     = mean(data$tNIE_m))
    
    results_se <- data.frame(mi_set = mi,
                             Y_nat      = sd(data$Y_nat),          YA_nat     = sd(data$YA_nat),    
                             Ya0_nat    = sd(data$Ya0_nat),        Ya1_nat    = sd(data$Ya1_nat),
                             YA0_nat    = sd(data$YA0_nat),        YA1_nat    = sd(data$YA1_nat),
                             Ya0g0      = sd(data$Ya0g0),          Ya0g1      = sd(data$Ya0g1),
                             Ya1g0      = sd(data$Ya1g0),          Ya1g1      = sd(data$Ya1g1),
                             TE         = sd(data$TE),             TE_m       = sd(data$TE_m),
                             pNDE       = sd(data$pNDE),           tNDE       = sd(data$tNDE),
                             pNIE       = sd(data$pNIE),           tNIE       = sd(data$tNIE),
                             pNDE_m     = sd(data$pNDE_m),         tNDE_m     = sd(data$tNDE_m),
                             pNIE_m     = sd(data$pNIE_m),         tNIE_m     = sd(data$tNIE_m))
    
    results <- rbind(results, results_se)
    int_results <- rbind(int_results, results)
    rm(results, results_se)
  }
  assign(paste0("int_", int, "_results_lin"), int_results)
  int_point <- int_results[c(seq(from = 1, to = 79, by = 2)), c(2:21)]
  int_se    <- int_results[c(seq(from = 2, to = 80, by = 2)), c(2:21)]
  int_point_summary <- data.frame(Y_nat          = mean(int_point$Y_nat),          YA_nat     = mean(int_point$YA_nat),     
                                  Ya0_nat        = mean(int_point$Ya0_nat),        Ya1_nat    = mean(int_point$Ya1_nat),
                                  YA0_nat        = mean(int_point$YA0_nat),        YA1_nat    = mean(int_point$YA1_nat),
                                  Ya0g0          = mean(int_point$Ya0g0),          Ya0g1      = mean(int_point$Ya0g1),
                                  Ya1g0          = mean(int_point$Ya1g0),          Ya1g1      = mean(int_point$Ya1g1),
                                  TE             = mean(int_point$TE),             TE_m       = mean(int_point$TE_m),
                                  pNDE           = mean(int_point$pNDE),           tNDE       = mean(int_point$tNDE),
                                  pNIE           = mean(int_point$pNIE),           tNIE       = mean(int_point$tNIE),
                                  pNDE_m         = mean(int_point$pNDE_m),         tNDE_m     = mean(int_point$tNDE_m),
                                  pNIE_m         = mean(int_point$pNIE_m),         tNIE_m     = mean(int_point$tNIE_m))
  
  int_within_var_summary <- data.frame(Y_nat      = mean((int_se$Y_nat^2)),       YA_nat     = mean((int_se$YA_nat^2)),     
                                       Ya0_nat    = mean((int_se$Ya0_nat^2)),     Ya1_nat    = mean((int_se$Ya1_nat^2)),
                                       YA0_nat    = mean((int_se$YA0_nat^2)),     YA1_nat    = mean((int_se$YA1_nat^2)),
                                       Ya0g0      = mean((int_se$Ya0g0^2)),       Ya0g1  = mean((int_se$Ya0g1^2)),
                                       Ya1g0      = mean((int_se$Ya1g0^2)),       Ya1g1  = mean((int_se$Ya1g1^2)),
                                       TE         = mean((int_se$TE^2)),          TE_m   = mean((int_se$TE_m^2)),
                                       pNDE       = mean((int_se$pNDE^2)),        tNDE   = mean((int_se$tNDE^2)),
                                       pNIE       = mean((int_se$pNIE^2)),        tNIE   = mean((int_se$tNIE^2)),
                                       pNDE_m     = mean((int_se$pNDE_m^2)),      tNDE_m = mean((int_se$tNDE_m^2)),
                                       pNIE_m     = mean((int_se$pNIE_m^2)),      tNIE_m = mean((int_se$tNIE_m^2)))
  
  int_between_var_summary <- data.frame(Y_nat = sum((int_point$Y_nat - int_point_summary$Y_nat)^2)/39,
                                        YA_nat = sum((int_point$YA_nat - int_point_summary$YA_nat)^2)/39,
                                        Ya0_nat = sum((int_point$Ya0_nat - int_point_summary$Ya0_nat)^2)/39,
                                        Ya1_nat = sum((int_point$Ya1_nat - int_point_summary$Ya1_nat)^2)/39,
                                        YA0_nat = sum((int_point$YA0_nat - int_point_summary$YA0_nat)^2)/39,
                                        YA1_nat = sum((int_point$YA1_nat - int_point_summary$YA1_nat)^2)/39,
                                        Ya0g0 = sum((int_point$Ya0g0 - int_point_summary$Ya0g0)^2)/39,
                                        Ya0g1 = sum((int_point$Ya0g1 - int_point_summary$Ya0g1)^2)/39,
                                        Ya1g0 = sum((int_point$Ya1g0 - int_point_summary$Ya1g0)^2)/39,
                                        Ya1g1 = sum((int_point$Ya1g1 - int_point_summary$Ya1g1)^2)/39,
                                        TE   = sum((int_point$TE - int_point_summary$TE)^2)/39,
                                        TE_m = sum((int_point$TE_m - int_point_summary$TE_m)^2)/39,
                                        pNDE = sum((int_point$pNDE - int_point_summary$pNDE)^2)/39,
                                        tNDE = sum((int_point$tNDE - int_point_summary$tNDE)^2)/39,
                                        pNIE = sum((int_point$pNIE - int_point_summary$pNIE)^2)/39,
                                        tNIE = sum((int_point$tNIE - int_point_summary$tNIE)^2)/39,
                                        pNDE_m = sum((int_point$pNDE_m - int_point_summary$pNDE_m)^2)/39,
                                        tNDE_m = sum((int_point$tNDE_m - int_point_summary$tNDE_m)^2)/39,
                                        pNIE_m = sum((int_point$pNIE_m - int_point_summary$pNIE_m)^2)/39,
                                        tNIE_m = sum((int_point$tNIE_m - int_point_summary$tNIE_m)^2)/39)
  int_se_summary <- data.frame(Y_nat       = sqrt(sum(int_within_var_summary$Y_nat + int_between_var_summary$Y_nat + (int_between_var_summary$Y_nat/40))),
                               YA_nat      = sqrt(sum(int_within_var_summary$YA_nat + int_between_var_summary$YA_nat + (int_between_var_summary$YA_nat/40))),
                               Ya0_nat     = sqrt(sum(int_within_var_summary$Ya0_nat + int_between_var_summary$Ya0_nat + (int_between_var_summary$Ya0_nat/40))),
                               Ya1_nat     = sqrt(sum(int_within_var_summary$Ya1_nat + int_between_var_summary$Ya1_nat + (int_between_var_summary$Ya1_nat/40))),
                               YA0_nat     = sqrt(sum(int_within_var_summary$YA0_nat + int_between_var_summary$YA0_nat + (int_between_var_summary$YA0_nat/40))),
                               YA1_nat     = sqrt(sum(int_within_var_summary$YA1_nat + int_between_var_summary$YA1_nat + (int_between_var_summary$YA1_nat/40))),
                               Ya0g0  = sqrt(sum(int_within_var_summary$Ya0g0 + int_between_var_summary$Ya0g0 + (int_between_var_summary$Ya0g0/40))),
                               Ya0g1  = sqrt(sum(int_within_var_summary$Ya0g1 + int_between_var_summary$Ya0g1 + (int_between_var_summary$Ya0g1/40))),
                               Ya1g0  = sqrt(sum(int_within_var_summary$Ya1g0 + int_between_var_summary$Ya1g0 + (int_between_var_summary$Ya1g0/40))),
                               Ya1g1  = sqrt(sum(int_within_var_summary$Ya1g1 + int_between_var_summary$Ya1g1 + (int_between_var_summary$Ya1g1/40))),
                               TE     = sqrt(sum(int_within_var_summary$TE + int_between_var_summary$TE + (int_between_var_summary$TE/40))),
                               TE_m   = sqrt(sum(int_within_var_summary$TE_m + int_between_var_summary$TE_m + (int_between_var_summary$TE_m/40))),
                               pNDE   = sqrt(sum(int_within_var_summary$pNDE + int_between_var_summary$pNDE + (int_between_var_summary$pNDE/40))),
                               tNDE   = sqrt(sum(int_within_var_summary$tNDE + int_between_var_summary$tNDE + (int_between_var_summary$tNDE/40))),
                               pNIE   = sqrt(sum(int_within_var_summary$pNIE + int_between_var_summary$pNIE + (int_between_var_summary$pNIE/40))),
                               tNIE   = sqrt(sum(int_within_var_summary$tNIE + int_between_var_summary$tNIE + (int_between_var_summary$tNIE/40))),
                               pNDE_m = sqrt(sum(int_within_var_summary$pNDE_m + int_between_var_summary$pNDE_m + (int_between_var_summary$pNDE_m/40))),
                               tNDE_m = sqrt(sum(int_within_var_summary$tNDE_m + int_between_var_summary$tNDE_m + (int_between_var_summary$tNDE_m/40))),
                               pNIE_m = sqrt(sum(int_within_var_summary$pNIE_m + int_between_var_summary$pNIE_m + (int_between_var_summary$pNIE_m/40))),
                               tNIE_m = sqrt(sum(int_within_var_summary$tNIE_m + int_between_var_summary$tNIE_m + (int_between_var_summary$tNIE_m/40))))
  int_lambda_summary <- data.frame(Y_nat = NA,
                                   YA_nat = NA,
                                   Ya0_nat = NA,
                                   Ya1_nat = NA,
                                   YA0_nat = NA,
                                   YA1_nat = NA,
                                   Ya0g0 = NA,
                                   Ya0g1 = NA,
                                   Ya1g0 = NA,
                                   Ya1g1 = NA,
                                   TE   = (int_between_var_summary$TE + (int_between_var_summary$TE/40))/int_se_summary$TE,
                                   TE_m = (int_between_var_summary$TE_m + (int_between_var_summary$TE_m/40))/int_se_summary$TE_m,
                                   pNDE = (int_between_var_summary$pNDE + (int_between_var_summary$pNDE/40))/int_se_summary$pNDE,
                                   tNDE = (int_between_var_summary$tNDE + (int_between_var_summary$tNDE/40))/int_se_summary$tNDE,
                                   pNIE = (int_between_var_summary$pNIE + (int_between_var_summary$pNIE/40))/int_se_summary$pNIE,
                                   tNIE = (int_between_var_summary$tNIE + (int_between_var_summary$tNIE/40))/int_se_summary$tNIE,
                                   pNDE_m = (int_between_var_summary$pNDE_m + (int_between_var_summary$pNDE_m/40))/int_se_summary$pNDE_m,
                                   tNDE_m = (int_between_var_summary$tNDE_m + (int_between_var_summary$tNDE_m/40))/int_se_summary$tNDE_m,
                                   pNIE_m = (int_between_var_summary$pNIE_m + (int_between_var_summary$pNIE_m/40))/int_se_summary$pNIE_m,
                                   tNIE_m = (int_between_var_summary$tNIE_m + (int_between_var_summary$tNIE_m/40))/int_se_summary$tNIE_m)
  int_lci_summary <- data.frame(Y_nat       = int_point_summary$Y_nat - 1.96*int_se_summary$Y_nat,
                                YA_nat      = int_point_summary$YA_nat - 1.96*int_se_summary$YA_nat,
                                Ya0_nat     = int_point_summary$Ya0_nat - 1.96*int_se_summary$Ya0_nat,
                                Ya1_nat     = int_point_summary$Ya1_nat - 1.96*int_se_summary$Ya1_nat,
                                YA0_nat     = int_point_summary$YA0_nat - 1.96*int_se_summary$YA0_nat,
                                YA1_nat     = int_point_summary$YA1_nat - 1.96*int_se_summary$YA1_nat,
                                Ya0g0  = int_point_summary$Ya0g0 - 1.96*int_se_summary$Ya0g0,
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
  int_uci_summary <- data.frame(Y_nat       = int_point_summary$Y_nat + 1.96*int_se_summary$Y_nat,
                                YA_nat      = int_point_summary$YA_nat + 1.96*int_se_summary$YA_nat,
                                Ya0_nat     = int_point_summary$Ya0_nat + 1.96*int_se_summary$Ya0_nat,
                                Ya1_nat     = int_point_summary$Ya1_nat + 1.96*int_se_summary$Ya1_nat,
                                YA0_nat     = int_point_summary$YA0_nat + 1.96*int_se_summary$YA0_nat,
                                YA1_nat     = int_point_summary$YA1_nat + 1.96*int_se_summary$YA1_nat,
                                Ya0g0  = int_point_summary$Ya0g0 + 1.96*int_se_summary$Ya0g0,
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
  int_summary <- rbind(int_point_summary, int_se_summary, int_lci_summary, int_uci_summary, int_lambda_summary)
  int_summary <- cbind(data.frame(Measure = c("Point", "SE", "LCI", "UCI", "lambda")), int_summary)
  int_summary <- as.data.frame(t(int_summary))
  names(int_summary) <- int_summary[1,]
  int_summary <- int_summary[2:21, ]
  int_summary$Measure <- row.names(int_summary)
  int_summary$Period <- int*2
  int_summary$T_val <- ifelse((grepl("Y", int_summary$Measure)==T), NA, 
                              ifelse(((grepl("_m", int_summary$Measure)==T) & as.numeric(int_summary$Point) >= 1), ((as.numeric(int_summary$Point)-1)/as.numeric(int_summary$SE)), 
                                     ifelse(((grepl("_m", int_summary$Measure)==T) & as.numeric(int_summary$Point) <1), ((1/(as.numeric(int_summary$Point)))-1/as.numeric(int_summary$SE)), 
                                            ((as.numeric(int_summary$Point)/as.numeric(int_summary$SE))))))
  int_summary$df_old <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               (39/as.numeric(int_summary$lambda)^2))
  dataset <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_1.csv")) #Needed for estimating observed df via Barnard and Rubin (1999) method
  sample <- NROW(dataset)
  rm(dataset)
  int_summary$df_obs <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               ((((sample - 10)+1)/((sample - 10)+3))*(sample - 10)*(1 - as.numeric(int_summary$lambda)))) #I.e. [(n-k)+1/(n-k)+3]*(n-k)*(1-lambda)
  int_summary$df_adj <- ifelse((grepl("Y", int_summary$Measure)==T), NA,
                               ((int_summary$df_old*int_summary$df_obs)/(int_summary$df_old+int_summary$df_obs)))
  int_summary$p_val <- (2*(1 - pt(q = abs(int_summary$T_val), df = int_summary$df_adj)))
  int_summary <- int_summary[,c("Period", "Measure", "Point", "SE", "LCI", "UCI", "p_val", "T_val", "df_old", "df_obs", "df_adj")]
  assign(paste0("int_", int, "_summary_lin"), int_summary)
  write.csv(eval(parse(text = paste0("int_", int, "_summary_lin"))), paste0(directory, "MI_data/pooled_summary_10112021_int_", int, "_lin.csv"), row.names = FALSE, na = "")
  summary_all <- rbind(summary_all, int_summary)
  write.csv(summary_all, paste0(directory, "MI_data/pooled_summary_10112021_lin.csv"), row.names = FALSE, na = "")
  rm(int_results, int_point, int_se, 
     int_point_summary, int_within_var_summary, int_between_var_summary, 
     int_se_summary, int_lci_summary, int_uci_summary, int_lambda_summary)
}
rm(int, mi)






#### Setting seed for replicability ####
set.seed(2021)

#### Performing mediational g-formula using linear EQ with CC datasets ####
for (int in 3:6) {
  start <- Sys.time()
  print(paste0("Int ", int, " Start time: ", start))
  
  #Loading multiply imputed dataset
  medgf_wide <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_cc.csv"))
  
  #Subsetting to necessary variables for analysis
  medgf_wide %>% select(unique_id, base_ed, base_year, base_age, female, poc_hisp, 
                        nativity_not_US, base_disability, parents_poor, base_region, base_occupation, 
                        starts_with("srh_vgood_exc"), starts_with("married"), starts_with("linear_eq"),
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
  
  #Performing mediational g-formula in parallel and exporting output
  clusters <- parallel::makeCluster((detectCores()-1))
  doParallel::registerDoParallel(clusters)
  if (int == 3) {
    assign(paste0("int_3_cc_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_3_cc(dataset))
  }
  if (int == 4) {
    assign(paste0("int_4_cc_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_4_cc(dataset))
  }
  if (int == 5) {
    assign(paste0("int_5_cc_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_5_cc(dataset))
  }
  if (int == 6) {
    assign(paste0("int_6_cc_lin"), foreach(boot = 0:200, .combine = rbind) %dopar% medgf_6_cc(dataset))
  }
  parallel::stopCluster(clusters)
  write.csv(eval(parse(text = paste0("int_", int, "_cc_lin"))), 
            paste0(directory, "MI_data/boot_results_101121_", int, "_cc_lin.csv"), row.names = FALSE, na = "")
  end <- Sys.time()
  print(paste0("Int ", int, " completed. Time elapsed: ", (end - start)))
  print("")
  rm(dataset, start, end, clusters)
}
rm(int)


#### Extracting estimates from each linear CC dataset where using the observed data ####
for (int in 3:6) {
  data <- eval(parse(text = paste0("int_", int, "_cc_lin")))
  assign(paste0("int_", int, "_cc_obs_lin"), data[1, ])
  assign(paste0("int_", int, "_cc_lin"), data[2:201, ]) 
}
rm(int)

#### Combining findings to give summary findings for linear CC analyses above ####
#This step provides:
# - Point estimates (calculated as the mean of bootstrap-specific estimates)
# - SE estimates (calculated as SD of bootstrap estimates)
# - 95% CI (calculated using bootstrap SE)
# - p-value (calculated using dataset size, bootstrap SE for T-value, )
for (int in 3:6) {
  data <- eval(parse(text = paste0("int_", int, "_cc_lin")))
  
  int_results <- data.frame(Y_nat      = mean(data$Y_nat),      YA_nat     = mean(data$YA_nat),
                            Ya0_nat    = mean(data$Ya0_nat),    Ya1_nat    = mean(data$Ya1_nat),
                            YA0_nat    = mean(data$YA0_nat),    YA1_nat    = mean(data$YA1_nat),
                            Ya0g0      = mean(data$Ya0g0),      Ya0g1      = mean(data$Ya0g1),
                            Ya1g0      = mean(data$Ya1g0),      Ya1g1      = mean(data$Ya1g1),
                            TE         = mean(data$TE),         TE_m       = mean(data$TE_m),
                            pNDE       = mean(data$pNDE),       tNDE       = mean(data$tNDE),
                            pNIE       = mean(data$pNIE),       tNIE       = mean(data$tNIE),
                            pNDE_m     = mean(data$pNDE_m),     tNDE_m     = mean(data$tNDE_m),
                            pNIE_m     = mean(data$pNIE_m),     tNIE_m     = mean(data$tNIE_m))
  
  int_stderror <- data.frame(Y_nat      = sd(data$Y_nat),      YA_nat     = sd(data$YA_nat),
                             Ya0_nat    = sd(data$Ya0_nat),    Ya1_nat    = sd(data$Ya1_nat),
                             YA0_nat    = sd(data$YA0_nat),    YA1_nat    = sd(data$YA1_nat),
                             Ya0g0      = sd(data$Ya0g0),      Ya0g1      = sd(data$Ya0g1),
                             Ya1g0      = sd(data$Ya1g0),      Ya1g1      = sd(data$Ya1g1),
                             TE         = sd(data$TE),         TE_m       = sd(data$TE_m),
                             pNDE       = sd(data$pNDE),       tNDE       = sd(data$tNDE),
                             pNIE       = sd(data$pNIE),       tNIE       = sd(data$tNIE),
                             pNDE_m     = sd(data$pNDE_m),     tNDE_m     = sd(data$tNDE_m),
                             pNIE_m     = sd(data$pNIE_m),     tNIE_m     = sd(data$tNIE_m))
  
  results <- rbind(int_results, int_stderror)
  results <- as.data.frame(t(results))
  results$Measure <- rownames(results)
  rownames(results) <- NULL
  names(results) <- c("Point", "SE", "Measure")
  
  dataset <- read.csv(paste0(directory, "MI_data/dataset_100121_wide_", int, "_cc.csv"))
  sample <- NROW(dataset)
  
  results$T_val <- ifelse((grepl("Y", results$Measure)==T), NA, 
                          ifelse(((grepl("_m", results$Measure)==T) & as.numeric(results$Point) >= 1), ((as.numeric(results$Point)-1)/(as.numeric(results$SE)/sqrt(sample))), 
                                 ifelse(((grepl("_m", results$Measure)==T) & as.numeric(results$Point) <1), ((1/(as.numeric(results$Point)))-1/(as.numeric(results$SE)/sqrt(sample))), 
                                        ((as.numeric(results$Point)/(as.numeric(results$SE)/sqrt(sample)))))))
  results$p_val <- (2*(1 - pt(q = abs(results$T_val), df = (sample-10))))
  results %>% mutate(LCI = Point - 1.96*SE, UCI = Point + 1.96*SE) %>%
    select(Measure, Point, SE, LCI, UCI,T_val,  p_val) -> results
  
  assign(paste0("int_", int, "_cc_results_lin"), results)
  write.csv(eval(parse(text = paste0("int_", int, "_cc_results_lin"))), paste0(directory, "boot_summary_10112021_int_", int, "_cc_lin.csv"), row.names = FALSE, na = "")
  rm(int_results, int_stderror, results)
}

results_all <- cbind(data.frame(Time = (rep(3:6, each = 20)*2)), (rbind(int_3_cc_results_lin, int_4_cc_results_lin, int_5_cc_results_lin, int_6_cc_results_lin)))
write.csv(results_all, paste0(directory, "boot_summary_10112021_cc_lin.csv"), row.names = FALSE, na = "")

