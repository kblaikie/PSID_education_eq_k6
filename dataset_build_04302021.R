#### Creating working mediational g-formula dataset based on all i ####
#Author - Kieran Blaikie
#Date - 30/4/2021

#Changes - Same as dataset_build_04262021 but creates srh_vgood_exc not srh_poor_fair

#Loading libraries
library(foreign)
library(haven)
library(tidyverse)
library(zoo)

#Loading EQ dataset, Jerzy dataset
directory <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"
eq <- read.csv(paste0(directory, "eq_imputed.csv"))
demog <- read.csv(paste0(directory, "R_data_5_18_20.csv"))

#Subsetting to variables needed
eq <- eq[,c("unique_id", "year", "age_unique", "eq_score", "eq_score_imp_pct_age")]
demog <- demog[,c("unique_id", "year", "age_unique", "gender", "race", "ethnicity", 
                  "race_hisp", "nativity", "education", "srh", "disabl_limits_work", 
                  "employment_status", "Region", "ind1990_broad", 
                  "occ1990dd_broad", "self_employed", "parents_poor", 
                  "marital_status_revised", "k6")] 

#Merging datasets, subsetting to observations:
#    1) With i between 30-60 years old
#    2) In or after 2001 (first time assessing K-6)
#    3) For i not permanently self-employed or unemployed for full follow-up
medgf <- merge(eq, demog, by=c("unique_id", "year", "age_unique"), all.x = T)
rm(eq, demog)

#Data manipulation
#    Creates
#       - 'is_eq_imputed' (=1 is.na(eq_score), =0 otherwise)
#       - 'eq_score_trim' (eq_score trimmed to >=1% eCDF, <=99% eCDF)
#       - 'eq_score_imp_pct_age_trim' (eq_score_imp_pct_age trimmed as above)
#       - 'female' (=1 gender-female, =0 gender-male)
#       - 'poc' (=1 race-Black,Other, =0 race-White)
#       - 'hispanic' (=1 ethnicity-Hispanic, =0 ethnicity-Non-Hispanic)
#       - 'poc_hisp' (=1 poc-1|hispanic-1, =0 poc-0 & hispanic-0)
#       - 'foreign' (=1 nativity-Grew up foreign, =0 nativity-Grew up US)
#       - 'hs_or_less' (=1 education-<HS|HS, =0 education-Some college|College+)
#       - 'srh_vgood_exc' (=1 srh-Very good|Excellent, =0 srh-Good|Fair|Poor)
#       - 'disabled_work' (=1 disabl_limits_work-Yes, =0 disabl_limits_work-No)
#       - 'unemployed' (=1 employment_status-Unemployed, =0 otherwise not NA)
#       - 'low_eq_occ' (=1 eq_score~occ1990dd_broad < mean eq, =0 otherwise not NA)
#              - medgf %>% group_by(occ1990dd_broad) %>% summarise(low_eq = mean(eq_score_imp_pct_age_trim) < -0.005278194)
#       - 'any_self_employed' (=1 self_employed-Self only|Both someone else and 
#                              self, =0 otherwise)
#       - 'parents_wealth' (=parents_poor variable)
#       - 'parents_poor' (=1 parents_wealth-Poor, =0 otherwise)
#       - 'married' (=1 marital_status_revised-Married/cohabitating, 
#                    =0 otherwise not NA)
medgf %>% mutate(is_eq_imputed = ifelse(is.na(medgf$eq_score), 1, 0),
                 eq_score_trim = ifelse(eq_score < quantile(eq_score, c(0.01), na.rm=T), quantile(eq_score, c(0.01), na.rm=T),
                                 ifelse(eq_score > quantile(eq_score, c(0.99), na.rm=T), quantile(eq_score, c(0.99), na.rm=T),
                                 ifelse(is.na(eq_score), NA, eq_score))),
                 eq_score_imp_pct_age_trim = ifelse(eq_score_imp_pct_age < quantile(eq_score_imp_pct_age, c(0.01), na.rm=T), quantile(eq_score_imp_pct_age, c(0.01), na.rm=T),
                                             ifelse(eq_score_imp_pct_age > quantile(eq_score_imp_pct_age, c(0.99), na.rm=T), quantile(eq_score_imp_pct_age, c(0.99), na.rm=T),
                                             ifelse(is.na(eq_score_imp_pct_age), NA, eq_score_imp_pct_age))),
                 female = ifelse(gender == "Female", 1, 0),
                 poc = ifelse(race == "White", 0, 1),
                 hispanic = ifelse(ethnicity == "Hispanic", 1, 0),
                 poc_hisp = ifelse(poc == 1 | hispanic == 1, 1, 0),
                 foreign = ifelse(nativity == "Grew up foreign", 1, 0),
                 hs_or_less = ifelse(education == "HS" | education == "<HS", 1, 0),
                 srh_vgood_exc = ifelse(srh == "Very good" | srh == "Excellent", 1, 0),
                 disabled_work = ifelse(disabl_limits_work=="Yes", 1, 0),
                 unemployed = ifelse(employment_status == "Unemployed", 1, 0),
                 low_eq_occ = ifelse((occ1990dd_broad == "Farming, forestry, and fishing" | 
                                      occ1990dd_broad == "Inapp." | 
                                      occ1990dd_broad == "Operators, fabricators, and laborers" | 
                                      occ1990dd_broad == "Services" | 
                                      occ1990dd_broad == "Technical, sales, and admin support"), 1, 0), 
                 any_self_employed = ifelse((self_employed == "Self only" | self_employed == "Both someone else and self"), 1, 0), 
                 parents_wealth = parents_poor,
                 parents_poor = ifelse(parents_wealth == "Poor", 1, 0),
                 married = ifelse(marital_status_revised == "Married/cohabiting", 1, 0)) -> medgf

#Restricting to only variables needed
medgf %>% select(unique_id, year, age_unique, female, race_hisp, poc, hispanic, 
                 poc_hisp, foreign, disabled_work, srh, srh_vgood_exc, parents_wealth, 
                 parents_poor, Region, occ1990dd_broad, low_eq_occ, 
                 any_self_employed, unemployed, married, k6, eq_score_trim, 
                 eq_score_imp_pct_age_trim, is_eq_imputed, hs_or_less) -> medgf

#Reformatting/reordering variables
medgf %>% mutate(race_hisp = as.factor(race_hisp),
                 race_hisp = relevel(race_hisp, "NH white"),
                 srh = as.factor(srh),
                 srh = relevel(srh, "Good"),
                 parents_wealth = as.factor(parents_wealth),
                 Region = as.factor(Region)) -> medgf

#Data manipulation
#    K-6 imputing 2005 K-6 based on mean of 2003-2007 if applicable
medgf %>% mutate(k6 = as.numeric(ifelse(k6 == "Not asked", NA, k6))) -> medgf
for (i in unique(medgf$unique_id)) {
  ind <- medgf[medgf$unique_id == i & medgf$year >= 2003 & medgf$year <=2007, c("year", "k6")]
  if (NROW(ind) == 3) {
    if ((!is.na(ind$k6[ind$year==2003])) & !is.na(ind$k6[ind$year==2007])) {
      medgf$k6[medgf$unique_id==i & medgf$year==2005] <- mean(medgf$k6[medgf$unique_id==i & (medgf$year==2003|medgf$year==2007)])
    }
  }
} 
rm(i, ind)

#Subsetting such that first and last K-6 & EQ score are not imputed
#Subsetting to observations per i >=first-k6/eq, <= last k6/eq
#NOTE - Takes ~10 minutes 
pb <- txtProgressBar(min=0, max=NROW(unique(medgf$unique_id)), style=3)
n <- 0
medgf_cut <- medgf[0, ]
for (i in unique(medgf$unique_id)) {
  n <- n+1
  ind <- medgf[medgf$unique_id==i, ]
  ind <- ind[ind$year >= min(ind$year[!is.na(ind$k6) & !is.na(ind$eq_score_trim)]), ]
  ind <- ind[ind$year <= max(ind$year[!is.na(ind$k6) & !is.na(ind$eq_score_trim)]), ]
  medgf_cut <- rbind(medgf_cut, ind)
  setTxtProgressBar(pb, n)
}  
rm(n, i, ind, pb)

#Replacing original with reduced dataset
medgf <- medgf_cut
rm(medgf_cut)

#Carry-forward then carry-backward for K-6 (not first or last K-6)
medgf %>% group_by(unique_id) %>%
  mutate(k6_impute = k6,
         k6_impute=na.locf(k6, na.rm=FALSE),
         k6_impute=na.locf(k6, na.rm=FALSE, fromLast=TRUE)) %>%
  ungroup() -> medgf

#Carry-forward then carry-backward for multiple vars
medgf %>% group_by(unique_id) %>%
  mutate(race_hisp=na.locf(race_hisp, na.rm=FALSE),
         race_hisp=na.locf(race_hisp, na.rm=FALSE, fromLast=TRUE),
         poc=na.locf(poc, na.rm=FALSE),
         poc=na.locf(poc, na.rm=FALSE, fromLast=TRUE),
         hispanic=na.locf(hispanic, na.rm=FALSE),
         hispanic=na.locf(hispanic, na.rm=FALSE, fromLast=TRUE),
         poc_hisp=na.locf(poc_hisp, na.rm=FALSE),
         poc_hisp=na.locf(poc_hisp, na.rm=FALSE, fromLast=TRUE),
         foreign=na.locf(foreign, na.rm=FALSE),
         foreign=na.locf(foreign, na.rm=FALSE, fromLast=TRUE),
         disabled_work=na.locf(disabled_work, na.rm=FALSE),
         disabled_work=na.locf(disabled_work, na.rm=FALSE, fromLast=TRUE),
         srh=na.locf(srh, na.rm=FALSE),
         srh=na.locf(srh, na.rm=FALSE, fromLast=TRUE),
         srh_vgood_exc=na.locf(srh_vgood_exc, na.rm=FALSE),
         srh_vgood_exc=na.locf(srh_vgood_exc, na.rm=FALSE, fromLast=TRUE),
         parents_wealth=na.locf(parents_wealth, na.rm=FALSE),
         parents_wealth=na.locf(parents_wealth, na.rm=FALSE, fromLast=TRUE),
         parents_poor=na.locf(parents_poor, na.rm=FALSE),
         parents_poor=na.locf(parents_poor, na.rm=FALSE, fromLast=TRUE),
         Region=na.locf(Region, na.rm=FALSE),
         Region=na.locf(Region, na.rm=FALSE, fromLast=TRUE),
         low_eq_occ=na.locf(low_eq_occ, na.rm=FALSE),
         low_eq_occ=na.locf(low_eq_occ, na.rm=FALSE, fromLast=TRUE),
         any_self_employed=na.locf(any_self_employed, na.rm=FALSE),
         any_self_employed=na.locf(any_self_employed, na.rm=FALSE, fromLast=TRUE),
         unemployed=na.locf(unemployed, na.rm=FALSE),
         unemployed=na.locf(unemployed, na.rm=FALSE, fromLast=TRUE),
         hs_or_less=na.locf(hs_or_less, na.rm=FALSE),
         hs_or_less=na.locf(hs_or_less, na.rm=FALSE, fromLast=TRUE)) %>%
  ungroup() -> medgf

#Removing those individuals still with NA for vars (removes 696)
medgf_cut <- na.omit(medgf[,c(1:20,23:26)])
na_names <- medgf$unique_id[(medgf$unique_id %in% medgf_cut$unique_id) == FALSE] 
medgf <- medgf[(medgf$unique_id %in% na_names) == FALSE, ]
rm(na_names, medgf_cut)

#Data manipulation
#    max_obs = maximum number of observations that could be used
medgf %>% group_by(unique_id) %>% mutate(max_obs = n()) %>% ungroup() -> medgf

#Data manipulation
#    Creates
#       - base_year (= first year for each person)
#       - base_disabled_work (= first year disability status)
#       - base_Region (= first year Region)
#       - employment_type (categ. employed/self-employed/unemployed)
#       - k6_l1 (1-lag K-6 based on k6_impute)
#       - eq_l1 (1-lag EQ based on eq_score_imp_pct_age_trim)
#       - hs_or_less_majority (=1 mean(hs_or_less)>=0.5, =0 otherwise)
#       - k6_final (= last year K-6 from k6_impute)
medgf %>% group_by(unique_id) %>% 
  mutate(base_year = min(year),
         last_year = max(year),
         base_age = age_unique[year == base_year],
         base_disabled_work = disabled_work[year == base_year],
         base_Region = Region[year == base_year],
         employment_type = ifelse(any_self_employed == 1, "self-employed",
                                  ifelse(unemployed == 1, "unemployed",
                                         "employed-notself")),
         employment_type = as.factor(employment_type),
         k6_l1 = lag(k6_impute),
         eq_l1 = lag(eq_score_imp_pct_age_trim),
         hs_or_less_majority = ifelse(mean(hs_or_less)>=0.5, 1, 0),
         k6_final = k6_impute[year == last_year],
         k6_impute_mod = ifelse(k6_impute >=5, 1, 0),
         k6_final_mod = ifelse(k6_final >=5, 1, 0),
         k6_l1_mod = ifelse(k6_l1 >=5, 1, 0),
         last_year = NULL) %>%
  ungroup() -> medgf

#Creating a time indicator variable 'seq_n'
#Note - Takes ~10 minutes
medgf$seq_n <- as.numeric(NA)
for (i in unique(medgf$unique_id)) {
  for (j in 1:NROW(medgf[medgf$unique_id==i,])) {
    medgf[medgf$unique_id==i,]$seq_n[j] <- j 
  }
}
medgf$seq_n <- medgf$seq_n-1
rm(i,j)

#We should now have all variables needed for analysis under different specifications
#Saving long-format copy of data
write.csv(medgf, paste0(directory, "r_dataset_long_04302021.csv"), row.names = FALSE, na = "")

#Reformatting data to wide-format for use in analyses
medgf %>% select(unique_id, seq_n, base_year, base_age, female, race_hisp, poc,
                 hispanic, poc_hisp, foreign, base_disabled_work, parents_wealth, 
                 parents_poor, base_Region, hs_or_less_majority, k6_final, 
                 k6_final_mod, age_unique, disabled_work, srh, srh_vgood_exc, 
                 Region, occ1990dd_broad, low_eq_occ, employment_type, married, 
                 hs_or_less, eq_score_trim, eq_score_imp_pct_age_trim, eq_l1, 
                 k6_impute, k6_impute_mod, k6_l1, k6_l1_mod, is_eq_imputed) -> medgf_wide
medgf_wide <- as.data.frame(medgf_wide) #reshape() needs data as data.frame, not tibble

medgf_wide <- reshape(medgf_wide, 
                      direction = "wide",
                      idvar = c("unique_id", "base_year", "base_age", 
                                "female", "race_hisp", "poc", "hispanic", 
                                "poc_hisp", "foreign", "base_disabled_work", 
                                "parents_wealth", "parents_poor", "base_Region", 
                                "hs_or_less_majority", "k6_final", "k6_final_mod"),
                      timevar = "seq_n",
                      sep = "_")

#Saving wide-format copy of data
write.csv(medgf_wide, paste0(directory, "r_dataset_wide_04302021.csv"), row.names = FALSE, na = "")
