#### Creating datasets for mediational g-formula incorporating 2019 data ####
#Author - Kieran Blaikie
#Date started - 01 Oct 21
#Total Runtime - >5 hours

#Overview - This script does the following:
#            1 - Create an EQ variable in Stata as in 'EQ Script 2.17.21.R', 
#                as well as a linear EQ score more similar to the PES
#            2 - Create maximum and complete-case samples for 3-6 wave analyses
#            3 - Impute EQ in those unemployed, self-employed, NILF making same 
#                assumptions as in 'dataset_build_07222021.R'

#Changes - Compared to the script 'dataset_build_090921.R', the following changes are made:
#           1 - Annual hours worked and unemployment duration are left unassigned for those
#               unemployed or NILF before multiple imputation
#           2 - In multiple imputation, the dummy variable 'unemp_or_nilf' is separated
#               into two dummy variables for 'unemployed' and 'nilf' separately
#           3 - A linear EQ score variable is created
#        - Compared to the script 'dataset_build_091221.R', the following changes are made:
#           1 - This script excludes all those fully or partially self-employed
#           2 - This script uses the 091621 dataset, which leaves annual hours worked
#               for those unemployed or NILF as-reported, instead of coercing to 0
#        - Compared to the script 'dataset_build_092921.R', the following changes are made:
#           1 - In creating a linear EQ score, the material rewards income component
#               is now based on >= median age-standardised income instead of >= mean
#           2 - In creating a linear EQ score, the working time arrangements component
#               is now based on working >=32 hours per week instead of >=35 hours

#Loading needed libraries
library(tidyverse)   # For data-wrangling
library(mice)        # For multiple imputation
library(miceadds)    # For person-clustered imputation
library(haven)       # For saving .dta files
library(foreign)
library(readstata13) # For reading in EQ dataset made in Part 1
library(zoo)         #For carry-forward, carry-backward

#Setting directory
directory <- "R:/Project/precarityR01/PSID/data_management_and_descriptives/"

#Loading needed datasets
data <- read.csv(paste0(directory, "psid_dataset_091621.csv")) #1984-2019 Data
cpi <- read.csv(paste0(directory, "cpi.csv")) #Fed Reserve CPI

#Subsetting to eligible age sample (those 30-60 years old)
data <- data[data$age %in% c(30:60) & !is.na(data$age), ]

#Subsetting the eligible employment sample (unemployed, NILF, employed - not self)
data <- data[data$employment_status_cat == "NILF" | data$employment_status_cat == "Unemployed" | 
               (data$employment_status_cat == "Working or temp. layoff" & data$self_employed == "Someone else"), ]
data <- data[!is.na(data$employment_status_cat), ]

#### Part 1 ####

#### Reformatting employment variables to format needed for PCA ####
#PCA Variables include : Union (0/1), Employer-paid HI (0/1), Salaried (0/1),
#                        Overtime pay (0/1), Pension contribution (0/1), 
#                        Annual hours worked (N), Unemployment duration (N), 
#                        Z-Score Labor Income, Z-Score Employment tenure
#Adjusting 'total_labor_income' to reflect 2019 dollars before creating Z-Score
cpi %>% right_join(data, cpi, by="year") %>% 
  mutate(ratio = cpi/255.7, #Ratio of CPI in X year to CPI in 2019 (255.7)
         total_labor_income = total_labor_income/ratio) %>%
  select(-ratio) -> data
rm(cpi)

#Creating all needed PCA variables
#Notes - For employer health insurance, if a person is covered by employer health insurance
#        but it isn't known whether it is their employer (i.e. they are coded 'Yes' for 
#        'employer_health_insurance_1319_any' but NA for employer_health_insurance_1319')
#        these persons are coded as '1/Yes', assuming it's their own employer unless unemployed
#      - For unemployment duration, all those employed with unknown unemployment duration are 
#        recoded as having 0 months unemployed in the previous year
data %>% mutate(union = ifelse(current_job_union_member == "Yes", 1, 0),
                emp_hi = ifelse(employer_health_insurance_84_9919 == "Yes", 1, 0), 
                salaried = ifelse(wage_salary == "Salaried", 1, 0),
                overtime = case_when(wage_salary == "Salaried" ~ 1, wage_salary == "Other" ~ 0,
                                     (wage_salary == "Hourly" & (paid_overtime_waged == "Other" | 
                                                                   paid_overtime_waged == "Exact amount" | 
                                                                   paid_overtime_waged == "1.5x" |
                                                                   paid_overtime_waged == "2x")) ~ 1, 
                                     (wage_salary == "Hourly" & (paid_overtime_waged == "Other" | 
                                                                   paid_overtime_waged == "Comp time" | 
                                                                   paid_overtime_waged == "IAP" |
                                                                   paid_overtime_waged == "1x")) ~ 0,
                                     TRUE ~ NA_real_),
                pension = ifelse(current_job_pension == "Yes", 1, 0)) -> data

#Annual hours worked, unemployment duration variables don't need recoded
data %>% group_by(age) %>%
  summarise(mean_age_inc = mean(total_labor_income, na.rm = T),
            sd_age_inc = sd(total_labor_income, na.rm = T),
            mean_age_emp_tenure = mean(emp_tenure_months, na.rm = T),
            sd_age_emp_tenure = sd(emp_tenure_months,na.rm = T)) %>%
  right_join(data, by = "age") %>%
  mutate(z_total_labor_income = (total_labor_income - mean_age_inc)/sd_age_inc,
         z_emp_tenure_months = (emp_tenure_months - mean_age_emp_tenure)/sd_age_emp_tenure) %>%
  select(-mean_age_inc, -sd_age_inc, -mean_age_emp_tenure, -sd_age_emp_tenure) -> data

#PCA-based EQ was created in the dataset_build_092921.R script, so is loaded in here
#Loading in EQ dataset
#Note - Comparing EQ created here against that created in EQ Script 2.17.21.R, 
#       there is high correlation (Pearson R = 0.97)
eq_data <- read.dta13(paste0(directory, "eq_092921.dta"))
eq_data %>% select(unique_id, year, age, eq) -> eq_data

#Adding EQ into larger dataset
data <- merge(data, eq_data, by = c("unique_id", "year", "age"), all = T)
rm(eq_data)

#### Creating a linear EQ Score ####
#This is a 0-5 point scale, with the following scoring
#1 - Employment stability - Unemployed = 0, employed < 12 months = 0.5, > 12 months = 1
#2 - Material rewards - > Average Income = 1/3, Employer HI = 1/3, Employer pension = 1/3
#3 - Rights and protections - Overtime pay = 0.5, Salaried employment = 0.5
#4 - Working time arrangements - Considered full time (35 hrs) = 1
#5 - Collective organisation - Part of a union = 1
data$unemployed <- ifelse(data$employment_status_cat == "Unemployed", 1, 0)
data$nilf <- ifelse(data$employment_status_cat == "NILF", 1, 0)
data$unemp_or_nilf <- ifelse(data$unemployed == 1 | data$nilf == 1, 1, 0)

data$eq_dim_1 <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$emp_tenure_months >=12, 1, 0.5))
data$eq_dim_2_a <- ifelse(data$z_total_labor_income >= quantile(data$z_total_labor_income, c(0.5)), (1/3), 0)
data$eq_dim_2_b <- ifelse(data$emp_hi == 1, (1/3), 0)
data$eq_dim_2_c <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$pension == 1, (1/3), 0))
data$eq_dim_3_a <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$overtime == 1, 0.5, 0))
data$eq_dim_3_b <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$salaried == 1, 0.5, 0))
data$eq_dim_4 <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$avg_hours_worked_week >= 32, 1, 0)) 
data$eq_dim_5 <- ifelse(data$unemp_or_nilf == 1, 0, ifelse(data$union == 1, 1, 0))
data$linear_eq <- data$eq_dim_1 + data$eq_dim_2_a + data$eq_dim_2_b + data$eq_dim_2_c +
  data$eq_dim_3_a + data$eq_dim_3_b + data$eq_dim_4 + data$eq_dim_5
#NOTE - For each person, they'll only have a recorded linear EQ if all dependent variables are known
#       However, they'll have a known maximum and minimum linear EQ based on the known linear EQ components
#       After multiple imputation of linear_eq, imputed scores < linear_eq_min or > linear_eq_max will
#       be recoded to equal the minimum or maximum respectively
data %>% mutate(linear_eq_min = ifelse(is.na(eq_dim_1), 0, eq_dim_1) +
                  ifelse(is.na(eq_dim_2_a), 0, eq_dim_2_a) +
                  ifelse(is.na(eq_dim_2_b), 0, eq_dim_2_b) +
                  ifelse(is.na(eq_dim_2_c), 0, eq_dim_2_c) +
                  ifelse(is.na(eq_dim_3_a), 0, eq_dim_3_a) + 
                  ifelse(is.na(eq_dim_3_b), 0, eq_dim_3_b) +
                  ifelse(is.na(eq_dim_4), 0, eq_dim_4) +
                  ifelse(is.na(eq_dim_5), 0, eq_dim_5),
                linear_eq_na = ifelse(!is.na(eq_dim_1), 0, 1) +
                  ifelse(!is.na(eq_dim_2_a), 0, (1/3)) +
                  ifelse(!is.na(eq_dim_2_b), 0, (1/3)) +
                  ifelse(!is.na(eq_dim_2_c), 0, (1/3)) +
                  ifelse(!is.na(eq_dim_3_a), 0, 0.5) + 
                  ifelse(!is.na(eq_dim_3_b), 0, 0.5) +
                  ifelse(!is.na(eq_dim_4), 0, 1) +
                  ifelse(!is.na(eq_dim_5), 0, 1),
                linear_eq_max = linear_eq_min + linear_eq_na) -> data
data %>% select(unique_id, year, linear_eq_min, linear_eq_max) -> linear_eq

#### Part 2 ####

#### Creating variables needed to determine eligibility ####
#At its simplest, the following are needed for eligiblity:
#     1 - Known Baseline C, Educational attainment (A), EQ (M), and K6 (Y) at 1 wave
#     2 - A set of N consecutive-wave records including a record meeting 1.
#Notes - Conservatively, missing educational attainment can be imputed through carry-forward
#        (unless earliest known dichotomised education is <= HS where carry-backward is fine),
#        while Baseline C not already imputed in psid_download_initial_dataset_build_091621
#        (disability, region, occupation) can be imputed through carry-forward and backward 
#      - Race and ethnicity are included as a dichotomised 'POC or Hispanic vs. Neither'
#        variable. As persons could be 'POC or Hispanic' without having known information
#        on both race and ethnicity, we create a new variable before sample restriction 
data %>% 
  mutate(
    #Creating 'hs_or_less' dichotomous educational attainment variable
    hs_or_less = case_when((education == "HS" | education == "<HS") ~ 1,
                           (education == "Some college" | education == "College") ~ 0,
                           TRUE ~ NA_real_),
    #Creating indicators for known EQ, K6, both EQ and K6
    eq_known = ifelse(!is.na(eq), 1, 0),
    k6_known = ifelse(!is.na(k6), 1, 0),
    eq_k6_known = ifelse(!is.na(eq) & !is.na(k6), 1, 0),
    #Creating a 'POC or Hisanic' variable
    poc_hisp = case_when((is.na(race) & is.na(ethnicity)) ~ NA_real_,
                         (race == "White" & ethnicity == "Non-Hispanic") ~ 0,
                         ethnicity == "Hispanic" ~ 1,
                         race != "White" ~ 1,
                         TRUE ~ NA_real_)) %>%
  #Imputing educational attainment, baseline confounders through CF, CB
  group_by(unique_id) %>%
  mutate(education = na.locf(education, na.rm = F),
         hs_or_less = na.locf(hs_or_less, na.rm = F),
         hs_or_less = case_when((is.na(hs_or_less) & year < min(year[which(hs_or_less == 0)])) ~ 0, 
                                TRUE ~ hs_or_less),
         disabl_limits_work = na.locf(disabl_limits_work, na.rm = F),
         disabl_limits_work = na.locf(disabl_limits_work, na.rm = F, fromLast = T),
         region = na.locf(region, na.rm = F),
         region = na.locf(region, na.rm = F, fromLast = T),
         occupation_broad = na.locf(occupation_broad, na.rm = F),
         occupation_broad = na.locf(occupation_broad, na.rm = F, fromLast = T)) %>% ungroup() -> data

#Recoding those whose imputed occupation via CF/CB was 'Unemployed or NILF' while employed to NA
#This step affects 37 records
data$occupation_broad[data$occupation_broad == "Unemployed or NILF" & data$employment_status_cat == "Working or temp. layoff"] <- NA

#Restricting dataset to the set of persons meeting 1 in at least 1 wave
#This step removes 98,191 records, 15,709 individuals
data %>% 
  select(unique_id, year, hs_or_less, age, gender, poc_hisp, nativity, disabl_limits_work, marital_status,
         parents_poor, region, occupation_broad, eq, eq_known, k6, k6_known, eq_k6_known, SRH) %>%
  na.omit() %>% select(unique_id, year) %>% 
  mutate(elig_in_wave = 1) -> eligible_ids
data %>% filter(unique_id %in% eligible_ids$unique_id) -> data

#Creating indicator for meeting 1 in a given wave
data <- merge(data, eligible_ids, by = c("unique_id", "year"), all = T)
data$elig_in_wave[is.na(data$elig_in_wave)] <- 0
rm(eligible_ids)

#Restricting to data from 2001 onwards, as K6 will only be known from this point on
data <- data[data$year >=2001, ]

#Notes - In the primary analysis, time-varying C will include unemployment status, marital status,
#        and SRH. As PCA-based EQ requires employment, however, a complete-case analysis
#        will only include marital status and SRH as time-varying C, as all those with a known EQ 
#        must be either employed or primarily NILF (though employed) to have known EQ. 
#      - To be eligible for a complete-case analysis, persons must have known
#        time-varying C (marital status, SRH) in each wave.

#### Creating a sequence number for each person from their first wave visit >= 2001 #### 
data %>% group_by(unique_id) %>% 
  mutate(base_year = max(2001, min(year))) %>% ungroup() %>%
  mutate(seq = (year - base_year)/2) %>% group_by(unique_id) %>%
  mutate(base_year = ifelse(min(year) < 2001, NA, base_year)) %>% ungroup() %>%
  mutate(seq = ifelse(seq < 0, NA, seq)) %>% group_by(unique_id) %>% 
  mutate(waves = paste(sort(unique(seq, na.rm = T)), collapse = ""),
         cc_waves = paste0(sort(unique(seq[elig_in_wave == 1], na.rm = T)), collapse = "")) %>% ungroup() -> data

#Creating strings for each person for earliest consecutive sequences of length N (overall and within waves with complete data)
#where at least one wave in the sequence has complete exposure, mediator, outcome, confounder information
for (n in 3:6) { 
  data[[paste0("earliest_", n)]] <- NA_character_ 
  data[[paste0("earliest_cc_", n)]] <- NA_character_
}

for (i in 7:0) {
  data$earliest_3 <- ifelse(grepl(paste0(c(0:2)+i, collapse = ""), data$waves), 
                            ifelse((grepl((c(0:2)+i)[1], data$cc_waves) | 
                                      grepl((c(0:2)+i)[2], data$cc_waves) |  
                                      grepl((c(0:2)+i)[3], data$cc_waves)), paste0(c(0:2)+i, collapse = ""), data$earliest_3), data$earliest_3)
  data$earliest_cc_3 <- ifelse(grepl(paste0(c(0:2)+i, collapse = ""), data$cc_waves), paste0(c(0:2)+i, collapse = ""), data$earliest_cc_3)
}

for (i in 6:0) {
  data$earliest_4 <- ifelse(grepl(paste0(c(0:3)+i, collapse = ""), data$waves), 
                            ifelse(grepl((c(0:3)+i)[1], data$cc_waves) | 
                                     grepl((c(0:3)+i)[2], data$cc_waves) |  
                                     grepl((c(0:3)+i)[3], data$cc_waves) | 
                                     grepl((c(0:3)+i)[4], data$cc_waves), paste0(c(0:3)+i, collapse = ""), data$earliest_4), data$earliest_4)
  data$earliest_cc_4 <- ifelse(grepl(paste0(c(0:3)+i, collapse = ""), data$cc_waves), paste0(c(0:3)+i, collapse = ""), data$earliest_cc_4)
}
for (i in 5:0) {
  data$earliest_5 <- ifelse(grepl(paste0(c(0:4)+i, collapse = ""), data$waves), 
                            ifelse(grepl((c(0:4)+i)[1], data$cc_waves) | 
                                     grepl((c(0:4)+i)[2], data$cc_waves) |  
                                     grepl((c(0:4)+i)[3], data$cc_waves) | 
                                     grepl((c(0:4)+i)[4], data$cc_waves) |
                                     grepl((c(0:4)+i)[5], data$cc_waves), paste0(c(0:4)+i, collapse = ""), data$earliest_5), data$earliest_5)
  data$earliest_cc_5 <- ifelse(grepl(paste0(c(0:4)+i, collapse = ""), data$cc_waves), paste0(c(0:4)+i, collapse = ""), data$earliest_cc_5)
}
for (i in 4:0) {
  data$earliest_6 <- ifelse(grepl(paste0(c(0:5)+i, collapse = ""), data$waves), 
                            ifelse(grepl((c(0:5)+i)[1], data$cc_waves) | 
                                     grepl((c(0:5)+i)[2], data$cc_waves) |  
                                     grepl((c(0:5)+i)[3], data$cc_waves) | 
                                     grepl((c(0:5)+i)[4], data$cc_waves) |
                                     grepl((c(0:5)+i)[5], data$cc_waves) |
                                     grepl((c(0:5)+i)[6], data$cc_waves), paste0(c(0:5)+i, collapse = ""), data$earliest_6), data$earliest_6)
  data$earliest_cc_6 <- ifelse(grepl(paste0(c(0:5)+i, collapse = ""), data$cc_waves), paste0(c(0:5)+i, collapse = ""), data$earliest_cc_6)
}
rm(n, i)

#Creating indicators for whether each record is included within a person's N-length sequence, for easy subsetting
data$seq_3 <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_3[x]), 1, 0) } )
data$seq_3_cc <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_cc_3[x]), 1, 0) } )
data$seq_4 <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_4[x]), 1, 0) } )
data$seq_4_cc <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_cc_4[x]), 1, 0) } )
data$seq_5 <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_5[x]), 1, 0) } )
data$seq_5_cc <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_cc_5[x]), 1, 0) } )
data$seq_6 <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_6[x]), 1, 0) } )
data$seq_6_cc <- sapply(1:NROW(data), FUN = function(x) { ifelse(grepl(as.character(data$seq[x]), data$earliest_cc_6[x]), 1, 0) } )

#Creating record subsets for each analysis
#Primary analyses
data %>% filter(year >= 2001 & seq_3 == 1) %>% select(unique_id, year) -> analysis_3
data %>% filter(year >= 2001 & seq_4 == 1) %>% select(unique_id, year) -> analysis_4
data %>% filter(year >= 2001 & seq_5 == 1) %>% select(unique_id, year) -> analysis_5
data %>% filter(year >= 2001 & seq_6 == 1) %>% select(unique_id, year) -> analysis_6
#Complete case analyses
data %>% filter(year >= 2001 & seq_3_cc == 1) %>% select(unique_id, year) -> cc_analysis_3
data %>% filter(year >= 2001 & seq_4_cc == 1) %>% select(unique_id, year) -> cc_analysis_4
data %>% filter(year >= 2001 & seq_5_cc == 1) %>% select(unique_id, year) -> cc_analysis_5
data %>% filter(year >= 2001 & seq_6_cc == 1) %>% select(unique_id, year) -> cc_analysis_6

#### Part 3 ####
#Notes - Following the White et al. (2011) rule of thumb, the number of imputations
#        should be at least as large as the % of incomplete cases. Equally, based
#        on Sullivan et al. (2018), multiple imputation should be performed separately
#        by exposure. Likewise, as multiple imputation performs best when using the maximum
#        available sample, MI should be performed on the entire eligible sample before
#        subsetting to the analysis-specific samples.
#      - Based on 'elig_in_wave', for each primary analysis sample, the percentage of
#        incomplete cases ranges from 36.8% to 37.7%. Therefore, 40 imputations will be used.
#      - Though persons could have different visit subsets used in different analyses
#        (e.g. 0,1,2 for seq_3, and 4,5,6,7 for seq_4) where baseline visits (e.g. 0,4 for left)
#        had different exposure levels (e.g. <= HS at 0, but >HS at 4), hs_or_less status
#        changes across analysis baseline visits for 42 persons (0.7%), so it may be safe to perform 
#        multiple imputation once stratifying by baseline education for the 3-timepoint analysis
#        instead of separately for each length analysis. 

#Creating indicators for baseline exposure in each analysis
data %>% group_by(unique_id) %>% 
  mutate(base_education_3 = ifelse(is.na(seq_3), NA, hs_or_less[year == min(year[which(seq_3 == 1)])]),
         base_education_4 = ifelse(is.na(seq_4), NA, hs_or_less[year == min(year[which(seq_4 == 1)])]),
         base_education_5 = ifelse(is.na(seq_5), NA, hs_or_less[year == min(year[which(seq_5 == 1)])]),
         base_education_6 = ifelse(is.na(seq_6), NA, hs_or_less[year == min(year[which(seq_6 == 1)])])) %>% ungroup() -> data

#Checking differences in baseline education across primary analysis samples
data %>% filter((base_education_3 != base_education_4 & !is.na(base_education_4)) | 
                  (base_education_3 != base_education_5 & !is.na(base_education_5)) | 
                  (base_education_3 != base_education_6 & !is.na(base_education_6))) %>% 
  select(unique_id, year, base_education_3, base_education_4, base_education_5, base_education_6) %>% View()

#Theory-based imputation of EQ components for those unemployed, self-employed before MI
#Notes - We make the following assumptions about EQ components for those:
#        Unemployed or NILF - 
#            Union (No), Employer HI (No assumptions), Salaried (No), Overtime (No), 
#            Pension (No), Annual Hours (No assumptions), Unemployment Duration (No assumptions), 
#            Z-Score Total Labor Income (No assumptions),
#            Z-Score Emp Tenure (assign z_emp_tenure_months corresponding to emp_tenure_months = 0)
#        Partially self-employed - 
#            No assumptions made, as all EQ components can depend on their non-self-employed role 
#        Fully self-employed - 
#            Union (No), Employer HI (No assumptions), Salaried (No), Overtime (No), Pension (No assumptions), 
#            Annual Hours (No assumptions), Unemployment Duration (No assumptions), 
#            Z-Score Total Labor Income (No assumptions), Z-Score Emp Tenure (No assumptions)

data %>% mutate(union = case_when(employment_status_cat == "Unemployed" | employment_status_cat == "NILF" ~ 0,
                                  TRUE ~ union),
                salaried = case_when(employment_status_cat == "Unemployed" | employment_status_cat == "NILF"  ~ 0,
                                     TRUE ~ salaried),
                overtime = case_when(employment_status_cat == "Unemployed" | employment_status_cat == "NILF"  ~ 0,
                                     TRUE ~ overtime),
                pension = case_when(employment_status_cat == "Unemployed" | employment_status_cat == "NILF"  ~ 0,
                                    TRUE ~ pension),
                z_emp_tenure_months = case_when(employment_status_cat == "Unemployed" | employment_status_cat == "NILF"  ~ -1.35,
                                                TRUE ~ z_emp_tenure_months)) -> data

#### Recoding variables as they will be used in analyses ####
#Notes - Only 209 records have US Territory or Foreign Country for region, which could cause convergence issues. 
#        Checking reported EQ by region, this group most closely resembles the South, so are recoded as 'South'
#      - From prior situations, string variables with commas or slashes can affect the 'miceadds' model functions.
#        To prevent this, occupation_broad is recoded to remove spaces, dashes, commas
#      - Variables on different scales can also affect model convergence. For this reason, 'annual_hours_worked'
#        and 'year' are recoded 1) annual hours worked in 1000s, and 2) years since 2001

data %>% mutate(female = ifelse(gender == "Female", 1, 0),
                married = case_when(marital_status == "Married/Cohabiting" ~ 1, !is.na(marital_status) ~ 0, TRUE ~ NA_real_),
                nativity_not_US = ifelse(nativity == "US", 0, 1),
                srh_vgood_exc = case_when(SRH == "Very good" | SRH == "Excellent" ~ 1, !is.na(SRH) ~ 0, TRUE ~ NA_real_),
                parents_poor = ifelse(parents_poor == "Poor", 1, 0),
                region_revised = case_when(region == "US Territory or Foreign Country" ~ "South", TRUE ~ region),
                occupation_broad = str_replace_all(occupation_broad, "[ ,/]", ""),
                year_2001 = year - 2001,
                annual_hours_worked_revised = annual_hours_worked/1000) -> data

#Creating a baseline education variable to split multiple imputation by
#The value of base_ed will be hs_or_less from the first year where seq_3 == 1
data %>% group_by(unique_id) %>% mutate(min_year = ifelse(mean(seq_3) == 0, min(year), min(year[seq_3==1]))) %>% ungroup() -> data
data %>% group_by(unique_id) %>% mutate(base_ed = hs_or_less[year == min_year]) %>% ungroup() -> data

#Restricting to necessary variables for imputation and/or mediational g-formula
data %>% select(#Demographic factors
  unique_id, year_2001, age, female, poc_hisp, hs_or_less, base_ed, 
  nativity_not_US, parents_poor, region_revised, hs_or_less, married, 
  srh_vgood_exc, disabl_limits_work, occupation_broad, k6,
  #Employment-related factors
  eq, linear_eq, unemployed, nilf, union, emp_hi, 
  salaried, overtime, pension, annual_hours_worked_revised, 
  unemp_duration_months, z_total_labor_income, z_emp_tenure_months) -> data

#### Saving complete case analysis samples before multiple imputation ####
#Setting new directory to save created datasets
directory_analysis <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"

#Reformatting and saving complete case analysis samples
#Creating samples
data$year <- data$year_2001+2001 #Needed for merging
cc_data_3 <- merge(data, cc_analysis_3, by = c("unique_id", "year"), all.x = F, all.y = T)
cc_data_4 <- merge(data, cc_analysis_4, by = c("unique_id", "year"), all.x = F, all.y = T)
cc_data_5 <- merge(data, cc_analysis_5, by = c("unique_id", "year"), all.x = F, all.y = T)
cc_data_6 <- merge(data, cc_analysis_6, by = c("unique_id", "year"), all.x = F, all.y = T)
data$year <- NULL

#Creating baseline variables as needed, binary K6, single 'unemp_or_nilf' dummy
cc_data_3 %>% mutate(k6_bin = ifelse(k6 >= 5, 1, 0)) -> cc_data_3
cc_data_4 %>% mutate(k6_bin = ifelse(k6 >= 5, 1, 0)) -> cc_data_4
cc_data_5 %>% mutate(k6_bin = ifelse(k6 >= 5, 1, 0)) -> cc_data_5
cc_data_6 %>% mutate(k6_bin = ifelse(k6 >= 5, 1, 0)) -> cc_data_6

cc_data_3 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                             base_age = age[year == min(year)],
                                             base_region = region_revised[year == min(year)],
                                             base_disability = disabl_limits_work[year == min(year)],
                                             base_occupation = occupation_broad[year == min(year)],
                                             base_ed = hs_or_less[year == min(year)],
                                             seq = (year - base_year)/2) %>% ungroup() %>%
  select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad) -> cc_data_3
cc_data_4 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                             base_age = age[year == min(year)],
                                             base_region = region_revised[year == min(year)],
                                             base_disability = disabl_limits_work[year == min(year)],
                                             base_occupation = occupation_broad[year == min(year)],
                                             base_ed = hs_or_less[year == min(year)],
                                             seq = (year - base_year)/2) %>% ungroup() %>%
  select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad) -> cc_data_4
cc_data_5 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                             base_age = age[year == min(year)],
                                             base_region = region_revised[year == min(year)],
                                             base_disability = disabl_limits_work[year == min(year)],
                                             base_occupation = occupation_broad[year == min(year)],
                                             base_ed = hs_or_less[year == min(year)],
                                             seq = (year - base_year)/2) %>% ungroup() %>%
  select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad) -> cc_data_5
cc_data_6 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                             base_age = age[year == min(year)],
                                             base_region = region_revised[year == min(year)],
                                             base_disability = disabl_limits_work[year == min(year)],
                                             base_occupation = occupation_broad[year == min(year)],
                                             base_ed = hs_or_less[year == min(year)],
                                             seq = (year - base_year)/2) %>% ungroup() %>%
  select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad) -> cc_data_6

#Restricting to necessary variables in each case
cc_data_3 %>% select(unique_id, seq, base_year, base_age, female, poc_hisp, base_ed, 
                     nativity_not_US, parents_poor, base_region, base_disability, 
                     base_occupation, married, srh_vgood_exc, k6, k6_bin, eq, linear_eq) -> cc_data_3
cc_data_4 %>% select(unique_id, seq, base_year, base_age, female, poc_hisp, base_ed, 
                     nativity_not_US, parents_poor, base_region, base_disability, 
                     base_occupation, married, srh_vgood_exc, k6, k6_bin, eq, linear_eq) -> cc_data_4
cc_data_5 %>% select(unique_id, seq, base_year, base_age, female, poc_hisp, base_ed, 
                     nativity_not_US, parents_poor, base_region, base_disability, 
                     base_occupation, married, srh_vgood_exc, k6, k6_bin, eq, linear_eq) -> cc_data_5
cc_data_6 %>% select(unique_id, seq, base_year, base_age, female, poc_hisp, base_ed, 
                     nativity_not_US, parents_poor, base_region, base_disability, 
                     base_occupation, married, srh_vgood_exc, k6, k6_bin, eq, linear_eq) -> cc_data_6

#Turning each dataset back from 'tibble' to 'data.frame'
cc_data_3 <- as.data.frame(cc_data_3)
cc_data_4 <- as.data.frame(cc_data_4)
cc_data_5 <- as.data.frame(cc_data_5)
cc_data_6 <- as.data.frame(cc_data_6)

#Saving long-format datasets
write.csv(cc_data_3, paste0(directory_analysis, "dataset_100121_long_3_cc.csv"), row.names = F, na = "")
write.csv(cc_data_4, paste0(directory_analysis, "dataset_100121_long_4_cc.csv"), row.names = F, na = "")
write.csv(cc_data_5, paste0(directory_analysis, "dataset_100121_long_5_cc.csv"), row.names = F, na = "")
write.csv(cc_data_6, paste0(directory_analysis, "dataset_100121_long_6_cc.csv"), row.names = F, na = "")

#Creating wide-format datasets
cc_data_3_wide <- reshape(cc_data_3, direction = "wide", timevar = "seq", sep = "_",
                          idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                    "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
cc_data_4_wide <- reshape(cc_data_4, direction = "wide", timevar = "seq", sep = "_",
                          idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                    "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
cc_data_5_wide <- reshape(cc_data_5, direction = "wide", timevar = "seq", sep = "_",
                          idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                    "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
cc_data_6_wide <- reshape(cc_data_6, direction = "wide", timevar = "seq", sep = "_",
                          idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                    "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
#Saving wide-format datasets
write.csv(cc_data_3_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_3_cc.csv"), row.names = F, na = "")
write.csv(cc_data_4_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_4_cc.csv"), row.names = F, na = "")
write.csv(cc_data_5_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_5_cc.csv"), row.names = F, na = "")
write.csv(cc_data_6_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_6_cc.csv"), row.names = F, na = "")

rm(cc_data_3, cc_data_4, cc_data_5, cc_data_6, cc_data_3_wide, cc_data_4_wide, cc_data_5_wide, cc_data_6_wide)

#Converting occupation to factor to avoid imputing empty strings
data$occupation_broad <- as.factor(data$occupation_broad)

#### Creating 40 Multiply Imputed Datasets ####
pred_matrix <- make.predictorMatrix(data)
pred_matrix[, c(1,7)] <- 0      #Excluding ID and baseline education from predictor list 
pred_matrix[c(15:17), 1] <- -2  #Clustering EQ and K6 by ID

#Performing imputation by exposure then remerging exposure strata to create complete datasets
mi_type <- mice(data, maxit = 0)$method
mi_type[14] <- "lda"
mi_type[c(15:17)] <- "2l.pmm"

data <- as.data.frame(data)
data_hs_less <- data[data$base_ed == 1, ]
data_more_hs <- data[data$base_ed == 0, ]
mi_hs_less <- mice(data_hs_less, method = mi_type, pred = pred_matrix, seed = 2021, maxit = 25, m = 40, printFlag = TRUE)
mi_more_hs <- mice(data_more_hs, method = mi_type, pred = pred_matrix, seed = 2021, maxit = 25, m = 40, printFlag = TRUE)

#Saving each mids object
saveRDS(mi_hs_less, file = paste0(directory_analysis, "MI_data/mids_hsless_100121.rds"))
saveRDS(mi_more_hs, file = paste0(directory_analysis, "MI_data/mids_morehs_100121.rds"))

# Combining multiply imputed strata datasets, creating binary K6
for (i in 1:40) {
  hs_less <- complete(mi_hs_less, i)
  more_hs <- complete(mi_more_hs, i)
  imp <- rbind(hs_less, more_hs)
  imp <- imp[, c(1:19)]
  imp$k6_bin <- ifelse(imp$k6 >= 5, 1, 0)
  assign(paste0("imp_", i), imp)
  rm(hs_less, more_hs, imp)
}
rm(i, data_hs_less, data_more_hs, pred_matrix, mi_type)

#### Final reformatting post-imputation for each dataset ####
#Steps - 1 - For each analysis length, the imputed dataset is restricted to records 
#            used in that analysis (merging with e.g. analysis_3 for 6-year analyses)
#        2 - Create an 'unemp_or_nilf' dummy variable
#        3 - Recode linear_eq to min or max (if beyond) based on 'linear_eq' dataset
#        4 - Creating baseline variables for year, age, region, disability, occupation, and education are created
#        5 - Long- and wide-format versions of the dataset are saved

count <- 0
for (i in grep(paste0("imp_"), ls(), value = T)) {
  #Loading data
  data <- eval(parse(text = i))
  n <- gsub(paste0("imp_"), "", i)
  
  #Counter output
  count <- count + 1
  print(paste0("Starting ", count, " - ", round(count/40*100, 2), "% Complete"))
  
  #Recreating 'year' variable
  data$year <- data$year_2001 + 2001
  
  #Restricting to needed records per analysis
  data_3 <- merge(data, analysis_3, by = c("unique_id", "year"), all.x = F, all.y = T)
  data_4 <- merge(data, analysis_4, by = c("unique_id", "year"), all.x = F, all.y = T)
  data_5 <- merge(data, analysis_5, by = c("unique_id", "year"), all.x = F, all.y = T)
  data_6 <- merge(data, analysis_6, by = c("unique_id", "year"), all.x = F, all.y = T)
  
  #Create 'unemp_or_nilf' dummy variable
  data_3 %>% mutate(unemp_or_nilf = ifelse(unemployed == 1 | nilf == 1, 1, 0)) -> data_3
  data_4 %>% mutate(unemp_or_nilf = ifelse(unemployed == 1 | nilf == 1, 1, 0)) -> data_4
  data_5 %>% mutate(unemp_or_nilf = ifelse(unemployed == 1 | nilf == 1, 1, 0)) -> data_5
  data_6 %>% mutate(unemp_or_nilf = ifelse(unemployed == 1 | nilf == 1, 1, 0)) -> data_6
  
  #Recoding linear_eq to min or max linear_eq if outside range
  data_3 <- merge(data_3, linear_eq, by = c("unique_id", "year"), all.x = T, all.y = F)
  data_4 <- merge(data_4, linear_eq, by = c("unique_id", "year"), all.x = T, all.y = F)
  data_5 <- merge(data_5, linear_eq, by = c("unique_id", "year"), all.x = T, all.y = F)
  data_6 <- merge(data_6, linear_eq, by = c("unique_id", "year"), all.x = T, all.y = F)
  
  data_3 %>% mutate(linear_eq = case_when(linear_eq < linear_eq_min ~ linear_eq_min,
                                          linear_eq > linear_eq_max ~ linear_eq_max,
                                          TRUE ~ linear_eq)) -> data_3
  data_4 %>% mutate(linear_eq = case_when(linear_eq < linear_eq_min ~ linear_eq_min,
                                          linear_eq > linear_eq_max ~ linear_eq_max,
                                          TRUE ~ linear_eq)) -> data_4
  data_5 %>% mutate(linear_eq = case_when(linear_eq < linear_eq_min ~ linear_eq_min,
                                          linear_eq > linear_eq_max ~ linear_eq_max,
                                          TRUE ~ linear_eq)) -> data_5
  data_6 %>% mutate(linear_eq = case_when(linear_eq < linear_eq_min ~ linear_eq_min,
                                          linear_eq > linear_eq_max ~ linear_eq_max,
                                          TRUE ~ linear_eq)) -> data_6
  
  #Creating baseline variables per dataset
  data_3 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                            base_age = age[year == min(year)],
                                            base_region = region_revised[year == min(year)],
                                            base_disability = disabl_limits_work[year == min(year)],
                                            base_occupation = occupation_broad[year == min(year)],
                                            base_ed = hs_or_less[year == min(year)],
                                            seq = (year - base_year)/2) %>% ungroup() %>%
    select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad, -linear_eq_min, -linear_eq_max) -> data_3
  data_4 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                            base_age = age[year == min(year)],
                                            base_region = region_revised[year == min(year)],
                                            base_disability = disabl_limits_work[year == min(year)],
                                            base_occupation = occupation_broad[year == min(year)],
                                            base_ed = hs_or_less[year == min(year)],
                                            seq = (year - base_year)/2) %>% ungroup() %>%
    select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad, -linear_eq_min, -linear_eq_max) -> data_4
  data_5 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                            base_age = age[year == min(year)],
                                            base_region = region_revised[year == min(year)],
                                            base_disability = disabl_limits_work[year == min(year)],
                                            base_occupation = occupation_broad[year == min(year)],
                                            base_ed = hs_or_less[year == min(year)],
                                            seq = (year - base_year)/2) %>% ungroup() %>%
    select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad, -linear_eq_min, -linear_eq_max) -> data_5
  data_6 %>% group_by(unique_id) %>% mutate(base_year = year[year == min(year)],
                                            base_age = age[year == min(year)],
                                            base_region = region_revised[year == min(year)],
                                            base_disability = disabl_limits_work[year == min(year)],
                                            base_occupation = occupation_broad[year == min(year)],
                                            base_ed = hs_or_less[year == min(year)],
                                            seq = (year - base_year)/2) %>% ungroup() %>%
    select(-year_2001, -year, -age, -hs_or_less, -region_revised, -disabl_limits_work, -occupation_broad, -linear_eq_min, -linear_eq_max) -> data_6
  
  #Turning each dataset back from 'tibble' to 'data.frame'
  data_3 <- as.data.frame(data_3)
  data_4 <- as.data.frame(data_4)
  data_5 <- as.data.frame(data_5)
  data_6 <- as.data.frame(data_6)
  
  #Saving long-format datasets
  write.csv(data_3, paste0(directory_analysis, "MI_data/dataset_100121_long_3_", n, ".csv"), row.names = F, na = "")
  write.csv(data_4, paste0(directory_analysis, "MI_data/dataset_100121_long_4_", n, ".csv"), row.names = F, na = "")
  write.csv(data_5, paste0(directory_analysis, "MI_data/dataset_100121_long_5_", n, ".csv"), row.names = F, na = "")
  write.csv(data_6, paste0(directory_analysis, "MI_data/dataset_100121_long_6_", n, ".csv"), row.names = F, na = "")
  
  #Creating wide-format datasets
  data_3_wide <- reshape(data_3, direction = "wide", timevar = "seq", sep = "_",
                         idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                   "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
  data_4_wide <- reshape(data_4, direction = "wide", timevar = "seq", sep = "_",
                         idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                   "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
  data_5_wide <- reshape(data_5, direction = "wide", timevar = "seq", sep = "_",
                         idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                   "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
  data_6_wide <- reshape(data_6, direction = "wide", timevar = "seq", sep = "_",
                         idvar = c("unique_id", "base_year", "base_age", "female", "poc_hisp", "nativity_not_US", 
                                   "parents_poor", "base_ed", "base_region", "base_disability", "base_occupation"))
  #Saving wide-format datasets
  write.csv(data_3_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_3_", n, ".csv"), row.names = F, na = "")
  write.csv(data_4_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_4_", n, ".csv"), row.names = F, na = "")
  write.csv(data_5_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_5_", n, ".csv"), row.names = F, na = "")
  write.csv(data_6_wide, paste0(directory_analysis, "MI_data/dataset_100121_wide_6_", n, ".csv"), row.names = F, na = "")
  
  rm(data, n, data_3, data_4, data_5, data_6, data_3_wide, data_4_wide, data_5_wide, data_6_wide)
}
rm(list = c(grep(paste0("imp_"), ls(), value = T), i, count))
