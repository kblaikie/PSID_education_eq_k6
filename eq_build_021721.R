### Build Script - EQ Score Dataset 2/17/21 ####
#Libraries #
library(tidyverse)
library(haven)
library(foreign)
library(zoo) #For carry forward/backward NA

#Loading data #
directory <- "R:/Project/precarityR01/PSID/analyses/Mediation/Kieran/Data/"
data <- read.csv(paste0(directory, "R_data_5_18_20.csv"))

#Restricting to those of working age #
data <- data[data$age_unique >=30 & data$age_unique <=60 & !is.na(data$unique_id), ]

#Reformatting/creating variables as needed for PCA #
data$union <- ifelse(data$current_job_union_member == "Yes", 1, 0)
data$empl_HI_bin <- ifelse(data$empl_hlth_ins_99_17 == "Yes", 1, 0)
data$salaried <- ifelse(data$wage_vs_salary_revised == "Salaried", 1, 
                        ifelse(data$wage_vs_salary_revised == "Waged"|data$wage_vs_salary_revised == "Other", 0, NA))
data$paid_extra_overtime <- NA #Those wage_vs_salary_revised ('Inapp.', NA) and paid_for_overtime_main_job_waged ('Not asked', NA) 
data$paid_extra_overtime[data$wage_vs_salary_revised=="Salaried"] <- 1 #All salaried counted as paid extra
data$paid_extra_overtime[data$wage_vs_salary_revised=="Other"] <- 0 #All 'Other' counted as not paid extra
data$paid_extra_overtime[data$wage_vs_salary_revised=="Waged" & (data$paid_for_overtime_main_job_waged=="1.5x"|
                                                                 data$paid_for_overtime_main_job_waged=="2x"|
                                                                 data$paid_for_overtime_main_job_waged=="Exact amount"|
                                                                 data$paid_for_overtime_main_job_waged=="Other")] <- 1
data$paid_extra_overtime[data$wage_vs_salary_revised=="Waged" & (data$paid_for_overtime_main_job_waged=="1x"|
                                                                 data$paid_for_overtime_main_job_waged=="IAP"|
                                                                 data$paid_for_overtime_main_job_waged=="Comp time")] <- 0
data$employment_status_tri <- ifelse(data$employment_status=="Employed/temp. layoff", "Employed/temp. layoff",
                              ifelse(data$employment_status=="Unemployed", "Unemployed", 
                              ifelse(data$employment_status %in% c("Other", "Keeping house", "Retired or disabled", "Student"), "Keeping house/student/retired/disabled/other", NA)))
data$pension <- ifelse(data$have_pension_retirement_plan_current_job == "Yes", 1, 0)
data$pension <- ifelse(data$have_pension_retirement_plan_current_job == "Inapp."|data$have_pension_retirement_plan_current_job == "Not asked", NA, data$pension)
data$unempl_duration_mnths_revised <- as.numeric(data$unempl_duration_mnths_revised) #Coerces 'Inapp.' to NA
data$empl_tenure_mnths_revised <- as.numeric(data$empl_tenure_mnths_revised) #Coerces 'Inapp.: self-employed only' to NA
data$family_income <- data$family_income/10000
data$occ1990dd_broad <- ifelse(data$employment_status_tri == "Keeping house/student/retired/disabled/other"|data$employment_status_tri=="Unemployed", data$employment_status_tri,
                                   ifelse(data$occ1990dd_broad == "Inapp.", NA, data$occ1990dd_broad))
data$ind1990_broad <- ifelse(data$employment_status_tri == "Keeping house/student/retired/disabled/other"|data$employment_status_tri=="Unemployed", data$employment_status_tri,
                                 ifelse(data$ind1990_broad == "Inapp.", NA, data$ind1990_broad))
cpi <- read.csv(paste0(directory, "cpi.csv")) #Adjustment for income and wealth to 2017 dollars 
cpi %>%
  right_join(data, cpi, by="year") %>%
  mutate(ratio=cpi/245.1, #ratio of cpi in current year to cpi in 2017
         total_labor_income=total_labor_income/ratio,
         family_income=family_income/ratio,
         wealth_no_home_equity=as.numeric(as.character(wealth_no_home_equity))/ratio,
         wealth_with_home_equity=as.numeric(as.character(wealth_with_home_equity))/ratio) -> data
rm(cpi)

#Carry-forward then carry-backward imputation 
data %>%  
  arrange(unique_id, age_unique) %>%
  group_by(unique_id) %>%
  mutate(occ1990dd_broad=na.locf(occ1990dd_broad, na.rm=FALSE),
         occ1990dd_broad=na.locf(occ1990dd_broad, na.rm=FALSE, fromLast=TRUE),
         ind1990_broad=na.locf(ind1990_broad, na.rm=FALSE),
         ind1990_broad=na.locf(ind1990_broad, na.rm=FALSE, fromLast=TRUE)) %>%
ungroup() -> data

data %>% #Working out year first employed within sample
  filter(employment_status_tri=="Employed/temp. layoff" & self_employed=="Someone else") %>%
  group_by(unique_id) %>%  
  summarise(first_worker_year=min(year)) %>%
  ungroup() -> first_worker
data <- data[data$unique_id %in% unique(first_worker$unique_id), ] #3906 removed either through never working or being self-employed

data %>% #age-standardizing labor income and employment tenure
  group_by(age_unique) %>%
  summarise(mean_age_inc=mean(total_labor_income,na.rm=T),
            sd_age_inc=sd(total_labor_income,na.rm=T),
            mean_age_emp_tenure=mean(empl_tenure_mnths_revised,na.rm=T),
            sd_age_emp_tenure=sd(empl_tenure_mnths_revised,na.rm=T)) %>%
  right_join(data, by='age_unique') %>%
  mutate(z_score_total_labor_income=(total_labor_income-mean_age_inc)/sd_age_inc,
         z_score_empl_tenure_mnths_revised=(empl_tenure_mnths_revised-mean_age_emp_tenure)/sd_age_emp_tenure) -> data
data <- data[, c("unique_id",                        "year",                       "age_unique", 
                 "union",                            "empl_HI_bin",                "salaried", 
                 "paid_extra_overtime",              "pension",                    "annual_hours_worked", 
                 "unempl_duration_mnths_revised",    "z_score_total_labor_income", "z_score_empl_tenure_mnths_revised",
                 "occ1990dd_broad",                  "ind1990_broad",              "wealth_no_home_equity", 
                 "employment_status_tri",            "self_employed")]

data %>% #Cuts out those not employed for the full period 
  filter(occ1990dd_broad != "Unemployed" & occ1990dd_broad != "Keeping house/student/retired/disabled/other" & 
           ind1990_broad != "Unemployed" & ind1990_broad != "Keeping house/student/retired/disabled/other") -> data_full_emp
data_full_emp %>% select(-c("occ1990dd_broad", "ind1990_broad", "wealth_no_home_equity", "employment_status_tri", "age_unique", "self_employed")) -> data_full_emp
data_full_emp <- na.omit(data_full_emp) #Removes all those with missing data for PCA

write.dta(data_full_emp, paste0(directory, "eq_build.dta"))

#Creating an EQ score for all #
#   use "R:\Project\precarityR01\PSID\analyses\Mediation\Kieran\Data\eq_build.dta", clear
#   pca union empl_HI_bin salaried paid_extra_overtime pension annual_hours_worked unempl_dur z_score_total_labor_income z_score_empl_tenure_mnths_revsd, comp(3)
#   rotate, promax horst
#   predict comp1 comp2 comp3, score
#   gen eq_score = (comp1*.1941) + (comp2*.1784) + (comp3*.1609)
#   save "R:\Project\precarityR01\PSID\analyses\Mediation\Kieran\Data\eq_all.dta", replace

#Adding employment quality back into larger dataset
eq_all <- read.dta13(paste0(directory, "eq_all.dta"))
eq_all <- eq_all[,c("unique_id","year","eq_score")] 
data <- merge(data, eq_all, by=c("unique_id", "year"), all=T)
data$eq_known <- ifelse(!is.na(data$eq_score), 1, 0)
data <- data[data$year >=2001, ] #Restricting to study years of interest

#Instead of assigning a single bad EQ or good EQ value for all those unemployed based on wealth, 
#I considered assigning EQ based on % in wealth eCDF (e.g. 10% wealth, assign 10% EQ)
data$wealth_percentile <- as.numeric(NA)
for (i in unique(data$age_unique)) {
  data$wealth_percentile <- ecdf(data$wealth_no_home_equity)(data$wealth_no_home_equity)
}
data$eq_score_imp_pct <- ifelse(!is.na(data$eq_score), data$eq_score, quantile(data$eq_score,c(data$wealth_percentile), na.rm=T))

#Should people 'keeping house/students/retired/disabled/other' have their last EQ carried forward?
data$eq_score_imp_lcf <- data$eq_score
data %>%  
  arrange(unique_id, age_unique) %>%
  group_by(unique_id) %>%
  fill(eq_score_imp_lcf, .direction = "down") %>%
ungroup() -> data
data$eq_score_imp_lcf[data$occ1990dd_broad == "Unemployed"|data$ind1990_broad == "Unemployed"|
                      data$self_employed=="Self only"|data$self_employed=="Both someone else and self"] <- NA
table(!is.na(data$eq_score_imp_lcf), data$eq_known)


#I also made this age-specific, given wealth needs when e.g. 25 are likely not equal to those e.g. 50
data$wealth_age_percentile <- as.numeric(NA)
for (i in unique(data$age_unique)) {
  data[data$age_unique==i,]$wealth_age_percentile <- ecdf(data[data$age_unique==i,]$wealth_no_home_equity)(data[data$age_unique==i,]$wealth_no_home_equity)
}
data$eq_score_imp_pct_age <- ifelse(!is.na(data$eq_score), data$eq_score, quantile(data$eq_score,c(data$wealth_age_percentile), na.rm=T))
View(data[,c("unique_id","year", "age_unique", "employment_status_tri", "wealth_percentile","wealth_age_percentile","eq_score","eq_score_imp_pct", "eq_score_imp_pct_age")])

#Thoughts - %wealth and %EQ won't be a 1-for-1, especially if unemployed. Worth altering some way?
#Saving as current
write.csv(data, paste0(directory, "eq_imputed.csv"), col.names = T, row.names = F, na = "")
