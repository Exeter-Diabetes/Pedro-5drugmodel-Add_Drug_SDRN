
# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval)
library(viridis)
library(patchwork)


# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

source("Closedtest_functions/closedtest_continuous_function.R")
source("CPRD/03.impute_missingness.R")
source("CPRD/04.model_predictions.R")


# load model
load("fivedrugmodel_5knot_share_20230823.Rdata")


# set up dataset
analysis_post_2019_raw <- analysis_post_2019_raw %>%
  analysis$cached("analysis_post_2019") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)



# Pre-processing dataset ########################################

## Post 2019-01-01 ----
analysis_post_2019 <- analysis_post_2019_raw %>%
  mutate(
    pated = paste(patid, dstartdate, drug_class, sep = "."),
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    agetx = dstartdate_age,
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),
    
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
    smoke = factor(smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),
    imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
    imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5), labels = c("1 (least)", "2", "3", "4", "5 (most)")),
    
    ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
    ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
    ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),
    drugline = ifelse(drugline_all > 5, 5, drugline_all),
    drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),
    drugclass = drug_class
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  select(all_of(c(
    "pated", "agetx", "sex", "t2dmduration", "ethnicity",
    "drug_substance", "drug_class",
    "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  )
  )) %>%
  rename("drugclass" = "drug_class") %>%
  as.data.frame()



# Missing data imputation ######################################

## Post 2019-01-01 ----
### By mice
analysis_post_2019_mice_imputation <- imputation_methods(data = analysis_post_2019,
                                                         method = "mice", 
                                                         mice.m = 20, 
                                                         mice.maxit = 10,
                                                         mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))




# Predictions from datasets #########################################

## Post 2019-01-01 ----
### Imputation by mice
analysis_post_2019_prediction_mice <- predict_5drugmodel(analysis_post_2019_mice_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"))


### merge impute columns into main dataset
analysis_post_2019 <- analysis_post_2019 %>%
  select(-c("prebmi", "preegfr", "pretotalcholesterol", "prehdl", "prealt")) %>%
  cbind(
    analysis_post_2019_prediction_mice %>%
      select(contains("mice_impute")) %>%
      rename_with(~ str_replace(., "_mice_impute", "")),
    analysis_post_2019_prediction_mice %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred."))
  )


# Split cohort into old vs new drugs #########################################

analysis_standard <- analysis_post_2019 %>%
  filter(!(drug_substance %in% c("Oral semaglutide", "Low-dose semaglutide", "Semaglutide, dose unclear", "High-dose semaglutide")))

analysis_oral <- analysis_post_2019 %>%
  filter(drug_substance %in% c("Oral semaglutide"))

analysis_injectable <- analysis_post_2019 %>%
  filter(drug_substance %in% c("Low-dose semaglutide", "Semaglutide, dose unclear", "High-dose semaglutide"))

n_train <- floor(0.2 * nrow(analysis_injectable))
set.seed(123)

analysis_injectable <- analysis_injectable %>%
  mutate(split = ifelse(row_number() %in% sample(nrow(analysis_injectable), n_train), "training", "testing"))



# Closed loop test #################################################

## Injectable Semaglutide ----

### mice variables
#### GLP1
closed_loop_test_injectable_semaglutide <- closedtest_continuous_function(
  cohort = "Injectable semaglutide cohort",
  dataset = analysis_injectable %>%
    filter(split == "training"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

closed_loop_test_injectable_semaglutide_adjustment <- 7.3

## Adjustment used: https://doi.org/10.1016/j.diabres.2021.108904
## Difference between Fig5 Dulaglutide vs Semaglutide 0.48% (0.17,0.8)
injectable_semaglutide_hba1c_percent <- 0.48
injectable_semaglutide_hba1c_mmol <- injectable_semaglutide_hba1c_percent * 10.929



### Make predictions ----
analysis_standard <- analysis_standard %>%
  mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide_adjustment)

analysis_oral <- analysis_oral %>%
  mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide_adjustment)

analysis_injectable <- analysis_injectable %>%
  mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide_adjustment)



# Unified validation ########################################

## Standard ----
unified_analysis_standard_adj <- unified_validation(
  data = analysis_standard, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


## Oral semaglutide ----
unified_analysis_standard_oral_adj <- unified_validation(
  data = analysis_standard %>% filter(drugclass != "GLP1") %>% rbind(analysis_oral), 
  drug_var = "drugclass",
  drugs = c("GLP1", "SGLT2", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("GLP1", "SGLT2", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


## Injectable semaglutide ----
unified_analysis_standard_inje_adj <- unified_validation(
  data = analysis_standard %>% filter(drugclass != "GLP1") %>% rbind(analysis_injectable %>% filter(split == "testing") %>% select(-split)), 
  drug_var = "drugclass",
  drugs = c("GLP1", "DPP4", "SGLT2", "SU", "TZD"),
  prediction_vars = paste0("pred.", c("Inje", "DPP4", "SGLT2", "SU", "TZD")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)



# Plots ############################################

drug_names <- c(
  "SGLT2" = "SGLT2i",
  "GLP1"= "GLP-1RA",
  "SU" = "SU",
  "DPP4" = "DPP4i",
  "TZD" = "TZD"
)

## Standard ----
plot_standard_drugpairs <- unified_analysis_standard_adj %>%
  mutate(
    drug1 = str_replace_all(drug1, drug_names),
    drug2 = str_replace_all(drug2, drug_names)
  ) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")


pdf("Outputs/CPRD/06.drug_pair_standard.pdf", width = 12, height = 5)
plot_standard_drugpairs +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 13)
  )
dev.off()

## Oral semaglutide ----

drug_names <- c(
  "SGLT2" = "SGLT2i",
  "GLP1"= "Oral Semaglutide",
  "SU" = "SU",
  "DPP4" = "DPP4i",
  "TZD" = "TZD"
)

plot_oral_drugpairs <- unified_analysis_standard_oral_adj %>%
  filter(drug1 == "GLP1") %>%
  mutate(
    drug1 = str_replace_all(drug1, drug_names),
    drug2 = str_replace_all(drug2, drug_names)
  ) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")

pdf("Outputs/CPRD/06.drug_pair_oral_semaglutide.pdf", width = 6, height = 6)
plot_oral_drugpairs +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 13)
  )
dev.off()


## Injectable semaglutide ----

drug_names <- c(
  "SGLT2" = "SGLT2i",
  "GLP1"= "Injectable Semaglutide",
  "SU" = "SU",
  "DPP4" = "DPP4i",
  "TZD" = "TZD"
)

plot_inje_drugpairs <- unified_analysis_standard_inje_adj %>%
  filter(drug1 == "GLP1") %>%
  mutate(
    drug1 = str_replace_all(drug1, drug_names),
    drug2 = str_replace_all(drug2, drug_names)
  ) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")


pdf("Outputs/CPRD/06.drug_pair_injectable_semaglutide.pdf", width = 6, height = 6)
plot_inje_drugpairs +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 12)
  )
dev.off()




















