# Author: pcardoso
###############################################################################

# Overall benefit calibration

###############################################################################

# load 
library(tableone)
# remotes::install_github("harrelfe/rms")
library(rms)
library(qpdf)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval)
library(viridis)
library(patchwork)


# load function for generating dataset
source("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/SDRN/01.cohort_script.R")
source("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/Closedtest_functions/closedtest_continuous_function.R")
source("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/SDRN/03.impute_missingness.R")
source("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/SDRN/04.model_predictions.R")

# load dataset
analysis_cohort_raw <- set_up_data(data = "mm_20250506_t2d_1stinstance")

# load model
load("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/fivedrugmodel_5knot_share_20230823.Rdata")



# Pre-processing dataset ########################################

analysis_cohort <- analysis_cohort_raw %>%
		mutate(
				pated = paste(serialno, dstartdate, drug_class, sep = "."),
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

## All drug initiations ----
### By mice
analysis_cohort_mice_imputation <- imputation_methods(data = analysis_cohort,
		method = "mice", 
		mice.m = 20, 
		mice.maxit = 10,
		mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))




# Predictions from datasets #########################################

## All drug initiations ----
### Imputation by mice
analysis_cohort_prediction_mice <- predict_5drugmodel(analysis_cohort_mice_imputation,
		model = m1.5.final,
		drug_var = "drugclass",
		drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"))


### merge impute columns into main dataset
analysis_cohort <- analysis_cohort %>%
		select(-c("prebmi", "preegfr", "pretotalcholesterol", "prehdl", "prealt")) %>%
		cbind(
				analysis_cohort_prediction_mice %>%
						select(contains("mice_impute")) %>%
						rename_with(~ str_replace(., "_mice_impute", "")),
				analysis_cohort_prediction_mice %>%
						select(contains("pred.mice")) %>%
						rename_with(~ str_replace(., "pred\\.mice\\.", "pred."))
		)


# Split cohort into old vs new drugs #########################################

analysis_standard <- analysis_cohort %>%
		filter(!(drug_substance %in% c("Oral semaglutide", "Low-dose semaglutide", "Semaglutide, dose unclear", "High-dose semaglutide")))

analysis_standard_post_2019 <- analysis_cohort %>%
		filter(!(drug_substance %in% c("Oral semaglutide", "Low-dose semaglutide", "Semaglutide, dose unclear", "High-dose semaglutide"))) %>%
		separate(pated, into = c("serialno", "dstartdate", "drug_class"), sep = "\\.") %>%
		mutate(dstartdate = as.Date(dstartdate)) %>%
		filter(dstartdate >= "2019-01-01") %>%
		mutate(pated = paste(serialno, dstartdate, drug_class, sep = ".")) %>%
		select(all_of(colnames(analysis_standard)))

analysis_oral <- analysis_cohort %>%
		filter(drug_substance %in% c("Oral semaglutide"))

analysis_injectable <- analysis_cohort %>%
		filter(drug_substance %in% c("Low-dose semaglutide", "Semaglutide, dose unclear", "High-dose semaglutide"))



# Closed loop test #################################################

## Injectable Semaglutide ----

### mice variables
#### GLP1
closed_loop_test_injectable_semaglutide <- 7.33 # change this value to the value from CPRD

## Adjustment used: https://doi.org/10.1016/j.diabres.2021.108904
## Difference between Fig5 Dulaglutide vs Semaglutide 0.48% (0.17,0.8)
injectable_semaglutide_hba1c_percent <- 0.48
injectable_semaglutide_hba1c_mmol <- injectable_semaglutide_hba1c_percent * 10.929



### Make predictions ----
analysis_standard <- analysis_standard %>%
		mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide)

analysis_standard_post_2019 <- analysis_standard_post_2019 %>%
		mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide)

analysis_oral <- analysis_oral %>%
		mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide)

analysis_injectable <- analysis_injectable %>%
		mutate(pred.Inje = pred.GLP1 - closed_loop_test_injectable_semaglutide)


# Calibration slope ----
calibration_summary <- analysis_standard %>%
		rowwise() %>%
		mutate(pred = get(paste0("pred.", drugclass))) %>%
		ungroup() %>%
		group_by(drugclass) %>%
		do(tidy(lm(posthba1cfinal ~ pred, data = .), conf.int = TRUE)) %>%
		filter(term %in% c("(Intercept)", "pred")) %>%
		select(drugclass, term, estimate, conf.low, conf.high) %>%
		pivot_wider(
				names_from = term,
				values_from = c(estimate, conf.low, conf.high),
				names_glue = "{term}_{.value}"
		) %>%
		rbind(
				analysis_standard %>%
						rowwise() %>%
						mutate(pred = get(paste0("pred.", drugclass))) %>%
						ungroup() %>%
						do(tidy(lm(posthba1cfinal ~ pred, data = .), conf.int = TRUE)) %>%
						filter(term %in% c("(Intercept)", "pred")) %>%
						mutate(drugclass = "Overall") %>%
						select(drugclass, term, estimate, conf.low, conf.high) %>%
						pivot_wider(
								names_from = term,
								values_from = c(estimate, conf.low, conf.high),
								names_glue = "{term}_{.value}"
						)
		) %>%
		as.data.frame()


# Find optimal drug ----

standard_optimal <- analysis_standard %>%
		cateval::get_best_drugs(
				rank = 1,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
				final_var_name = "pred."
		) %>%
		select(matches("rank1|drugclass")) %>%
		mutate(concordant = drugclass == pred.rank1_drug_name)


standard_post_2019_optimal <- analysis_standard_post_2019 %>%
		cateval::get_best_drugs(
				rank = 1,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
				final_var_name = "pred."
		) %>%
		select(matches("rank1|drugclass")) %>%
		mutate(concordant = drugclass == pred.rank1_drug_name)

## Rank 1 ----

rank1_optimal_standard <- analysis_standard_post_2019 %>%
#		rbind(analysis_oral, analysis_injectable) %>%
		cateval::get_best_drugs(
				rank = 1,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
				final_var_name = "pred."
		) %>%
		select(matches("rank1|sex"))

rank1_optimal_standard_inje <- analysis_standard_post_2019 %>%
#		rbind(analysis_oral, analysis_injectable) %>%
		cateval::get_best_drugs(
				rank = 1,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje")),
				final_var_name = "pred."
		) %>%
		select(matches("rank1|sex"))

## Tolerance 3 ----
tolerance3_optimal_standard <- analysis_standard_post_2019 %>%
#		rbind(analysis_oral, analysis_injectable) %>%
		cateval::get_best_drugs(
				tolerance = 3,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
				final_var_name = "pred."
		) %>%
		select(matches("within|sex"))

tolerance3_optimal_standard_inje <- analysis_standard_post_2019 %>%
#		rbind(analysis_oral, analysis_injectable) %>%
		cateval::get_best_drugs(
				tolerance =3,
				column_names = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje")),
				final_var_name = "pred."
		) %>%
		select(matches("within|sex"))






# Overall calibration #########################################

## Rank 1 ----

### Post 2019-01-01 ----
### Exact matching only on best drug / sex / hba1c_10
overall_benefit_calibration_rank1_summary <- cateval::compute_overall_benefit_performance(
		data = analysis_standard %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var = "posthba1cfinal",
		pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

overall_benefit_calibration_rank1_one_group <- cateval::compute_overall_benefit(
		data = analysis_standard %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var = "posthba1cfinal",
		cal_groups = 1,
		pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

overall_benefit_calibration_rank1 <- cateval::compute_overall_benefit(
		data = analysis_standard %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var = "posthba1cfinal",
		cal_groups = 10,
		pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

overall_benefit_calibration_tolerance3_one_group <- cateval::compute_overall_benefit(
		data = analysis_standard %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var = "posthba1cfinal",
		cal_groups = 1,
		conc_tolerance = 3,
		pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

overall_benefit_calibration_tolerance3 <- cateval::compute_overall_benefit(
		data = analysis_standard %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var = "posthba1cfinal",
		cal_groups = 10,
		conc_tolerance = 3,
		pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)



# Plots ----


### Tolerance 3 optimal therapies ----
groups <- list(
		"1-drug combination" = 1,
		"2-drug combinations" = 2,
		"3-drug combinations" = 3,
		"4/5-drug combinations" = 4:5
)
drug_names <- c(
		"SGLT2" = "SGLT2i",
		"GLP1"= "GLP-1RA",
		"SU" = "SU",
		"DPP4" = "DPP4i",
		"TZD" = "TZD",
		"Inje" = "Injectable Semaglutide"
)



plot_tol3_optimal_standard <- tolerance3_optimal_standard %>%
		select(pred.within_3_of_best_drug_name) %>%
		mutate(
				pred.within_3_of_best_drug_name = str_replace_all(pred.within_3_of_best_drug_name, drug_names)
		) %>%
		unlist() %>%
		cateval::optimal_drug_comparison_plot(groups = groups, plot = TRUE)

plot_tol3_optimal_standard_values <- tolerance3_optimal_standard %>%
		select(pred.within_3_of_best_drug_name) %>%
		mutate(
				pred.within_3_of_best_drug_name = str_replace_all(pred.within_3_of_best_drug_name, drug_names)
		) %>%
		unlist() %>%
		cateval::optimal_drug_comparison_plot(groups = groups, plot = FALSE)


plot_tol3_optimal_standard_inje <- tolerance3_optimal_standard_inje %>%
		select(pred.within_3_of_best_drug_name) %>%
		mutate(
				pred.within_3_of_best_drug_name = str_replace_all(pred.within_3_of_best_drug_name, drug_names)
		) %>%
		unlist() %>%
		cateval::optimal_drug_comparison_plot(groups = groups, plot = TRUE)

plot_tol3_optimal_standard_inje_values <- tolerance3_optimal_standard_inje %>%
		select(pred.within_3_of_best_drug_name) %>%
		mutate(
				pred.within_3_of_best_drug_name = str_replace_all(pred.within_3_of_best_drug_name, drug_names)
		) %>%
		unlist() %>%
		cateval::optimal_drug_comparison_plot(groups = groups, plot = FALSE)




pdf("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/Outputs/SDRN/05.drug_combinations.pdf", width = 14, height = 6)
plot_tol3_optimal_standard
plot_tol3_optimal_standard_inje
dev.off()


### Rank 1 optimal therapy ----

# colours for drugs
drug_colours_viridis <- viridis::viridis(5)
drug_colours_viridis <- c(drug_colours_viridis, "grey")
names(drug_colours_viridis) <- c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje")

drug_colours_github <- c(
		"SGLT2" = "#E69F00",
		"GLP1"= "#56B4E9",
		"SU" = "#CC79A7",
		"DPP4" = "#0072B2",
		"TZD" = "#D55E00",
		"Inje" = "grey"
)

plot_best_drug_github <- patchwork::wrap_plots(
				
				# A: standard
				rank1_optimal_standard %>%
						rbind(
								rank1_optimal_standard %>%
										mutate(sex = "Overall")
						) %>%
						count(sex, pred.rank1_drug_name) %>%
						mutate(
								sex = factor(sex, levels = rev(c("Overall", "Male", "Female"))),
								pred.rank1_drug_name = factor(pred.rank1_drug_name, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje"))
						) %>%
						group_by(sex) %>%
						mutate(
								percent = n / sum(n) * 100
						) %>%
						ggplot(aes(x = sex, y = percent, fill = pred.rank1_drug_name)) +
						geom_col(width = 0.7) +
						coord_flip() +
						scale_fill_manual(
								values = drug_colours_github, 
								labels = drug_names
						) +
						theme_void() +
						theme(
								axis.text.y = element_text(),
								legend.title = element_blank(),
								legend.position = "top"
						),
				
				# B: standard + injectable
				rank1_optimal_standard_inje %>%
						rbind(
								rank1_optimal_standard_inje %>%
										mutate(sex = "Overall")
						) %>%
						count(sex, pred.rank1_drug_name) %>%
						mutate(
								sex = factor(sex, levels = rev(c("Overall", "Male", "Female"))),
								pred.rank1_drug_name = factor(pred.rank1_drug_name, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje"))
						) %>%
						complete(sex, pred.rank1_drug_name, fill = list(n = 0)) %>%
						group_by(sex) %>%
						mutate(
								percent = n / sum(n) * 100
						) %>%
						ggplot(aes(x = sex, y = percent, fill = pred.rank1_drug_name)) +
						geom_col(width = 0.7) +
						coord_flip() +
						scale_fill_manual(
								values = drug_colours_github, 
								labels = drug_names,
								drop = FALSE
						) +
						theme_void() +
						theme(
								axis.text.y = element_text(),
								legend.title = element_blank(),
								legend.position = "top"
						),
				
				ncol = 2, nrow = 1
		
		) &
		theme(
				legend.position = "top"
		)



plot_best_drug_viridis <- patchwork::wrap_plots(
				
				# A: standard
				rank1_optimal_standard %>%
						rbind(
								rank1_optimal_standard %>%
										mutate(sex = "Overall")
						) %>%
						count(sex, pred.rank1_drug_name) %>%
						mutate(
								sex = factor(sex, levels = rev(c("Overall", "Male", "Female"))),
								pred.rank1_drug_name = factor(pred.rank1_drug_name, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje"))
						) %>%
						group_by(sex) %>%
						mutate(
								percent = n / sum(n) * 100
						) %>%
						ggplot(aes(x = sex, y = percent, fill = pred.rank1_drug_name)) +
						geom_col(width = 0.7) +
						coord_flip() +
						scale_fill_manual(
								values = drug_colours_viridis, 
								labels = drug_names
						) +
						theme_void() +
						theme(
								axis.text.y = element_text(),
								legend.title = element_blank(),
								legend.position = "top"
						),
				
				# B: standard + injectable
				rank1_optimal_standard_inje %>%
						rbind(
								rank1_optimal_standard_inje %>%
										mutate(sex = "Overall")
						) %>%
						count(sex, pred.rank1_drug_name) %>%
						mutate(
								sex = factor(sex, levels = rev(c("Overall", "Male", "Female"))),
								pred.rank1_drug_name = factor(pred.rank1_drug_name, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD", "Inje"))
						) %>%
						complete(sex, pred.rank1_drug_name, fill = list(n = 0)) %>%
						group_by(sex) %>%
						mutate(
								percent = n / sum(n) * 100
						) %>%
						ggplot(aes(x = sex, y = percent, fill = pred.rank1_drug_name)) +
						geom_col(width = 0.7) +
						coord_flip() +
						scale_fill_manual(
								values = drug_colours_viridis, 
								labels = drug_names,
								drop = FALSE
						) +
						theme_void() +
						theme(
								axis.text.y = element_text(),
								legend.title = element_blank(),
								legend.position = "top"
						),
				
				ncol = 2, nrow = 1
		
		) &
		theme(
				legend.position = "top"
		)


pdf("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/Outputs/SDRN/05.optimal_therapy_rank1.pdf", width = 10, height = 4)
plot_best_drug_github
plot_best_drug_viridis
dev.off()


## Overall calibration ----

plot_calibration_rank1 <- overall_benefit_calibration_rank1 %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Rank 1 optimal")

plot_calibration_tol3 <- overall_benefit_calibration_tolerance3 %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "3mmol/mol optimal")

pdf("/home/pcardoso/workspace/Pedro-5drugmodel-Add_Drug_SDRN/Outputs/SDRN/05.post_overall_calibration.pdf", width = 7, height = 5)
plot_calibration_rank1
plot_calibration_tol3
dev.off()


