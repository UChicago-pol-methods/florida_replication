library(estimatr)
library(fastDummies)
library(grf)
library(modelsummary)
library(dplyr)
library(tidyr)
library(tibble)
library(kableExtra)
library(boot)
set.seed(60637)

## read in data 
df_analysis <- readRDS("~/Desktop/Uchicago/molly.nosync/DD-notes-memos/chatbot-replication/data/df_for_analysis_processed.rds")

all.dv.names.t1 <- c('florida_trans_policy_t1',
                     'florida_trans_policy2_t1',
                     'therm_trans_t1',
                     'gender_norm_sexchange_t1',
                     'gender_norm_moral_t1',
                     'gender_norm_abnormal_t1',
                     'gender_norm_trans_moral_wrong_t1',
                     'trans_teacher_t1',
                     'trans_bathroom_t1',
                     'gender_norm_dress_t1')
trans.law.dvs.t1 <- c('florida_trans_policy_t1', 'florida_trans_policy2_t1')

t0.covariate.names <- grep("t0$", names(df_analysis), value=TRUE)
t0.covariate.names <- t0.covariate.names[!t0.covariate.names %in% 
                                           c('healthcare_t0', 'abortion_t0')]
lin_cov <- formula(
  paste('~', 
        paste(t0.covariate.names, collapse = ' + '),
        ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))
### manski bounds data preparation #### 
# Possible Lower Bound Calculation
df_lower_bound <- df_analysis %>%
  mutate(
    therm_trans_t1 = ifelse(is.na(therm_trans_t1), 0, therm_trans_t1),
    gender_norm_sexchange_t1 = ifelse(is.na(gender_norm_sexchange_t1), -1, gender_norm_sexchange_t1),
    gender_norm_moral_t1 = ifelse(is.na(gender_norm_moral_t1), -1, gender_norm_moral_t1),
    gender_norm_abnormal_t1 = ifelse(is.na(gender_norm_abnormal_t1), -1, gender_norm_abnormal_t1),
    gender_norm_trans_moral_wrong_t1 = ifelse(is.na(gender_norm_trans_moral_wrong_t1), -1, gender_norm_trans_moral_wrong_t1),
    trans_teacher_t1 = ifelse(is.na(trans_teacher_t1), -1, trans_teacher_t1),
    trans_bathroom_t1 = ifelse(is.na(trans_bathroom_t1), -1, trans_bathroom_t1),
    gender_norm_dress_t1 = ifelse(is.na(gender_norm_dress_t1), -1, gender_norm_dress_t1),
    florida_trans_policy_t1 = ifelse(is.na(florida_trans_policy_t1), -3, florida_trans_policy_t1),
    florida_trans_policy2_t1 = ifelse(is.na(florida_trans_policy2_t1), -3, florida_trans_policy2_t1)
  )

# Possible Upper Bound Calculation
df_upper_bound <- df_analysis %>%
  mutate(
    therm_trans_t1 = ifelse(is.na(therm_trans_t1), 10, therm_trans_t1),
    gender_norm_sexchange_t1 = ifelse(is.na(gender_norm_sexchange_t1), 1, gender_norm_sexchange_t1),
    gender_norm_moral_t1 = ifelse(is.na(gender_norm_moral_t1), 1, gender_norm_moral_t1),
    gender_norm_abnormal_t1 = ifelse(is.na(gender_norm_abnormal_t1), 1, gender_norm_abnormal_t1),
    gender_norm_trans_moral_wrong_t1 = ifelse(is.na(gender_norm_trans_moral_wrong_t1), 1, gender_norm_trans_moral_wrong_t1),
    trans_teacher_t1 = ifelse(is.na(trans_teacher_t1), 1, trans_teacher_t1),
    trans_bathroom_t1 = ifelse(is.na(trans_bathroom_t1), 1, trans_bathroom_t1),
    gender_norm_dress_t1 = ifelse(is.na(gender_norm_dress_t1), 1, gender_norm_dress_t1),
    florida_trans_policy_t1 = ifelse(is.na(florida_trans_policy_t1), 3, florida_trans_policy_t1),
    florida_trans_policy2_t1 = ifelse(is.na(florida_trans_policy2_t1), 3, florida_trans_policy2_t1)
  )


compute.factor.dv <- function(dv.names,
                              data,
                              respondent.booleans,
                              print.loadings = TRUE){
  responders <- data[respondent.booleans,]
  # Factor analysis
  factor.obj <- princomp(responders[, dv.names], cor=TRUE)
  if(print.loadings) print(loadings(factor.obj))
  dv <- as.vector(factor.obj$scores[,1])
  # Rescale to mean 0 sd 1 in placebo group; treatment effects can then be
  # interpreted as the effect in standard deviations the treatment would have
  # among an untreated population.
  # if(cor(dv, responders$florida_trans_policy_t1) < 0) dv <- -1 * dv
  dv <- (dv - mean(dv[!data[,'treat_ind']],
                   na.rm=TRUE)) /
    sd(dv[!data[,'treat_ind']], na.rm=TRUE)
  return(as.vector(dv))
}


f1 <- function(x) format(round(x, 3), big.mark=",")

dataset_list <- list(
  df_upper_bound = df_upper_bound, 
  df_lower_bound = df_lower_bound
)


preprocess_data <- function(data) {
  data <- data %>%
    mutate(finished_dv_primary = ifelse(
      rowSums(is.na(data[, 2:9])) == 0 & rowSums(data[, 2:9] == 999) == 0,
      TRUE,
      FALSE
    )) %>%
    mutate(finished_dv_sec = ifelse(
      rowSums(is.na(data[, 10:11])) == 0 & rowSums(data[, 10:11] == 999) == 0,
      TRUE,
      FALSE
    )) %>%
    mutate(finished_dv_therm_thrans = ifelse(
      !is.na(data$therm_trans_t1) & data$therm_trans_t1 != 999,
      TRUE,
      FALSE
    ))
  
  #### tolerance scale 
  data$tolerance_index[data$finished_dv_primary] <- compute.factor.dv(
    all.dv.names.t1[-c(1,2)], 
    data,
    respondent.booleans = data$finished_dv_primary)
  #### policy acceptance
  data$laws_index <- rowMeans(data[, all.dv.names.t1[1:2]])
  
  return(data)
}

dataset_list <- lapply(dataset_list, preprocess_data)


# Function to apply bootstrapping
bootstrap_estimate <- function(data, indices) {
  # Create a bootstrap sample
  sampled_data <- data[indices, ]
  
  # Pre-compute covariates formula and dummy variables once, if they're constant across models
  covariate_formula <- formula(
    paste('~', 
          paste(t0.covariate.names, collapse = ' + '),
          ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))
  
  # Function to process data for attrition analysis
  # process_data_for_attrition <- function(data, finished_dv_col) {
  #   
  #   data$attrited <- 0 #in the manski bounds setting, there's no missing value
  #   data$treat_factor <- as.factor(ifelse(data$attrited == 0,
  #                                                 data$treat_ind,
  #                                                 2)+1)
  #   data_dummies <- dummy_cols(data, select_columns = 
  #                                c("healthcare_t0", 
  #                                  "abortion_t0"),
  #                              remove_first_dummy = TRUE,
  #                              remove_selected_columns = TRUE)
  #   return(list(data_dummies = data_dummies, data = data))
  #   
  # }

  # Function to calculate weighted estimates
  compute_weighted_estimate <- function(data, outcome_var, covariate_formula) {

    data$attrited <- 0 #in the manski bounds setting, there's no missing value
    data$treat_factor <- as.factor(ifelse(data$attrited == 0,
                                          data$treat_ind,
                                          2)+1)
    
    data_dummies <- dummy_cols(data, select_columns =
                                   c("healthcare_t0",
                                     "abortion_t0"),
                                 remove_first_dummy = TRUE,
                                 remove_selected_columns = TRUE)

    forest_probs <- probability_forest(Y = data$treat_factor,
                                       X = data_dummies[,grep('t0', names(data_dummies))])
    
    
    balwts <- 1/forest_probs$predictions[cbind(1:nrow(data), as.numeric(data$treat_factor))]
    

    # Filter based on `data` not `data_dummies` for the attrited and complete cases
    lm_lin(reformulate('treat_ind', outcome_var), covariates = covariate_formula, weights = balwts, data = data)
  }
  
  
  # Define a function to fit models and extract estimates
  fit_models_and_extract_estimates <- function(outcome_var, attrition, data) {
    mod1 <- lm_robust(reformulate('treat_ind', outcome_var), data = data)
    mod2 <- lm_lin(reformulate('treat_ind', outcome_var), covariates = covariate_formula, data = data)
    
    #data_processed <- process_data_for_attrition(data, sprintf("finished_dv_%s", attrition))

    #mod3 <- compute_weighted_estimate(data, outcome_var, covariate_formula)

    c(
      coef(mod1)[["treat_indTRUE"]], coef(mod2)[["treat_indTRUE"]]
      #coef(mod3)[["treat_indTRUE"]]
      )
  }
  
  # Apply the function for each outcome variable
  estimate_results_tlr <- fit_models_and_extract_estimates("tolerance_index","primary", sampled_data)
  estimate_results_law <- fit_models_and_extract_estimates("laws_index", "sec", sampled_data)
  estimate_results_therm <- fit_models_and_extract_estimates("therm_trans_t1","therm_thrans", sampled_data)
  
  # Combine all estimates
  c(estimate_results_tlr, estimate_results_law, estimate_results_therm)
}


# Bootstrapping (with reduced number of iterations for testing)

boot_results_upper <- boot(dataset_list$df_upper_bound, bootstrap_estimate, R = 10000)
boot_results_lower <- boot(dataset_list$df_lower_bound, bootstrap_estimate, R = 10000)

# 95% CI for the Upper Bound Estimates
# Initialize a list to store the CI results for each index
ci_upper_list = list()

# Loop through indices 1 to 6
for (i in 1:6) {
  ci_upper = boot.ci(boot_results_upper, type = 'bca', index = i)
  ci_upper_list[[i]] = ci_upper
}


# 95% CI for the Lower Bound Estimates
# Initialize a list to store the CI results for each index
ci_lower_list = list()

# Loop through indices 1 to 6
for (i in 1:6) {
  ci_lower = boot.ci(boot_results_lower, type = 'bca', index = i)
  ci_lower_list[[i]] = ci_lower
}

