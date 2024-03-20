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
df_analysis <- readRDS("../data/df_for_analysis_processed.rds")

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
  ))

data <- data %>%
  mutate(finished_dv_sec = ifelse(
    rowSums(is.na(data[, 10:11])) == 0 & rowSums(data[, 10:11] == 999) == 0,
    TRUE,
    FALSE
  ))

data <- data %>%
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
  # post-treatment response measures
  ## Dependent variable finished
  # Create a bootstrap sample
  sampled_data <- data[indices, ]
  
  ############# transgender tolerance scale #############

  # We provide simple difference in means estimates, with robust standard errors,
  mod1_tlr <- lm_robust(tolerance_index ~ treat_ind,
                        data = sampled_data)
  
  # And we estimate treatment effects using the Lin estimator
  lin_cov <- formula(
    paste('~', 
          paste(t0.covariate.names, collapse = ' + '),
          ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))
  
  mod2_tlr <- lm_lin(
    tolerance_index ~ treat_ind,
    covariates = lin_cov,
    data = sampled_data)
  
  # # estimation with weighted missing value
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  # # Convert to dummy variables and bind them to the original dataframe
  data_partial <- sampled_data
  data_partial$attrited <- 0
  data_partial$attrited[data_partial$finished_dv_primary == FALSE] <- 1
  data_partial$treat_factor <- as.factor(ifelse(data_partial$attrited == 0,
                                                data_partial$treat_ind,
                                                2)+1)
  data_partial_dummies <- dummy_cols(data_partial, select_columns = 
                                       c("healthcare_t0", 
                                         "abortion_t0"),
                                     remove_first_dummy = TRUE,
                                     remove_selected_columns = TRUE)
  forest_probs <- probability_forest(
    Y = data_partial_dummies$treat_factor,
    X = data_partial_dummies[,grep('t0', names(data_partial_dummies))])
  balwts <- 1/forest_probs$predictions[
    cbind(1:nrow(data_partial), as.numeric(data_partial$treat_factor))
  ]
  
  balwts <- balwts[data_partial$attrited == 0]
  data_complete <- data_partial[data_partial_dummies$attrited == 0,]
  
  mod3_tlr <- lm_lin(
    tolerance_index ~ treat_ind,
    covariates = lin_cov,
    data = data_complete,
    weights = balwts)
  
  
  # Store the models in a list
  estimate_results_tlr <- c(coef(mod1_tlr)[["treat_indTRUE"]], 
                     coef(mod2_tlr)[["treat_indTRUE"]], 
                     coef(mod3_tlr)[["treat_indTRUE"]])
  
  
  ############# transgender policy acceptance #############
  
  mod1_law <- lm_robust(laws_index ~ treat_ind,
                        data = sampled_data)
  
  # And we estimate treatment effects using the Lin estimator
  mod2_law <- lm_lin(
    laws_index ~ treat_ind,
    covariates = lin_cov,
    data = sampled_data)
  
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  
  data_partial <- sampled_data
  data_partial$attrited <- 0
  data_partial$attrited[data_partial$finished_dv_sec == FALSE] <- 1
  data_partial$treat_factor <- as.factor(ifelse(data_partial$attrited == 0,
                                                data_partial$treat_ind,
                                                2)+1)
  data_partial_dummies <- dummy_cols(data_partial, select_columns = 
                                       c("healthcare_t0", 
                                         "abortion_t0"),
                                     remove_first_dummy = TRUE,
                                     remove_selected_columns = TRUE)
  forest_probs <- probability_forest(
    Y = data_partial_dummies$treat_factor,
    X = data_partial_dummies[,grep('t0', names(data_partial_dummies))])
  balwts <- 1/forest_probs$predictions[
    cbind(1:nrow(data_partial), as.numeric(data_partial$treat_factor))
  ]
  
  balwts <- balwts[data_partial$attrited == 0]
  data_complete <- data_partial[data_partial$attrited == 0,]
  
  mod3_law <- lm_lin(
    laws_index ~ treat_ind,
    covariates = lin_cov,
    data = data_complete,
    weights = balwts)
  
  estimate_results_law = c(coef(mod1_law)[["treat_indTRUE"]], 
                       coef(mod2_law)[["treat_indTRUE"]], 
                       coef(mod3_law)[["treat_indTRUE"]])
  
  
  ############# feelings towards transgender #############
  
  # We provide simple difference in means estimates, with robust standard errors,
  mod1_therm <- lm_robust(therm_trans_t1 ~ treat_ind,
                          data = sampled_data)
  
  # And we estimate treatment effects using the Lin estimator
  lin_cov <- formula(
    paste('~', 
          paste(t0.covariate.names, collapse = ' + '),
          ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))
  
  mod2_therm <- lm_lin(
    therm_trans_t1 ~ treat_ind,
    covariates = lin_cov,
    data = sampled_data)
  
  # # estimation with weighted missing value
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  # # Convert to dummy variables and bind them to the original dataframe
  data_partial <- sampled_data
  data_partial$attrited <- 0
  data_partial$attrited[data_partial$finished_dv_therm_thrans == FALSE] <- 1
  data_partial$treat_factor <- as.factor(ifelse(data_partial$attrited == 0,
                                                data_partial$treat_ind,
                                                2)+1)
  data_partial_dummies <- dummy_cols(data_partial, select_columns = 
                                       c("healthcare_t0", 
                                         "abortion_t0"),
                                     remove_first_dummy = TRUE,
                                     remove_selected_columns = TRUE)
  forest_probs <- probability_forest(
    Y = data_partial_dummies$treat_factor,
    X = data_partial_dummies[,grep('t0', names(data_partial_dummies))])
  balwts <- 1/forest_probs$predictions[
    cbind(1:nrow(data_partial), as.numeric(data_partial$treat_factor))
  ]
  
  balwts <- balwts[data_partial$attrited == 0]
  data_complete <- data_partial[data_partial_dummies$attrited == 0,]
  
  mod3_therm <- lm_lin(
    therm_trans_t1 ~ treat_ind,
    covariates = lin_cov,
    data = data_complete,
    weights = balwts)
  
  estimate_results_therm <- c(coef(mod1_therm)[["treat_indTRUE"]], 
                     coef(mod2_therm)[["treat_indTRUE"]], 
                     coef(mod3_therm)[["treat_indTRUE"]])
  
  estimate_results = c(estimate_results_tlr,
                         estimate_results_law,
                         estimate_results_therm)

  return(estimate_results)
  
}


# 
# ## getting results 
# results <- lapply(names(dataset_list), function(dataset_name) {
#   apply_estimations(dataset_list[[dataset_name]], dataset_name)
# })



############ bootstrapping for 95% CI #############
# 
# # # Bootstrap function for estimating the treatment effect
# bootstrap_estimate <- function(data, indices, formula) {
#   sampled_data = data[indices, ]  # Create a bootstrap sample
#   model = lm_robust(formula, data = sampled_data)
#   estimate = coef(model)[["treat_indTRUE"]]  # Adjust for your specific variable name
#   return(estimate)
# }

# bootstrap_estimates <- function(data, indices, model) {
#   sampled_data = data[indices, ]  # Create a bootstrap sample
#   model = model
#   estimate = coef(model)[["treat_indTRUE"]]  # Adjust for your specific variable name
#   return(estimate)
# }

#model = lm_robust(tolerance_index ~ treat_ind, data = df_lower_bound)
# # Assuming your treatment effect model formula
# formula = tolerance_index ~ treat_ind  # Adjust the formula as necessary
# # 
# df_upper_bound$tolerance_index[df_upper_bound$finished_dv_primary] <- compute.factor.dv(
#   all.dv.names.t1[-c(1,2)],
#   df_upper_bound,
#   respondent.booleans = df_upper_bound$finished_dv_primary)
# # 
# 
# df_lower_bound <- df_lower_bound %>%
#   mutate(finished_dv_primary = ifelse(
#     rowSums(is.na(df_lower_bound[, 2:9])) == 0 & rowSums(data[, 2:9] == 999) == 0,
#     TRUE,
#     FALSE
#   ))
# 
# df_lower_bound$tolerance_index[df_lower_bound$finished_dv_primary] <- compute.factor.dv(
#   all.dv.names.t1[-c(1,2)], 
#   df_lower_bound,
#   respondent.booleans = df_lower_bound$finished_dv_primary)
# 
# 
# # Bootstrapping for Upper Bound
# set.seed(123)  # Ensure reproducibility

boot_results_upper = boot(dataset_list$df_upper_bound, bootstrap_estimate, R = 10000)

ci_list <- vector("list", ncol(boot_results_upper$t))

# Calculate the BCa confidence intervals for each estimate
for (i in 1:ncol(boot_results_upper$t)) {
  ci_list[[i]] <- boot.ci(boot_results_upper, type = "bca", index = i)
}


# Bootstrapping for Lower Bound
set.seed(123)  # Ensure reproducibility
boot_results_lower = boot(dataset_list$df_lower_bound, bootstrap_estimate, R = 10000)

# 95% CI for the Upper Bound Estimates
ci_upper = boot.ci(boot_results_upper,type='bca')


# 95% CI for the Lower Bound Estimates
ci_lower = boot.ci(boot_results_lower, type = "bca")
ci_lower[["bca"]][4]
ci_upper[["bca"]][5]

