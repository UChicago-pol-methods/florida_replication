# Load data, packages and functions
source('./code/utils.R')

# ------------------------------- Data Preparation -------------------------------
dataset_list <- list(
  df_upper_bound = df_upper_bound,
  df_lower_bound = df_lower_bound
)

dataset_list <- lapply(dataset_list, preprocess_data)

# ------------------------------- Function Definitions -------------------------------
# Function to calculate weighted estimates for manski bounds setting
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

# Function to fit models and extract treatment effect estimates
fit_models_and_extract_estimates <- function(outcome_var, attrition, data) {

  mod1 <- lm_robust(reformulate('treat_ind', outcome_var), data = data)
  mod2 <- lm_lin(reformulate('treat_ind', outcome_var), covariates = lin_cov, data = data)

  mod3 <- compute_weighted_estimate(data, outcome_var, lin_cov)

  c(
    coef(mod1)[["treat_indTRUE"]], coef(mod2)[["treat_indTRUE"]],
    coef(mod3)[["treat_indTRUE"]]
  )
}

# ------------------------------- Bootstrapping Analysis -------------------------------

# Function to perform bootstrapping and estimate effects
bootstrap_estimate <- function(data, indices) {
  # Create a bootstrap sample
  sampled_data <- data[indices, ]

  # Apply the function for each outcome variable
  estimate_results_tlr <- fit_models_and_extract_estimates("tolerance_index","primary", sampled_data)
  estimate_results_law <- fit_models_and_extract_estimates("laws_index", "sec", sampled_data)
  estimate_results_therm <- fit_models_and_extract_estimates("therm_trans_t1","therm_thrans", sampled_data)

  # Combine all estimates
  c(estimate_results_tlr, estimate_results_law, estimate_results_therm)
}


# Bootstrapping (with reduced number of iterations for testing)
set.seed(60637)
boot_results_upper <- boot(dataset_list$df_upper_bound, bootstrap_estimate, R = 1e4, parallel="multicore")
saveRDS(boot_results_upper, "../data/boot_UB.rds")

set.seed(60637)
boot_results_lower <- boot(dataset_list$df_lower_bound, bootstrap_estimate, R = 1e4, parallel="multicore")
saveRDS(boot_results_lower, "../data/boot_LB.rds")


# ------------------------------- Confidence Interval Calculation -------------------------------


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

