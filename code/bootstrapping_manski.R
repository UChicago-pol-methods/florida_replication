# Load data, packages and functions
source('./utils.R')

# ------------------------------- Data Preparation -------------------------------
dataset_list <- list(
  df_upper_bound = df_upper_bound,
  df_lower_bound = df_lower_bound
)

dataset_list <- lapply(dataset_list, preprocess_data)

# ------------------------------- Function Definitions -------------------------------
# function to calculate weighted estimates for manski bounds setting
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
                                     X = data_dummies[,grep('t0', names(data_dummies))],
                                     num.trees = 500,
                                     num.threads = 1,
                                     seed = 60637)

  balwts <- 1/forest_probs$predictions[cbind(1:nrow(data), as.numeric(data$treat_factor))]

  lm_lin(reformulate('treat_ind', outcome_var), covariates = covariate_formula, weights = balwts, data = data)
}

# function to fit models and extract treatment effect estimates
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
# function to perform bootstrapping and estimate effects
bootstrap_estimate <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  estimate_results_tlr <- fit_models_and_extract_estimates("tolerance_index","primary", sampled_data)
  estimate_results_law <- fit_models_and_extract_estimates("laws_index", "sec", sampled_data)
  estimate_results_therm <- fit_models_and_extract_estimates("therm_trans_t1","therm_thrans", sampled_data)
  
  c(estimate_results_tlr, estimate_results_law, estimate_results_therm)
}

# bootstrapping
set.seed(60637)
boot_results_upper <- boot(dataset_list$df_upper_bound, bootstrap_estimate, R = 1e4)
saveRDS(boot_results_upper, "../data/boot_UB.rds")
#cat("Upper bound bootstrapping completed.\n")

set.seed(60637)
boot_results_lower <- boot(dataset_list$df_lower_bound, bootstrap_estimate, R = 1e4)
saveRDS(boot_results_lower, "../data/boot_LB.rds")
#cat("Lower bound bootstrapping completed.\n")

# ------------------------------- Confidence Interval Calculation -------------------------------
ci_upper_list <- list()
for (i in 1:9) {
  ci_upper_list[[i]] <- boot.ci(boot_results_upper, type = "bca", index = i, conf = 0.975)
}
#cat("Upper bound CIs completed.\n")

ci_lower_list <- list()
for (i in 1:9) {
  ci_lower_list[[i]] <- boot.ci(boot_results_lower, type = "bca", index = i, conf = 0.975)
}
#cat("Lower bound CIs completed.\n")

# ------------------------------- Save CI Results to LaTeX Table -------------------------------

# function to extract CI values and calculate p-values from boot.ci object
extract_ci_values <- function(ci_obj, boot_result, index) {
  ci_vals <- ci_obj$bca[4:5]  # 97.5% CI bounds
  bootstrap_stats <- boot_result$t[, index]
  estimate <- ci_obj$t0
  
  # Calculate p-value
  if (estimate > 0) {
    p_value <- mean(bootstrap_stats <= 0)
    custom_p <- p_value  
  } else {
    p_value <- mean(bootstrap_stats >= 0)
    custom_p <- p_value  
  }
  

  if (custom_p < 0.001) {
    sig_stars <- "***"
  } else if (custom_p < 0.01) {
    sig_stars <- "**"
  } else if ((estimate > 0 & custom_p < 0.04) | (estimate < 0 & custom_p < 0.02)) {
    sig_stars <- "*"
  } else if (custom_p < 0.1) {
    sig_stars <- "+"
  } else {
    sig_stars <- ""
  }
  
  data.frame(
    Estimate = estimate,
    CI_Lower = ci_vals[1],
    CI_Upper = ci_vals[2],
    P_value = custom_p,
    Significance = sig_stars,
    stringsAsFactors = FALSE
  )
}

# Extract CI results for upper bound 
ci_upper_df <- do.call(rbind, lapply(1:length(ci_upper_list), function(i) {
  outcome_names <- c("Transgender tolerance scale", "Transgender policy acceptance", "Feelings towards transgender")
  model_names <- c("Difference-in-means", "Covariate adjusted", "Covariate adjusted and re-weighted")
  outcome_idx <- ((i - 1) %/% 3) + 1
  model_idx <- ((i - 1) %% 3) + 1
  
  result <- extract_ci_values(ci_upper_list[[i]], boot_results_upper, i)
  result$Outcome <- outcome_names[outcome_idx]
  result$Method <- model_names[model_idx]
  result
}))

# Extract CI results for lower bound 
ci_lower_df <- do.call(rbind, lapply(1:length(ci_lower_list), function(i) {
  outcome_names <- c("Transgender tolerance scale", "Transgender policy acceptance", "Feelings towards transgender")
  model_names <- c("Difference-in-means", "Covariate adjusted", "Covariate adjusted and re-weighted")
  outcome_idx <- ((i - 1) %/% 3) + 1
  model_idx <- ((i - 1) %% 3) + 1
  
  result <- extract_ci_values(ci_lower_list[[i]], boot_results_lower, i)
  result$Outcome <- outcome_names[outcome_idx]
  result$Method <- model_names[model_idx]
  result
}))

# Combine upper and lower bounds
combined_df <- merge(ci_lower_df, ci_upper_df, 
                     by = c("Outcome", "Method"),
                     suffixes = c("_LB", "_UB")) %>%
  mutate(
    Lower_Bound_Formatted = paste0(sprintf("%.4f", Estimate_LB), Significance_LB),
    Upper_Bound_Formatted = paste0(sprintf("%.4f", Estimate_UB), Significance_UB),
    CI_Formatted = paste0("(", sprintf("%.4f", CI_Lower_LB), ", ", sprintf("%.4f", CI_Upper_UB), ")")
  ) %>%
  select(Outcome, Method, Lower_Bound_Formatted, Upper_Bound_Formatted, CI_Formatted) %>%
  mutate(Outcome_Order = case_when(
    Outcome == "Feelings towards transgender" ~ 1,
    Outcome == "Transgender tolerance scale" ~ 2,
    Outcome == "Transgender policy acceptance" ~ 3
  )) %>%
  arrange(Outcome_Order,
          match(Method, c("Difference-in-means", 
                          "Covariate adjusted", 
                          "Covariate adjusted and re-weighted"))) %>%
  select(-Outcome_Order)

# Create LaTeX table
latex_table <- combined_df %>%
  select(-Outcome) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = c('l', 'c', 'c', 'c'),
        col.names = c("Method", "Lower Bound", "Upper Bound", "95% CI"),
        linesep = "") %>%
  pack_rows("Feelings towards transgender", 1, 3, bold = FALSE, italic = TRUE, indent = FALSE) %>%
  pack_rows("Transgender tolerance scale", 4, 6, bold = FALSE, italic = TRUE, indent = FALSE) %>%
  pack_rows("Transgender policy acceptance", 7, 9, bold = FALSE, italic = TRUE, indent = FALSE)

# Save to file
writeLines(latex_table, "../tables/manski_bounds_ci.tex")
cat("CI results saved to tables/manski_bounds_ci.tex\n")
