
# Load data, packages and functions
source('../code/utils.R')

# ------------------------------- Generate Balance Tables ------------------------

# Table 1: All respondents
cat("\n=== Balance Table: All Respondents ===\n")
create_balance_table(df_analysis, "../tables/balance_table_all.tex")

# Table 2: Completers only
cat("\n=== Balance Table: Completers Only ===\n")
df_completers <- df_analysis %>%
  filter(finished_dv_primary == TRUE & finished_dv_sec == TRUE)
create_balance_table(df_completers, "../tables/balance_table_completers.tex")


# ------------------------------- Placebo Analysis -------------------------------

## 1. Leave-one-out covariate placebo test
## Shows that reweighting procedure works correctly by demonstrating that when we
## reweight on all covariates EXCEPT one, there is no "treatment effect" on that
## left-out pre-treatment variable (since it's measured before treatment)

cat("\n=== Placebo Test 1: Leave-One-Out Covariate Test ===\n")
cat("Weighted treated-control differences in pre-treatment covariates (held-out from GRF)\n")
cat("Testing separately for each outcome module\n\n")

# IPW weights function from placebo.R (simpler implementation)
make_balancing_weights <- function(data, module = c("primary","secondary","therm"),
                                   drop_var = NULL) {
  module <- match.arg(module)
  data_partial <- data

  # module-specific attrition
  attr_flag <- switch(module,
                      "primary" = "finished_dv_primary",
                      "secondary" = "finished_dv_sec",
                      "therm" = "finished_dv_therm_trans")

  data_partial$attrited <- ifelse(data_partial[[attr_flag]], 0, 1)
  data_partial$treat_factor <- as.factor(ifelse(data_partial$attrited == 0,
                                                data_partial$treat_ind, 2) + 1)

  # dummies for healthcare_t0 & abortion_t0
  data_dummies <- fastDummies::dummy_cols(
    data_partial,
    select_columns = c("healthcare_t0","abortion_t0"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )

  # GRF features = all *_t0 columns, minus the held-out var
  x_cols <- grep('t0$', names(data_dummies), value = TRUE)
  if (!is.null(drop_var)) x_cols <- setdiff(x_cols, drop_var)

  # fit GRF and compute ipw
  forest_probs <- grf::probability_forest(
    Y = data_dummies$treat_factor,
    X = data_dummies[, x_cols, drop = FALSE]
  )
  balwts <- 1 / forest_probs$predictions[
    cbind(1:nrow(data_dummies), as.numeric(data_dummies$treat_factor))
  ]

  # keep non-attrition data
  keep <- data_dummies$attrited == 0
  list(weights = balwts[keep], data = data_partial[keep, , drop = FALSE])
}

# pretreatment placebo: leave-one-out balance on each baseline covariate
placebo_covariate_test <- function(holdout_var, module = "primary") {
  wdat <- make_balancing_weights(df_analysis, module = module, drop_var = holdout_var)
  f <- as.formula(paste0(holdout_var, " ~ treat_ind"))
  m <- estimatr::lm_robust(f, data = wdat$data, weights = wdat$weights, se_type = "HC2")

  est <- unname(coef(m)["treat_indTRUE"])
  se <- summary(m)$coef["treat_indTRUE","Std. Error"]
  p <- summary(m)$coef["treat_indTRUE","Pr(>|t|)"]
  n <- stats::nobs(m)
  data.frame(module = module, variable = holdout_var, est = est, se = se, p = p, N = n)
}

# run leave-one-out test for all baseline covariates
all_covs <- c(t0.covariate.names, "healthcare_t0","abortion_t0")

results_pre_primary <- do.call(rbind, lapply(all_covs, placebo_covariate_test, module = "primary"))
results_pre_secondary <- do.call(rbind, lapply(all_covs, placebo_covariate_test, module = "secondary"))
results_pre_therm <- do.call(rbind, lapply(all_covs, placebo_covariate_test, module = "therm"))

results_pre <- dplyr::bind_rows(results_pre_primary, results_pre_secondary, results_pre_therm) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "+", TRUE ~ ""
  ))

# LaTeX table for leave-one-out test (formatted like balance table)
results_pre_formatted <- results_pre %>%
  dplyr::mutate(
    p_formatted = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)),
    Covariate = recode(variable,
                       "age_t0" = "Age",
                       "gender_t0" = "Gender (Male)",
                       "ideology_t0" = "Ideology",
                       "pid_t0" = "Party ID",
                       "pol_interest_t0" = "Political Interest",
                       "healthcare_t0" = "Healthcare Importance",
                       "climate_t0" = "Climate Concern",
                       "religion_t0" = "Religiosity",
                       "abortion_t0" = "Abortion Importance",
                       "immigration_t0" = "Immigration Attitudes")
  )

# Write as tabular environment only
sink("../tables/placebo_leave_one_out.tex")
cat("\\begin{tabular}{llccc}\n")
cat("\\toprule\n")
cat("Module & Covariate & Estimate & p-value & N \\\\\n")
cat("\\midrule\n")

current_module <- ""
for (i in 1:nrow(results_pre_formatted)) {
  # Add spacing between modules
  if (results_pre_formatted$module[i] != current_module && current_module != "") {
    cat("\\addlinespace\n")
  }
  current_module <- results_pre_formatted$module[i]
  
  # Estimate row with significance stars
  cat(sprintf("%s & %s & %.3f%s & %s & %d \\\\\n",
              results_pre_formatted$module[i],
              results_pre_formatted$Covariate[i],
              results_pre_formatted$est[i],
              results_pre_formatted$sig[i],
              results_pre_formatted$p_formatted[i],
              results_pre_formatted$N[i]))
  
  # Standard error row
  cat(sprintf(" & & (%.3f) & & \\\\\n", results_pre_formatted$se[i]))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
sink()


## 2. Placebo outcomes test
## Shows that treatment doesn't affect unrelated post-treatment policy questions
## Using the full reweighting procedure (Model 3)


# Compute weighted estimates for placebo outcomes
compute_weighted_estimate_placebo <- function(data, outcome_var) {
  # Standard Model 3 reweighting procedure
  data$attrited <- is.na(data[[outcome_var]])
  data$treat_factor <- as.factor(ifelse(data$attrited == 0, data$treat_ind, 2) + 1)

  data_dummies <- dummy_cols(data, select_columns = c("healthcare_t0", "abortion_t0"),
                             remove_first_dummy = TRUE,
                             remove_selected_columns = TRUE)

  forest_probs <- probability_forest(Y = data$treat_factor,
                                     X = data_dummies[, grep('t0', names(data_dummies))])

  balwts <- 1 / forest_probs$predictions[cbind(1:nrow(data), as.numeric(data$treat_factor))]

  # Estimate on complete cases
  data_complete <- data[!data$attrited, ]
  balwts_complete <- balwts[!data$attrited]

  lm_lin(reformulate('treat_ind', outcome_var),
         covariates = lin_cov,
         weights = balwts_complete,
         data = data_complete)
}

# Test placebo outcomes
placebo_outcomes <- c("tax_policy_t1", "marijuana_policy_t1", "min_wage_policy_t1")

placebo_results <- data.frame(
  Outcome = character(),
  Model = character(),
  Estimate = numeric(),
  SE = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (outcome in placebo_outcomes) {
  # Model 1: Difference-in-means
  mod1 <- lm_robust(reformulate('treat_ind', outcome), data = df_analysis)

  # Model 2: Covariate adjusted
  mod2 <- lm_lin(reformulate('treat_ind', outcome),
                covariates = lin_cov,
                data = df_analysis)

  # Model 3: Covariate adjusted + reweighted
  mod3 <- tryCatch(
    compute_weighted_estimate_placebo(df_analysis, outcome),
    error = function(e) NULL
  )

  # Extract results
  placebo_results <- rbind(placebo_results,
    data.frame(
      Outcome = outcome,
      Model = "Difference-in-means",
      Estimate = round(coef(mod1)["treat_indTRUE"], 4),
      SE = round(mod1$std.error["treat_indTRUE"], 4),
      P_value = round(mod1$p.value["treat_indTRUE"], 4)
    ),
    data.frame(
      Outcome = outcome,
      Model = "Covariate adjusted",
      Estimate = round(coef(mod2)["treat_indTRUE"], 4),
      SE = round(mod2$std.error["treat_indTRUE"], 4),
      P_value = round(mod2$p.value["treat_indTRUE"], 4)
    )
  )

  if (!is.null(mod3)) {
    placebo_results <- rbind(placebo_results,
      data.frame(
        Outcome = outcome,
        Model = "Cov adjusted + reweighted",
        Estimate = round(coef(mod3)["treat_indTRUE"], 4),
        SE = round(mod3$std.error["treat_indTRUE"], 4),
        P_value = round(mod3$p.value["treat_indTRUE"], 4)
      )
    )
  }
}


# Generate LaTeX table for placebo outcomes
sink("../tables/placebo_outcomes.tex")
cat("\\begin{tabular}{lcccc}\n")
cat("\\toprule\n")
cat("Outcome & Model & Estimate & SE & P-value \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(placebo_results)) {
  if (i > 1 && placebo_results$Outcome[i] != placebo_results$Outcome[i-1]) {
    cat("\\midrule\n")
  }

  outcome_label <- switch(placebo_results$Outcome[i],
    "tax_policy_t1" = "Tax policy",
    "marijuana_policy_t1" = "Marijuana policy",
    "min_wage_policy_t1" = "Minimum wage policy"
  )

  cat(sprintf("%s & %s & %.4f & %.4f & %.4f \\\\\n",
              outcome_label,
              placebo_results$Model[i],
              placebo_results$Estimate[i],
              placebo_results$SE[i],
              placebo_results$P_value[i]))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
sink()

