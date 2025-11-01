
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

# All baseline covariates to test
all_covariates <- c("age_t0", "gender_t0", "ideology_t0", "pid_t0",
                    "pol_interest_t0", "climate_t0", "religion_t0",
                    "immigration_t0", "healthcare_t0", "abortion_t0")

# Function to test leave-one-out for a specific outcome module
test_leave_one_out_module <- function(left_out_var, data, module_name) {

  # Get list of all continuous covariates (exclude the left-out one)
  continuous_covs <- c("age_t0", "gender_t0", "ideology_t0", "pid_t0",
                       "pol_interest_t0", "climate_t0", "religion_t0",
                       "immigration_t0")
  continuous_remaining <- setdiff(continuous_covs, left_out_var)

  # Get list of categorical covariates (exclude the left-out one if applicable)
  categorical_covs <- c("healthcare_t0", "abortion_t0")
  categorical_remaining <- setdiff(categorical_covs, left_out_var)

  # Prepare data with dummies for probability forest
  # Use ignore_na = TRUE to handle missing categorical values
  if (length(categorical_remaining) > 0) {
    data_dummies <- dummy_cols(data,
                               select_columns = categorical_remaining,
                               remove_first_dummy = TRUE,
                               remove_selected_columns = FALSE,
                               ignore_na = TRUE)
  } else {
    data_dummies <- data
  }

  # Filter to only rows where left-out variable and treatment are non-missing
  data_analysis <- data_dummies %>%
    filter(!is.na(.data[[left_out_var]]) & !is.na(treat_ind))

  # Get all t0 variables for forest (excluding left-out)
  if (left_out_var %in% categorical_covs) {
    # If left-out is categorical, exclude all its dummies
    t0_vars_forest <- c(continuous_remaining,
                       names(data_analysis)[grepl(paste0("^(",
                                                        paste(categorical_remaining, collapse="|"),
                                                        ")_"), names(data_analysis))])
  } else {
    # If left-out is continuous, just exclude it
    t0_vars_forest <- c(continuous_remaining,
                       names(data_analysis)[grepl(paste0("^(",
                                                        paste(categorical_covs, collapse="|"),
                                                        ")_"), names(data_analysis))])
  }

  # Calculate propensity scores using GRF (excluding left-out variable)
  forest_probs <- probability_forest(Y = as.factor(data_analysis$treat_ind),
                                     X = data_analysis[, t0_vars_forest])

  balwts <- 1 / forest_probs$predictions[cbind(1:nrow(data_analysis),
                                               as.numeric(as.factor(data_analysis$treat_ind)))]

  # Test for "effect" on left-out variable using simple weighted regression
  # (not Lin-adjusted - cleaner test of whether reweighting balanced this variable)
  mod_placebo <- lm_robust(reformulate('treat_ind', left_out_var),
                          weights = balwts,
                          data = data_analysis)

  # Extract results
  coef_est <- coef(mod_placebo)["treat_indTRUE"]
  se <- mod_placebo$std.error["treat_indTRUE"]
  pval <- mod_placebo$p.value["treat_indTRUE"]

  # Add significance stars
  sig_stars <- ifelse(pval < 0.001, "***",
                     ifelse(pval < 0.01, "**",
                            ifelse(pval < 0.05, "*",
                                   ifelse(pval < 0.10, "+", ""))))

  return(data.frame(
    Module = module_name,
    Variable = left_out_var,
    Estimate = coef_est,
    SE = se,
    P_value = pval,
    Sig = sig_stars,
    N = nrow(data_analysis)
  ))
}

# Define outcome modules and their completion filters
modules <- list(
  primary = df_analysis %>% filter(finished_dv_primary == TRUE),
  secondary = df_analysis %>% filter(finished_dv_sec == TRUE),
  therm = df_analysis %>% filter(finished_dv_therm_trans == TRUE)
)

# Run leave-one-out tests for each module and each covariate
loo_results_all <- data.frame()

for (module_name in names(modules)) {
  cat(sprintf("Testing %s outcome module (N=%d)...\n",
              module_name, nrow(modules[[module_name]])))

  module_data <- modules[[module_name]] %>%
    select(all_of(c("treat_ind", all_covariates)))

  for (var in all_covariates) {
    result <- tryCatch(
      test_leave_one_out_module(var, module_data, module_name),
      error = function(e) {
        cat(sprintf("  Error testing %s: %s\n", var, e$message))
        data.frame(Module = module_name, Variable = var,
                  Estimate = NA, SE = NA, P_value = NA, Sig = "",
                  N = NA)
      }
    )
    loo_results_all <- rbind(loo_results_all, result)
  }
}

# Create formatted table
loo_latex <- loo_results_all %>%
  mutate(
    Diff_SE = sprintf("%.3f (%.3f)%s", Estimate, SE, Sig),
    P_value_fmt = sprintf("%.7f", P_value)
  ) %>%
  select(Module, Variable, Diff_SE, P_value_fmt, N)

# Write LaTeX table
sink("../tables/placebo_leave_one_out.tex")
cat("\\begin{table*}[!htbp]\n")
cat("\\centering\n")
cat("\\caption{Placebo: weighted treated--control differences in pre-treatment covariates (held-out from GRF)}\n")
cat("\\label{tab:placebo_weighting}\n")
cat("\\centering\n")
cat("\\begin{tabular}[t]{llccc}\n")
cat("\\toprule\n")
cat("Module & Variable & Diff (SE) & p-value & N\\\\\n")
cat("\\midrule\n")

current_module <- NULL
for (i in 1:nrow(loo_latex)) {
  # Add spacing between modules
  if (!is.null(current_module) && loo_latex$Module[i] != current_module) {
    cat("\\addlinespace\n")
  }
  current_module <- loo_latex$Module[i]

  cat(sprintf("%s & %s & %s & %s & %d\\\\\n",
              loo_latex$Module[i],
              gsub("_", "\\\\_", loo_latex$Variable[i]),
              loo_latex$Diff_SE[i],
              loo_latex$P_value_fmt[i],
              loo_latex$N[i]))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\\\\\n")
cat("{\\raggedright\n")
cat("\\small\n")
cat("\n")
cat("    Note: $^{***} p<0.001$, $^{**} p<0.01$, $^* p<0.05$, $^+ p<0.1$. Included are respondents for whom we collected complete post-test response. \n")
cat("    \n")
cat("}\n")
cat("\\end{table*}\n")
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

