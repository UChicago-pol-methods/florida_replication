
# Load data, packages and functions
source('./code/utils.R')

# ------------------------------- Generate Balance Tables ------------------------

# Table 1: All respondents
cat("\n=== Balance Table: All Respondents ===\n")
create_balance_table(df_analysis, "./tables/balance_table_all.tex")

# Table 2: Completers only
cat("\n=== Balance Table: Completers Only ===\n")
df_completers <- df_analysis %>%
  filter(finished_dv_primary == TRUE & finished_dv_sec == TRUE)
create_balance_table(df_completers, "./tables/balance_table_completers.tex")


# ------------------------------- Placebo Analysis -------------------------------

## 1. Leave-one-out covariate placebo test
## Shows that reweighting procedure works correctly by demonstrating that when we
## reweight on all covariates EXCEPT one, there is no "treatment effect" on that
## left-out pre-treatment variable (since it's measured before treatment)

cat("\n=== Placebo Test 1: Leave-One-Out Covariate Test ===\n")
cat("Weighted treated-control differences in pre-treatment covariates (held-out from GRF)\n")
cat("Testing separately for each outcome module\n\n")

# All baseline covariates to test (matches professor's list/order)

SEED <- 60637
set.seed(SEED)

all_covariates <- c("age_t0", "gender_t0", "ideology_t0", "pid_t0",
                    "pol_interest_t0", "climate_t0", "religion_t0",
                    "immigration_t0", "healthcare_t0", "abortion_t0")

normalize_module <- function(module) {
  allowed <- c("primary", "secondary", "therm")
  m <- as.character(module)[1]
  if (is.na(m) || !(m %in% allowed)) {
    stop(sprintf("`module` must be one of %s; got: %s", paste(allowed, collapse = ", "), deparse(module)))
  }
  m
}


make_balancing_weights <- function(data, module = "primary",
                                   drop_var = NULL, seed = SEED) {
  module <- normalize_module(module)
  dat <- data

  # Module-specific completion flag
  attr_flag <- switch(module,
                      "primary" = "finished_dv_primary",
                      "secondary" = "finished_dv_sec",
                      "therm" = "finished_dv_therm_trans")

  # Mark attrition
  dat$attrited <- ifelse(!is.na(dat[[attr_flag]]) & dat[[attr_flag]] == TRUE, 0L, 1L)

  # 3-class label: 1=control complete, 2=treated complete, 3=attrited (any arm)
  is_treated <- as.logical(dat$treat_ind)
  tf_base <- ifelse(dat$attrited == 0L, ifelse(is_treated, 1L, 0L), 2L)
  dat$treat_factor <- factor(tf_base + 1L, levels = c(1L, 2L, 3L))

  # Build dummies for categorical baseline vars
  has_ignore_na <- "ignore_na" %in% names(formals(fastDummies::dummy_cols))
  if (has_ignore_na) {
    dat_dum <- fastDummies::dummy_cols(
      dat,
      select_columns = c("healthcare_t0", "abortion_t0"),
      remove_first_dummy = TRUE,
      remove_selected_columns = FALSE,
      ignore_na = TRUE
    )
  } else {
    dat_dum <- fastDummies::dummy_cols(
      dat,
      select_columns = c("healthcare_t0", "abortion_t0"),
      remove_first_dummy = TRUE,
      remove_selected_columns = FALSE
    )
  }

  # Continuous vs. categorical covariates
  continuous_covs <- c("age_t0", "gender_t0", "ideology_t0", "pid_t0",
                       "pol_interest_t0", "climate_t0", "religion_t0",
                       "immigration_t0")
  categorical_covs <- c("healthcare_t0", "abortion_t0")

  # Exclude the held-out variable
  cont_keep <- setdiff(continuous_covs, drop_var)
  cat_keep  <- setdiff(categorical_covs, drop_var)

  # All dummies for remaining categorical covariates
  if (length(cat_keep) > 0) {
    pat <- paste0("^(", paste(cat_keep, collapse = "|"), ")_")
    cat_dummies <- names(dat_dum)[grepl(pat, names(dat_dum))]
  } else {
    cat_dummies <- character(0)
  }
  x_cols <- c(cont_keep, cat_dummies)
  if (length(x_cols) == 0) stop("No features left for GRF after exclusions.")

  # Rows available for forest: non-missing class label and continuous covariates observed
  keep_forest <- !is.na(dat$treat_factor)
  if (length(cont_keep) > 0) {
    # require complete cases on continuous-only covariates (dummies allow NA via ignore_na)
    keep_forest <- keep_forest & apply(dat_dum[, cont_keep, drop = FALSE], 1, function(r) all(!is.na(r)))
  }

  # Fit probability forest and get in-sample class probabilities
  set.seed(seed)
  pf <- grf::probability_forest(
    Y = dat$treat_factor[keep_forest],
    X = as.matrix(dat_dum[keep_forest, x_cols, drop = FALSE])
  )
  pred <- predict(pf)$predictions  # robust across GRF versions

  # Inverse probability weights for realized class
  realized_class <- as.integer(dat$treat_factor[keep_forest])
  w_sub <- 1 / pred[cbind(seq_len(nrow(pred)), realized_class)]

  balwts <- rep(NA_real_, nrow(dat))
  balwts[keep_forest] <- w_sub

  # Keep completed cases for this module (analysis sample)
  keep_analysis <- keep_forest & dat$attrited == 0L

  list(
    weights = balwts[keep_analysis],
    data = dat[keep_analysis, , drop = FALSE]
  )
}

# ------------------------------------------------------------------------------------------
# One leave-one-out test
# ------------------------------------------------------------------------------------------
test_leave_one_out_module <- function(left_out_var, module_name, seed = SEED) {
  wdat <- make_balancing_weights(df_analysis, module = module_name, drop_var = left_out_var, seed = seed)

  # Restrict to non-missing outcome var and treatment indicator
  nonmiss <- !is.na(wdat$data[[left_out_var]]) & !is.na(wdat$data$treat_ind)
  dat_use <- wdat$data[nonmiss, , drop = FALSE]
  w_use <- wdat$weights[nonmiss]

  if (nrow(dat_use) == 0) {
    return(data.frame(Module = module_name, Variable = left_out_var,
                      Estimate = NA_real_, SE = NA_real_, P_value = NA_real_,
                      Sig = "", N = 0L))
  }

  # Weighted regression (no Lin adjustment)
  f <- reformulate("treat_ind", left_out_var)
  mod <- estimatr::lm_robust(f, data = dat_use, weights = w_use)

  coef_names <- names(coef(mod))
  coef_name <- if ("treat_indTRUE" %in% coef_names) "treat_indTRUE" else "treat_ind"

  est <- unname(coef(mod)[coef_name])
  se <- mod$std.error[coef_name]
  p <- mod$p.value[coef_name]
  n <- stats::nobs(mod)

  sig <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.10) "+" else ""

  data.frame(Module = module_name, Variable = left_out_var,
             Estimate = est, SE = se, P_value = p, Sig = sig, N = n)
}

# ------------------------------------------------------------------------------------------
# Run tests across modules and covariates
# ------------------------------------------------------------------------------------------
test_leave_one_out_module <- function(left_out_var, module_name, seed = SEED) {
  wdat <- make_balancing_weights(df_analysis, module = module_name, drop_var = left_out_var, seed = seed)

  # Restrict to non-missing outcome var and treatment indicator
  nonmiss <- !is.na(wdat$data[[left_out_var]]) & !is.na(wdat$data$treat_ind)
  dat_use <- wdat$data[nonmiss, , drop = FALSE]
  w_use <- wdat$weights[nonmiss]

  if (nrow(dat_use) == 0) {
    return(data.frame(Module = module_name, Variable = left_out_var,
                      Estimate = NA_real_, SE = NA_real_, P_value = NA_real_,
                      Sig = "", N = 0L))
  }

  # Weighted regression (no Lin adjustment)
  f <- reformulate("treat_ind", left_out_var)
  mod <- estimatr::lm_robust(f, data = dat_use, weights = w_use)

  coef_names <- names(coef(mod))
  coef_name <- if ("treat_indTRUE" %in% coef_names) "treat_indTRUE" else "treat_ind"

  est <- unname(coef(mod)[coef_name])
  se <- mod$std.error[coef_name]
  p <- mod$p.value[coef_name]
  n <- stats::nobs(mod)

  sig <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.10) "+" else ""

  data.frame(Module = module_name, Variable = left_out_var,
             Estimate = est, SE = se, P_value = p, Sig = sig, N = n)
}

# ------------------------------------------------------------------------------------------
# Run tests across modules and covariates
# ------------------------------------------------------------------------------------------
modules <- list(
  primary = df_analysis %>% dplyr::filter(finished_dv_primary == TRUE),
  secondary = df_analysis %>% dplyr::filter(finished_dv_sec == TRUE),
  therm = df_analysis %>% dplyr::filter(finished_dv_therm_trans == TRUE)
)

loo_results_all <- data.frame()

for (module_name in names(modules)) {
  cat(sprintf("Testing %s outcome module (N=%d)...\n", module_name, nrow(modules[[module_name]])))
  for (var in all_covariates) {
    res <- tryCatch(
      test_leave_one_out_module(var, module_name, seed = SEED),
      error = function(e) {
        cat(sprintf("  Error testing %s in %s: %s\n", var, module_name, e$message))
        data.frame(Module = module_name, Variable = var,
                   Estimate = NA_real_, SE = NA_real_, P_value = NA_real_,
                   Sig = "", N = NA_integer_)
      }
    )
    loo_results_all <- dplyr::bind_rows(loo_results_all, res)
  }}

# Create formatted table
loo_latex <- loo_results_all %>%
    dplyr::mutate(
      Diff_SE     = sprintf("%.3f (%.3f)%s", Estimate, SE, Sig),
      P_value_fmt = sprintf("%.7f", P_value)
    ) %>%
    dplyr::select(Module, Variable, Diff_SE, P_value_fmt, N)

# Write LaTeX table
sink("./tables/placebo_leave_one_out.tex")
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
sink("./tables/placebo_outcomes.tex")
cat("\\begin{tabular}{lccc}\n")
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

