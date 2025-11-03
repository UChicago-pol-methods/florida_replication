# code/placebo_tests.R
# Placebo checks for the reweighting/attrition procedure

source('../code/utils.R')


# IPW weights
# drop one baseline covariate from the GRF feature set.
make_balancing_weights <- function(data, module = c("primary","secondary","therm"),
                                   drop_var = NULL) {
  module <- match.arg(module)
  data_partial <- data

  # module-specific attrition
  attr_flag <- switch(module,
                      "primary" = "finished_dv_primary",
                      "secondary" = "finished_dv_sec",
                      "therm" = "finished_dv_therm_thrans")

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

# 1) pretreatment placebo: leave-one-out balance on each baseline covariate
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

# LaTeX table
results_pre %>%
  dplyr::mutate(
    est_se = paste0(sprintf("%.3f", est), " (", sprintf("%.3f", se), ")", sig)
  ) %>%
  dplyr::select(Module = module, Variable = variable, `Diff (SE)` = est_se, `p-value` = p, N) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE, align = c("l","l","c","c","c"),
                    caption = "Placebo: weighted treated–control differences in pre-treatment covariates (held-out from GRF)") %>%
  kableExtra::save_kable("../tables/placebo_pre_treatment_balance.tex")

