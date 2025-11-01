# --------------------------------------- Load Packages ---------------------------------------
library(boot)
library(data.table)
library(dplyr)
library(estimatr)
library(fastDummies)
library(grf)
library(kableExtra)
library(modelsummary)
library(tibble)
library(tidyr)
set.seed(60637)

# --------------------------------------- Read Data ---------------------------------------
df_analysis <- readRDS("../data/df_for_analysis_processed.rds")

# --------------------------------------- Define Variables ---------------------------------------
# Post-treatment response measures
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

# Baseline covariates
t0.covariate.names <- grep("t0$", names(df_analysis), value=TRUE)
t0.covariate.names <- t0.covariate.names[!t0.covariate.names %in%
                                           c('healthcare_t0', 'abortion_t0')]

# Lin estimator covariates formula
lin_cov <- formula(
  paste('~',
        paste(t0.covariate.names, collapse = ' + '),
        ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))

# --------------------------------------- Helper Functions ---------------------------------------

# Replace NA with min/max observed values (for Manski bounds)
replace_na_with_min <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  return(ifelse(is.na(x), min_val, x))
}

replace_na_with_max <- function(x) {
  max_val <- max(x, na.rm = TRUE)
  return(ifelse(is.na(x), max_val, x))
}

# Compute factor analysis outcome
compute.factor.dv <- function(dv.names,
                              data,
                              respondent.booleans,
                              print.loadings = TRUE){
  responders <- data[respondent.booleans,]
  # Factor analysis
  factor.obj <- princomp(responders[, dv.names], cor=TRUE)
  if(print.loadings) print(loadings(factor.obj))
  dv <- as.vector(factor.obj$scores[,1])
  # Rescale to mean 0 sd 1 in placebo group
  dv <- (dv - mean(dv[!data[,'treat_ind']], na.rm=TRUE)) /
    sd(dv[!data[,'treat_ind']], na.rm=TRUE)
  return(as.vector(dv))
}

# Formatting helper
f1 <- function(x) format(round(x, 3), big.mark=",")

# Custom p-value and significance level (Kalla & Broockman approach)
create_custom_summary <- function(model) {
  # Extract coefficients, standard errors, and p-values
  coef_df <- as.data.frame(coef(summary(model)))
  names(coef_df) <- c("Estimate", "StdError", "t value", "Pr(>|t|)", "CI Lower", "CI Upper", "DF")

  # Apply custom significance criteria: one-sided p-values for positive effects
  coef_df$Custom_P_Value <- ifelse(coef_df$Estimate > 0, coef_df$`Pr(>|t|)` / 2, coef_df$`Pr(>|t|)`)
  coef_df$Significance <- ifelse(coef_df$Custom_P_Value < 0.001, "***",
                                 ifelse(coef_df$Custom_P_Value < 0.01, "**",
                                        ifelse((coef_df$Estimate > 0 & coef_df$Custom_P_Value < 0.04) |
                                               (coef_df$Estimate < 0 & coef_df$Custom_P_Value < 0.02), "*", "")))
  nobs_value <- model$nobs
  coef_df$Nobs <- nobs_value
  return(coef_df)
}

# Format model summary for table output
format_for_table <- function(model_summary, model_number) {
  # Create a data frame from the model summary
  model_data <- model_summary %>%
    as_tibble(rownames = "Term")

  # Filter for 'treat_indTRUE', '(Intercept)', and 'Nobs' terms
  term_data <- model_data %>%
    filter(Term %in% c("treat_indTRUE", "(Intercept)"))

  # Format the coefficients
  term_data <- term_data %>%
    mutate(
      Model = model_number,
      Formatted = paste0(
        round(Estimate, 4),
        Significance,
        " & (",
        "\\num{", round(StdError, 4), "}"
      )
    ) %>%
    select(Term, Model, Formatted)

  # Extract the Nobs value
  nobs_value <- model_data %>%
    select(Nobs) %>%
    slice(1)

  # Add the Nobs value as an additional row
  nobs_row <- tibble(Term = "Nobs", Model = model_number, Formatted = as.character(nobs_value))

  # Combine the term data with Nobs
  bind_rows(term_data, nobs_row)
}

# Extract results from model list for LaTeX table
extract_results <- function(model_list) {
  intercept_row <- model_list[[1]][rownames(model_list[[1]]) == "(Intercept)", ]
  treat_row <- model_list[[1]][rownames(model_list[[1]]) == "treat_indTRUE", ]

  intercept_adj <- model_list[[2]][rownames(model_list[[2]]) == "(Intercept)", ]
  treat_adj <- model_list[[2]][rownames(model_list[[2]]) == "treat_indTRUE", ]

  intercept_adj_rw <- model_list[[3]][rownames(model_list[[3]]) == "(Intercept)", ]
  treat_adj_rw <- model_list[[3]][rownames(model_list[[3]]) == "treat_indTRUE", ]

  N1 <- unique(model_list[[1]]$Nobs)
  N2 <- unique(model_list[[2]]$Nobs)
  N3 <- unique(model_list[[3]]$Nobs)

  # Create a data frame with proper LaTeX formatting
  df <- data.frame(
    Variable = c("Intercept", "", "Treatment", ""),
    DiffMeans = c(
      paste0("\\num{", sprintf("%.4f", intercept_row$Estimate), "}"),
      paste0("(", "\\num{", sprintf("%.4f", intercept_row$StdError), "})"),
      paste0("\\num{", sprintf("%.4f", treat_row$Estimate), "}", treat_row$Significance),
      paste0("(", "\\num{", sprintf("%.4f", treat_row$StdError), "})")
    ),
    CovAdj = c(
      paste0("\\num{", sprintf("%.4f", intercept_adj$Estimate), "}"),
      paste0("(", "\\num{", sprintf("%.4f", intercept_adj$StdError), "})"),
      paste0("\\num{", sprintf("%.4f", treat_adj$Estimate), "}", treat_adj$Significance),
      paste0("(", "\\num{", sprintf("%.4f", treat_adj$StdError), "})")
    ),
    CovAdjRW = c(
      paste0("\\num{", sprintf("%.4f", intercept_adj_rw$Estimate), "}"),
      paste0("(", "\\num{", sprintf("%.4f", intercept_adj_rw$StdError), "})"),
      paste0("\\num{", sprintf("%.4f", treat_adj_rw$Estimate), "}", treat_adj_rw$Significance),
      paste0("(", "\\num{", sprintf("%.4f", treat_adj_rw$StdError), "})")
    )
  )

  # Add N row separately with LaTeX formatting
  df <- rbind(df, data.frame(
    Variable = "N",
    DiffMeans = paste0("\\num{", sprintf("%.0f", N1), "}"),
    CovAdj = paste0("\\num{", sprintf("%.0f", N2), "}"),
    CovAdjRW = paste0("\\num{", sprintf("%.0f", N3), "}")
  ))

  return(df)
}

# Generate LaTeX table from data frame
generate_latex_table <- function(df, title) {
  cat("\\midrule\n\\textit{", title, "} \\\\\n", sep = "")
  for (i in seq_len(nrow(df))) {
    if (df$Variable[i] == "") {
      cat(" & ", df$DiffMeans[i], " & ", df$CovAdj[i], " & ", df$CovAdjRW[i],
          " \\\\\n", sep = "")
    } else {
      cat(df$Variable[i], " & ", df$DiffMeans[i], " & ", df$CovAdj[i], " & ",
          df$CovAdjRW[i], " \\\\\n", sep = "")
    }
  }
}

# Save LaTeX table to file
save_latex_table <- function(filename, therm_results, tlr_results, law_results) {
  sink(filename)  # Redirect output to file

  cat("\\begin{tabular}{lccc}\n")
  cat("\\toprule\n")
  cat(" & Difference-in-means & Covariate adjusted & Covariate adjusted and re-weighted \\\\\n")

  generate_latex_table(therm_results, "Feelings towards transgender individuals")
  generate_latex_table(tlr_results, "Transgender tolerance scale")
  generate_latex_table(law_results, "Transgender policy acceptance")

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")

  sink()  # Stop redirecting output
}

# Preprocess data for dependent variables
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

  # Tolerance scale
  data$tolerance_index[data$finished_dv_primary] <- compute.factor.dv(
    all.dv.names.t1[-c(1,2)],
    data,
    respondent.booleans = data$finished_dv_primary)

  # Policy acceptance
  data$laws_index <- rowMeans(data[, all.dv.names.t1[1:2]])

  return(data)
}


# --------------------------------------- Manski Bounds Data ---------------------------------------
# Starting Lower Bound Calculation
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

# Starting Upper Bound Calculation
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

# ------------------------------- Balance Table Function -------------------------------

create_balance_table <- function(data, filename) {
  # Get list of all t0 covariates
  all_covariates <- c(t0.covariate.names, 'healthcare_t0', 'abortion_t0')
  
  # Calculate balance statistics by treatment
  balance_stats <- data %>%
    group_by(treat_ind) %>%
    summarise(across(all_of(all_covariates),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                     .names = "{.col}_{.fn}"))
  
  # Extract control and treatment stats
  control_stats <- balance_stats %>% filter(treat_ind == FALSE)
  treatment_stats <- balance_stats %>% filter(treat_ind == TRUE)
  
  # Create mean and SE rows for the table
  balance_table <- data.frame(Covariate = character(), stringsAsFactors = FALSE)
  
  for (var in all_covariates) {
    # Get variable label
    var_label <- recode(var,
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
    
    # Get means and SEs
    control_mean <- control_stats[[paste0(var, "_mean")]]
    control_se <- control_stats[[paste0(var, "_se")]]
    treatment_mean <- treatment_stats[[paste0(var, "_mean")]]
    treatment_se <- treatment_stats[[paste0(var, "_se")]]
    
    # Calculate difference
    diff <- treatment_mean - control_mean
    
    # Run t-test for statistical significance
    control_vals <- data[[var]][data$treat_ind == FALSE]
    treatment_vals <- data[[var]][data$treat_ind == TRUE]
    
    if (var %in% c('healthcare_t0', 'abortion_t0')) {
      # For categorical variables, use chi-square test
      tbl <- table(data$treat_ind, data[[var]])
      test_result <- tryCatch({
        chisq.test(tbl)
      }, error = function(e) {
        list(p.value = NA)
      })
      p_val <- test_result$p.value
      # For categorical, SE of difference is approximated
      diff_se <- sqrt(control_se^2 + treatment_se^2)
    } else {
      # For continuous variables, use t-test
      test_result <- t.test(treatment_vals, control_vals)
      p_val <- test_result$p.value
      # SE of difference from t-test
      diff_se <- test_result$stderr
    }
    
    # Get significance stars
    sig_stars <- ifelse(p_val < 0.001, "***",
                        ifelse(p_val < 0.01, "**",
                               ifelse(p_val < 0.05, "*",
                                      ifelse(p_val < 0.1, "+", ""))))
    
    # Create mean row
    row_mean <- c(var_label,
                  sprintf("%.3f", control_mean),
                  sprintf("%.3f", treatment_mean),
                  paste0(sprintf("%.3f", diff), sig_stars),
                  sprintf("%.3f", p_val))
    
    # Create SE row
    row_se <- c("",
                sprintf("(%.3f)", control_se),
                sprintf("(%.3f)", treatment_se),
                sprintf("(%.3f)", diff_se),
                "")
    
    balance_table <- rbind(balance_table, row_mean, row_se)
  }
  
  # Add N row
  n_control <- sum(data$treat_ind == FALSE)
  n_treatment <- sum(data$treat_ind == TRUE)
  n_row <- c("N", sprintf("%.0f", n_control), sprintf("%.0f", n_treatment), "", "")
  balance_table <- rbind(balance_table, n_row)
  
  # Set column names
  colnames(balance_table) <- c("Covariate", "Control", "Treatment", "Difference", "P_value")
  
  # Create LaTeX table (just the tabular)
  latex_balance_table <- balance_table %>%
    kable(format = "latex", booktabs = TRUE, escape = FALSE,
          align = c('l', 'c', 'c', 'c', 'c'),
          linesep = "",
          col.names = c("Covariate", "Control", "Treatment", "Difference", "P-value")) %>%
    add_header_above(c(" " = 1, "Mean" = 2, " " = 1, " " = 1))
  
  # Save LaTeX table
  writeLines(latex_balance_table, filename)
  
  # Print the LaTeX table
  print(latex_balance_table)
}
