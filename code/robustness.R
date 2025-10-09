
# Load data, packages and functions
source('../code/utils.R')

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
    sig_stars <- get_significance_stars(sig_stars)

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

# ------------------------------- Generate Tables -------------------------------

# Table 1: All respondents
cat("\n=== Balance Table: All Respondents ===\n")
create_balance_table(df_analysis, "../tables/balance_table_all.tex")

# Table 2: Completers only
cat("\n=== Balance Table: Completers Only ===\n")
df_completers <- df_analysis %>%
  filter(finished_dv_primary == TRUE & finished_dv_sec == TRUE)
create_balance_table(df_completers, "../tables/balance_table_completers.tex")
