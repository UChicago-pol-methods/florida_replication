library(estimatr)
library(fastDummies)
library(grf)
library(modelsummary)
library(dplyr)
library(tidyr)
library(tibble)
library(kableExtra)
set.seed(60637)

## read in data 
df_analysis <- readRDS("~/Desktop/Uchicago/molly.nosync/DD-notes-memos/chatbot-replication/data/df_for_analysis_processed.rds")

# post-treatment response measures
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

## Actual Lower/Upper Bound Calculation

replace_na_with_min <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  return(ifelse(is.na(x), min_val, x))
}

# Function to replace NA with the maximum observed value
replace_na_with_max <- function(x) {
  max_val <- max(x, na.rm = TRUE)
  return(ifelse(is.na(x), max_val, x))
}

# Creating a lower bound dataframe
df_lower_bound_actual<- df_analysis
for (var in all.dv.names.t1) {
  df_lower_bound_actual[[var]] <- replace_na_with_min(df_analysis[[var]])
}

for (var in trans.law.dvs.t1) {
  df_lower_bound_actual[[var]] <- replace_na_with_min(df_analysis[[var]])
}

# Creating an upper bound dataframe
df_upper_bound_actual <- df_analysis
for (var in all.dv.names.t1) {
  df_upper_bound_actual[[var]] <- replace_na_with_max(df_analysis[[var]])
}

for (var in trans.law.dvs.t1) {
  df_upper_bound_actual[[var]] <- replace_na_with_max(df_analysis[[var]])
}

## compute factor analysis outcome
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

## Dependent variable
df_analysis$tolerance_index[df_analysis$finished_dv_primary] <- compute.factor.dv(
  all.dv.names.t1[-c(1,2)], 
  df_analysis,
  respondent.booleans = df_analysis$finished_dv_primary)
df_analysis$laws_index <- rowMeans(df_analysis[, all.dv.names.t1[1:2]])

## ---------------------------------- tolerance_index---------------------------------##
# We provide simple difference in means estimates, with robust standard errors,
mod1_tlr <- lm_robust(tolerance_index ~ treat_ind,
                      data = df_analysis)

# And we estimate treatment effects using the Lin estimator
lin_cov <- formula(
  paste('~', 
        paste(t0.covariate.names, collapse = ' + '),
        ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))

mod2_tlr <- lm_lin(
  tolerance_index ~ treat_ind,
  covariates = lin_cov,
  data = df_analysis)

# # estimation with weighted missing value
# Conditional on completing the pre-test response module and having been 
# assigned treatment, we estimate propensity to be in each condition AND to not 
# have attrited. We use inverse probability weights as sample weights in the Lin 
# estimator for this sample. 
# # Convert to dummy variables and bind them to the original dataframe
data_partial <- df_analysis
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

## ---------------------------------- laws_index---------------------------------##
# We provide simple difference in means estimates, with robust standard errors,
mod1_law <- lm_robust(laws_index ~ treat_ind,
                      data = df_analysis)

# And we estimate treatment effects using the Lin estimator
mod2_law <- lm_lin(
  laws_index ~ treat_ind,
  covariates = lin_cov,
  data = df_analysis)

# Conditional on completing the pre-test response module and having been 
# assigned treatment, we estimate propensity to be in each condition AND to not 
# have attrited. We use inverse probability weights as sample weights in the Lin 
# estimator for this sample. 

data_partial <- df_analysis
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


## ---------------------------------- therm_trans_t1---------------------------------##
# We provide simple difference in means estimates, with robust standard errors,
mod1_therm <- lm_robust(therm_trans_t1 ~ treat_ind,
                        data = df_analysis)

# And we estimate treatment effects using the Lin estimator
lin_cov <- formula(
  paste('~', 
        paste(t0.covariate.names, collapse = ' + '),
        ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))

mod2_therm <- lm_lin(
  therm_trans_t1 ~ treat_ind,
  covariates = lin_cov,
  data = df_analysis)

# # estimation with weighted missing value
# Conditional on completing the pre-test response module and having been 
# assigned treatment, we estimate propensity to be in each condition AND to not 
# have attrited. We use inverse probability weights as sample weights in the Lin 
# estimator for this sample. 
# # Convert to dummy variables and bind them to the original dataframe
data_partial <- df_analysis
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


# Tables

tlr_l <- list(mod1_tlr,
              mod2_tlr,
              mod3_tlr)

law_l <- list(mod1_law,
              mod2_law,
              mod3_law)

therm_l <- list(mod1_therm,
                mod2_therm,
                mod3_therm)

f1 <- function(x) format(round(x, 3), big.mark=",")

## customize p-value and significant level 
create_custom_summary <- function(model) {
  # Extract coefficients, standard errors, and p-values
  coef_df <- as.data.frame(coef(summary(model)))
  names(coef_df) <- c("Estimate", "StdError", "t value", "Pr(>|t|)", "CI Lower", "CI Upper", "DF")
  
  # Apply custom significance criteria and add significance stars
  coef_df$Custom_P_Value <- ifelse(coef_df$Estimate > 0, coef_df$`Pr(>|t|)` / 2, coef_df$`Pr(>|t|)`)
  coef_df$Significance <- ifelse(coef_df$Custom_P_Value < 0.001, "***",
                                 ifelse(coef_df$Custom_P_Value < 0.01, "**",
                                        ifelse((coef_df$Estimate > 0 & coef_df$Custom_P_Value < 0.04) | (coef_df$Estimate < 0 & coef_df$Custom_P_Value < 0.02), "*", "")))
  nobs_value <- model$nobs
  coef_df$Nobs<-nobs_value
  return(coef_df)
}

custom_summary_tlr_l <- lapply(tlr_l, create_custom_summary)
custom_summary_law_l <- lapply(law_l, create_custom_summary)
custom_summary_therm_l <- lapply(therm_l, create_custom_summary)


# Adjust the function to include Nobs
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
        get_significance_stars(Significance),
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

# Helper function to get significance stars
get_significance_stars <- function(sig) {
  ifelse(sig == "***", "***", 
         ifelse(sig == "**", "**", 
                ifelse(sig == "*", "*", "")))
}

# Apply the formatting function to each model summary
object_names <- c("custom_summary_therm_l","custom_summary_tlr_l", "custom_summary_law_l")

for (name in object_names) {
  model_list <- get(name)
  formatted_list <- lapply(seq_along(model_list), function(i) {
    format_for_table(model_list[[i]], i)
  })
  
  # Combine the formatted data from each model
  combined_formatted <- bind_rows(formatted_list)
  
  # Pivot the data to have models as columns and terms as rows
  final_table <- combined_formatted %>%
    pivot_wider(names_from = Model, values_from = Formatted, names_prefix = "Model_")
  
  final_table$Term <- recode(final_table$Term,
                             "(Intercept)" = "Intercept",
                             "treat_indTRUE" = "Treatment",
                             "Nobs" = "N")
  
  latex_table <- final_table %>%
    kable(format = "latex", booktabs = TRUE, escape = FALSE, align = c('l', 'c', 'c', 'c')) %>%
    add_header_above(c(" " = 1, "Difference-in-means" = 1, "Covariate adjusted" = 1, "Covariate adjusted and re-weighted" = 1))
  
  # Print the LaTeX table code
  print(latex_table)
}
# 
# modelsummary(custom_summary_tlr_l,
#              stars = TRUE,
#              escape = FALSE,
#              fmt = fmt_significant(digits=3),
#              coef_map = c('treat_indTRUE' = 'Treated',
#                           '(Intercept)' = 'Intercept'),
#              gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                  "fmt" = f1)))
# 
# out_tlr <- modelsummary(tlr_l,
#                         stars = TRUE,
#                         escape = FALSE,
#                         fmt = fmt_significant(digits=3),
#                         coef_map = c('treat_indTRUE' = 'Treated',
#                                      '(Intercept)' = 'Intercept'),
#                         gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                             "fmt" = f1)),
#                         output = 'latex_tabular')
# 
# writeLines(out_tlr, '../tables/tolerance_tables.tex')
# 
# modelsummary(law_l,
#              stars = TRUE,
#              escape = FALSE,
#              fmt = fmt_significant(digits=3),
#              coef_map = c('treat_indTRUE' = 'Treated',
#                           '(Intercept)' = 'Intercept'),
#              gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                  "fmt" = f1)))
# 
# out_law <- modelsummary(law_l,
#                         stars = TRUE,
#                         escape = FALSE,
#                         fmt = fmt_significant(digits=3),
#                         coef_map = c('treat_indTRUE' = 'Treated',
#                                      '(Intercept)' = 'Intercept'),
#                         gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                             "fmt" = f1)),
#                         output = 'latex_tabular')
# 
# 
# modelsummary(therm_l,
#              stars = TRUE,
#              escape = FALSE,
#              fmt = fmt_significant(digits=3),
#              coef_map = c('treat_indTRUE' = 'Treated',
#                           '(Intercept)' = 'Intercept'),
#              gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                  "fmt" = f1)))
# writeLines(out_therm, '../tables/law_tables.tex')
# 
# 
# out_therm <- modelsummary(therm_l,
#                         stars = TRUE,
#                         escape = FALSE,
#                         fmt = fmt_significant(digits=3),
#                         coef_map = c('treat_indTRUE' = 'Treated',
#                                      '(Intercept)' = 'Intercept'),
#                         gof_map = list(list("raw" = "nobs", "clean" = "N",
#                                             "fmt" = f1)),
#                         output = 'latex_tabular')
# 
# writeLines(out_therm, '../tables/therm_tables.tex')
# 
