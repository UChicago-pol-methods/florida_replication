# --------------------------------------- Load Packages ---------------------------------------
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(xtable)
library(estimatr)
library(fastDummies)
library(grf)
library(modelsummary)
library(tibble)
library(kableExtra)
library(boot)
set.seed(60637)

# --------------------------------------- Read Data ---------------------------------------
df_analysis <- readRDS("../data/df_for_analysis_processed.rds")

# --------------------------------------- Define Variables ---------------------------------------
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

## lin estimator covariates 
lin_cov <- formula(
  paste('~', 
        paste(t0.covariate.names, collapse = ' + '),
        ' + as.factor(healthcare_t0) + as.factor(abortion_t0)'))

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

f1 <- function(x) format(round(x, 3), big.mark=",")


## pre-processing data for dependent variables
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
  
  #### tolerance scale 
  data$tolerance_index[data$finished_dv_primary] <- compute.factor.dv(
    all.dv.names.t1[-c(1,2)], 
    data,
    respondent.booleans = data$finished_dv_primary)
  #### policy acceptance
  data$laws_index <- rowMeans(data[, all.dv.names.t1[1:2]])
  
  return(data)
}


# --------------------------------------- Manski Bounds Data Preparation ---------------------------------------
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


# ------------------------------- Results Compilation and Reporting Functions -------------------------------

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


# Helper function to get significance stars
get_significance_stars <- function(sig) {
  ifelse(sig == "***", "***",
         ifelse(sig == "**", "**",
                ifelse(sig == "*", "*",
                       ifelse(sig == "+", "+", ""))))
}


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

