
# Load data, packages and functions
source('../code/utils.R')

# ------------------------------- Data Preparation -------------------------------

# Computing the tolerance index based on certain variables
df_analysis$tolerance_index[df_analysis$finished_dv_primary] <- compute.factor.dv(
  all.dv.names.t1[-c(1,2)], 
  df_analysis,
  respondent.booleans = df_analysis$finished_dv_primary)

# Computing the laws index as the row mean of specified columns
df_analysis$laws_index <- rowMeans(df_analysis[, all.dv.names.t1[1:2]])

# ------------------------------- Analysis: Tolerance Index -------------------------------

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

#--------* laws_index---------------------------####
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


#--------* therm_trans_t1-----------------------####
# We provide simple difference in means estimates, with robust standard errors,
mod1_therm <- lm_robust(therm_trans_t1 ~ treat_ind,
                        data = df_analysis)

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
data_partial$attrited[data_partial$finished_dv_therm_trans == FALSE] <- 1
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


# Tables ####

tlr_l <- list(mod1_tlr,
              mod2_tlr,
              mod3_tlr)

law_l <- list(mod1_law,
              mod2_law,
              mod3_law)

therm_l <- list(mod1_therm,
                mod2_therm,
                mod3_therm)


custom_summary_tlr_l <- lapply(tlr_l, create_custom_summary)
custom_summary_law_l <- lapply(law_l, create_custom_summary)
custom_summary_therm_l <- lapply(therm_l, create_custom_summary)

therm_results <- extract_results(custom_summary_therm_l)
tlr_results <- extract_results(custom_summary_tlr_l)
law_results <- extract_results(custom_summary_law_l)


# Call function to save the table
save_latex_table("../tables/main_results.tex", therm_results, tlr_results, law_results)

# Save model objects for figures
save(mod1_tlr, mod2_tlr, mod3_tlr, mod1_law, mod2_law, mod3_law, mod1_therm, 
     mod2_therm, mod3_therm, file = "../data/models.RData")

#-------* Define groups-----------------------####

df_analysis_dems <- df_analysis |> 
  filter(pid_t0 %in% c(-3, -2, -1))

df_analysis_reps <- df_analysis |> 
  filter(pid_t0 %in% c(1, 2, 3))

df_analysis_inds <- df_analysis |> 
  filter(pid_t0 > -1 & pid_t0 < 1)

#-------* tolerance_index-----------------------####

## Democrats
mod1_tlrd <- lm_robust(tolerance_index ~ treat_ind, data = df_analysis_dems)
mod2_tlrd <- lm_lin(tolerance_index ~ treat_ind, covariates = lin_cov, data = df_analysis_dems)

## Republicans
mod1_tlrr <- lm_robust(tolerance_index ~ treat_ind, data = df_analysis_reps)
mod2_tlrr <- lm_lin(tolerance_index ~ treat_ind, covariates = lin_cov, data = df_analysis_reps)

## Independents
mod1_tlri <- lm_robust(tolerance_index ~ treat_ind, data = df_analysis_inds)
mod2_tlri <- lm_lin(tolerance_index ~ treat_ind, covariates = lin_cov, data = df_analysis_inds)

#-------* laws_index-----------------------####

## Democrats
mod1_lawd <- lm_robust(laws_index ~ treat_ind, data = df_analysis_dems)
mod2_lawd <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov, data = df_analysis_dems)

## Republicans
mod1_lawr <- lm_robust(laws_index ~ treat_ind, data = df_analysis_reps)
mod2_lawr <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov, data = df_analysis_reps)

## Independents
mod1_lawi <- lm_robust(laws_index ~ treat_ind, data = df_analysis_inds)
mod2_lawi <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov, data = df_analysis_inds)

#-------* therm_trans_t1-----------------------####

## Democrats
mod1_thermd <- lm_robust(therm_trans_t1 ~ treat_ind, data = df_analysis_dems)
mod2_thermd <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov, data = df_analysis_dems)

## Republicans
mod1_thermr <- lm_robust(therm_trans_t1 ~ treat_ind, data = df_analysis_reps)
mod2_thermr <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov, data = df_analysis_reps)

## Independents
mod1_thermi <- lm_robust(therm_trans_t1 ~ treat_ind, data = df_analysis_inds)
mod2_thermi <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov, data = df_analysis_inds)

#-------* Probability Weighting-----------------------####

apply_weighting <- function(data, dv_var) {
  data$attrited <- 0
  if (dv_var == "tolerance_index" | dv_var == "laws_index") {
    data$attrited[data$finished_dv_sec == FALSE] <- 1
  } else if (dv_var == "therm_trans_t1") {
    data$attrited[data$finished_dv_therm_trans == FALSE] <- 1
  }
  
  data$treat_factor <- as.factor(ifelse(data$attrited == 0, data$treat_ind, 2) + 1)
  
  data_dummies <- dummy_cols(data, select_columns = c("healthcare_t0", "abortion_t0"),
                             remove_first_dummy = TRUE, remove_selected_columns = TRUE)
  
  forest_probs <- probability_forest(
    Y = data_dummies$treat_factor,
    X = data_dummies[, grep('t0', names(data_dummies))]
  )
  
  balwts <- 1 / forest_probs$predictions[
    cbind(1:nrow(data), as.numeric(data$treat_factor))
  ]
  
  balwts <- balwts[data$attrited == 0]
  data_complete <- data[data$attrited == 0, ]
  
  return(list(data_complete = data_complete, weights = balwts))
}

# Apply weighting for each group and outcome
weighted_dems_tlr <- apply_weighting(df_analysis_dems, "tolerance_index")
mod3_tlrd <- lm_lin(tolerance_index ~ treat_ind, covariates = lin_cov,
                    data = weighted_dems_tlr$data_complete, weights = weighted_dems_tlr$weights)

weighted_reps_tlr <- apply_weighting(df_analysis_reps, "tolerance_index")
mod3_tlrr <- lm_lin(tolerance_index ~ treat_ind, covariates = ~age_t0 + gender_t0 + ideology_t0 + pid_t0 + pol_interest_t0 + 
                      climate_t0 + religion_t0 + immigration_t0 + as.factor(healthcare_t0) + 
                      as.factor(abortion_t0),
                    data = weighted_reps_tlr$data_complete, weights = weighted_reps_tlr$weights)

weighted_inds_tlr <- apply_weighting(df_analysis_inds, "tolerance_index")
mod3_tlri <- lm_lin(tolerance_index ~ treat_ind, covariates = lin_cov,
                    data = weighted_inds_tlr$data_complete, weights = weighted_inds_tlr$weights)

# Repeat for laws_index and therm_trans_t1
weighted_dems_law <- apply_weighting(df_analysis_dems, "laws_index")
mod3_lawd <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov,
                    data = weighted_dems_law$data_complete, weights = weighted_dems_law$weights)

weighted_reps_law <- apply_weighting(df_analysis_reps, "laws_index")
mod3_lawr <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov,
                    data = weighted_reps_law$data_complete, weights = weighted_reps_law$weights)

weighted_inds_law <- apply_weighting(df_analysis_inds, "laws_index")
mod3_lawi <- lm_lin(laws_index ~ treat_ind, covariates = lin_cov,
                    data = weighted_inds_law$data_complete, weights = weighted_inds_law$weights)

weighted_dems_therm <- apply_weighting(df_analysis_dems, "therm_trans_t1")
mod3_thermd <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov,
                      data = weighted_dems_therm$data_complete, weights = weighted_dems_therm$weights)

weighted_reps_therm <- apply_weighting(df_analysis_reps, "therm_trans_t1")
mod3_thermr <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov,
                      data = weighted_reps_therm$data_complete, weights = weighted_reps_therm$weights)

weighted_inds_therm <- apply_weighting(df_analysis_inds, "therm_trans_t1")
mod3_thermi <- lm_lin(therm_trans_t1 ~ treat_ind, covariates = lin_cov,
                      data = weighted_inds_therm$data_complete, weights = weighted_inds_therm$weights)

#------- Tables-----------------------####

save_latex_table("../tables/dem_results.tex", 
                 extract_results(lapply(list(mod1_thermd,
                                             mod2_thermd,
                                             mod3_thermd), 
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_tlrd,
                                             mod2_tlrd,
                                             mod3_tlrd),
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_lawd,
                                             mod2_lawd,
                                             mod3_lawd), 
                                        create_custom_summary)))

save_latex_table("../tables/rep_results.tex",
                 extract_results(lapply(list(mod1_thermr,
                                             mod2_thermr,
                                             mod3_thermr), 
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_tlrr,
                                             mod2_tlrr,
                                             mod3_tlrr),
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_lawr,
                                             mod2_lawr,
                                             mod3_lawr), 
                                        create_custom_summary)))

save_latex_table("../tables/ind_results.tex",
                 extract_results(lapply(list(mod1_thermi,
                                             mod2_thermi,
                                             mod3_thermi), 
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_tlri,
                                             mod2_tlri,
                                             mod3_tlri),
                                        create_custom_summary)),
                 extract_results(lapply(list(mod1_lawi,
                                             mod2_lawi,
                                             mod3_lawi), 
                                        create_custom_summary)))



#------- Model objects-----------------------####

save(
  mod1_tlrd, mod2_tlrd, mod3_tlrd, mod1_lawd, mod2_lawd, mod3_lawd, mod1_thermd, mod2_thermd, mod3_thermd,
  file = "../data/models_dems.RData"
)

save(
  mod1_tlrr, mod2_tlrr, mod3_tlrr, mod1_lawr, mod2_lawr, mod3_lawr, mod1_thermr, mod2_thermr, mod3_thermr,
  file = "../data/models_reps.RData"
)

save(
  mod1_tlri, mod2_tlri, mod3_tlri, mod1_lawi, mod2_lawi, mod3_lawi, mod1_thermi, mod2_thermi, mod3_thermi,
  file = "../data/models_inds.RData"
)


