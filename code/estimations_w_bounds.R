# Load data, packages and functions
source('../code/utils.R')

# ------------------------------- Data Preparation -------------------------------
dataset_list <- list(
  df_analysis = df_analysis,
  df_upper_bound = df_upper_bound, 
  df_lower_bound = df_lower_bound
)

# Preprocess data
dataset_list <- lapply(dataset_list, preprocess_data)

# ------------------------------- Function to apply statistical estimations on the datasets  -------------------------------
apply_estimations <- function(data, dataset_name) {

  # ------------------------------- Analysis: Tolerance Index -------------------------------
  # We provide simple difference in means estimates, with robust standard errors,
  
  mod1_tlr <- lm_robust(tolerance_index ~ treat_ind,
                        data = data)
  
  # And we estimate treatment effects using the Lin estimator
  mod2_tlr <- lm_lin(
    tolerance_index ~ treat_ind,
    covariates = lin_cov,
    data = data)
  
  # # estimation with weighted missing value
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  # # Convert to dummy variables and bind them to the original dataframe
  data_partial <- data
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
  
  
  # Store the models in a list
  models_tlr <- list(mod1 = mod1_tlr, mod2 = mod2_tlr, mod3 = mod3_tlr)
  
  
  # ------------------------------- Analysis: Laws Index -------------------------------
  
  # We provide simple difference in means estimates, with robust standard errors,

  mod1_law <- lm_robust(laws_index ~ treat_ind,
                        data = data)
  
  # And we estimate treatment effects using the Lin estimator
  mod2_law <- lm_lin(
    laws_index ~ treat_ind,
    covariates = lin_cov,
    data = data)
  
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  
  data_partial <- data
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
  
  models_law <- list(mod1 = mod1_law, mod2 = mod2_law, mod3 = mod3_law)
  

  # ------------------------------- Analysis: Therm Trans -------------------------------
  
  # We provide simple difference in means estimates, with robust standard errors,
  mod1_therm <- lm_robust(therm_trans_t1 ~ treat_ind,
                          data = data)
  
  # And we estimate treatment effects using the Lin estimator
  
  mod2_therm <- lm_lin(
    therm_trans_t1 ~ treat_ind,
    covariates = lin_cov,
    data = data)
  
  # # estimation with weighted missing value
  # Conditional on completing the pre-test response module and having been 
  # assigned treatment, we estimate propensity to be in each condition AND to not 
  # have attrited. We use inverse probability weights as sample weights in the Lin 
  # estimator for this sample. 
  # # Convert to dummy variables and bind them to the original dataframe
  data_partial <- data
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
  
  models_therm <- list(mod1 = mod1_therm, mod2 = mod2_therm, mod3 = mod3_therm)
  
  # ------------------------------- Generate Custom Summaries -------------------------------
  custom_summary_tlr_l <- lapply(models_tlr, create_custom_summary)
  custom_summary_law_l <- lapply(models_law, create_custom_summary)
  custom_summary_therm_l <- lapply(models_therm, create_custom_summary)
  
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
    
    # save the LaTeX table code
    filename<-paste0("../tables/",dataset_name,"_",name,".tex")
    writeLines(latex_table,filename)
    }
}

# ------------------------------- Apply Estimations to Each Dataset -------------------------------
results <- lapply(names(dataset_list), function(dataset_name) {
  apply_estimations(dataset_list[[dataset_name]], dataset_name)
})


