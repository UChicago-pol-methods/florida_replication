library(estimatr)
library(grf)
set.seed(60637)

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


# simulate data
n <- 1e3 # simulated sample size
treat_ind <- sample(rep(0:1, n/2)) # treatment indicator


data <- replicate(length(all.dv.names.t1), rnorm(n))
colnames(data) <- all.dv.names.t1

data <- cbind(id = 1:n, data, treat_ind)

## covariates 
t0.covariate.names <- c('age_t0', 'gender_t0',
                        'ideology_t0', 'pid_t0', 'pol_interest_t0',
                        'healthcare_t0', 'climate_t0', 'religion_t0',
                        'abortion_t0', 'immigration_t0')

xmat <- replicate(length(t0.covariate.names), rnorm(n))
colnames(xmat) <- t0.covariate.names

data <- cbind(data, xmat)

## modified directly from Broockman & Kalla (2016) SI code
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
  dv <- (dv - mean(dv[!data[,'treat_ind']], 
                   na.rm=TRUE)) /
    sd(dv[!data[,'treat_ind']], na.rm=TRUE)
  
  return(as.vector(dv))
}

data <- as.data.frame(data)

data$tolerance_index <- compute.factor.dv(all.dv.names.t1[-c(1,2)], 
                                          data,
                                          respondent.booleans = rep(TRUE, n))

data$laws_index <- rowMeans(data[, all.dv.names.t1[1:2]])

# We provide simple difference in means estimates, with robust standard errors, 
lm_robust(tolerance_index ~ treat_ind,
          data = data)

# And we estimate treatment effects using the Lin estimator
lm_lin(tolerance_index ~ treat_ind,
       covariates = formula(paste('~', paste(t0.covariate.names, collapse = ' + '))),
       data = data)


# For missing response data
# Assume 10% post-treatment missingness
data_partial$attrited <- 0
data_partial$attrited[sample(1:n, n/10, replace = FALSE)] <- 1
data_partial[data_partial$attrited == 1,all.dv.names.t1] <- NA
data_partial$treat_factor <- as.factor(ifelse(data_partial$attrited == 0,
                                              data_partial$treat_ind,
                                              2)+1)

forest_probs <- probability_forest(Y = data_partial$treat_factor,
                                   X = data_partial[, t0.covariate.names])

balwts <- 1/forest_probs$predictions[
  cbind(1:n, as.numeric(data_partial$treat_factor))
]
balwts <- balwts[data_partial$attrited == 0]

data_complete <- data_partial[data_partial$attrited == 0,]

# Conditional on completing the pre-test response module and having been 
# assigned treatment, we estimate propensity to be in each condition AND to not 
# have attrited. We use inverse probability weights as sample weights in the Lin 
# estimator for this sample. 


lm_lin(tolerance_index ~ treat_ind,
       covariates = 
         formula(paste('~', paste(t0.covariate.names, collapse = ' + '))),
       data = data_complete, 
       weights = balwts)

