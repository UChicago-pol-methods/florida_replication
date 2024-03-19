## load in pacakges 
library(data.table)
library(dplyr)
library(jsonlite)
library(tidyr)
library(xtable)

##--------------------------------------- description ---------------------------------------##
## this R code performs several data preprocessing tasks on sender_db, including concatenating user input data, 
## pivoting a long-format data frame into a wide-format data frame to better analyze the data, 
## and recoding variable names and values to improve clarity. 


###---------------------------------------sender_db---------------------------------------------## 
## read in data 
df <- readRDS("../data/df_for_analysis_final.rds")

## subset to users who responded after 19/Apr/2023 16:54:48 
df$timestamp_utc <- as.POSIXct(df$timestamp_utc)
df <- df %>%
  filter(timestamp_utc > "2023-04-19 16:55:00")


## make sure each user only answer the survey once (filtered by the panel question)
sender_id_fileter<-df %>%
  filter(section == "Panel1") %>%
  count(sender_id, section, section_id) %>%
  filter(n!=1) %>%
  select(sender_id) %>%
  unique()

sender_id_unique<-df %>%
  filter(sender_id %in% sender_id_fileter$sender_id) %>%
  arrange(sender_id, timestamp_utc) %>%
  distinct(sender_id, section, section_id, .keep_all = TRUE)

sender_id_unique$sender_id<-as.character(sender_id_unique$sender_id)

df<-df %>%
  filter(!sender_id %in% sender_id_unique$sender_id)
  
df<-rbind(df, sender_id_unique)

## ---------------------------------------concatenation----------------------------------------## 
## concatenate open response 
datalist <- list()
for (i in list(c("Video", "4","5","1"),c("Video","6","7","2"),c("Define","4","5","3"),
               c("Define","8","9","4"),c("Flexible","2","3","5"),c("Address","7","8","6")) ){
  subset_df <- subset(df, section == i[1] & section_id %in% c(i[2], i[3]))
  agg_df <- aggregate(sender_text ~ sender_id, subset_df, paste, collapse = ";")
  agg_df$section = i[1]
  agg_df$section_id = i[2]
  datalist[[i[4]]] <- agg_df
}

big_data <- do.call(rbind, datalist)
## replace old value with concatenated value
big_data <- setDT(big_data)
setDT(df)[big_data, sender_text:= i.sender_text,  on = .(sender_id, section,section_id)]


##---------------------------------------pivot from long to wide--------------------------------------- ##
df_wide <- df %>%
  select(-c(chatbot_text_seq,chatbot_text,redirection_id,redirection_section, timestamp_utc)) %>%
  group_by(sender_id) %>%
  pivot_wider(names_from = c(section, section_id), values_from = sender_text)

## ---------------------------------------re-order column--------------------------------------- ##
## column name
column_name <- fromJSON("../supporting_data/linear_dict.json")
results <- character()
for (section in names(column_name)) {
  for (subsection in names(column_name[[section]])) {
    results <- c(results, paste0(section, "_", subsection))
  }
}
results = c("sender_id", results)
df_wide <- df_wide[,results]
df_wide <- as.data.frame(df_wide)
## replace null with NA
df_wide[df_wide == "NULL"] <- NA
df_wide[df_wide == "NA"] <- NA

## --------------------------------------- re-code variable name ---------------------------------------##
## re-code variable name - primary outcome
old_primary_colname <- c("Positions_4", "Positions_8","Positions_9","Positions_10","Positions_11","Positions_12", "Positions_13","Positions_14")
new_primary_colname <- c("therm_trans_t1","gender_norm_sexchange_t1" ,"gender_norm_moral_t1","gender_norm_abnormal_t1","gender_norm_trans_moral_wrong_t1","trans_teacher_t1",
                       "trans_bathroom_t1","gender_norm_dress_t1")
## re-code variable name - secondary outcome
old_secondary_colname <- c("Positions_1", "Positions_2")
new_secondary_colname <- c("florida_trans_policy_t1", "florida_trans_policy2_t1")
## re-code variable name -covariates
old_covariates_colname <- c("Panel1_5", "Panel1_6","Panel1_2","Panel1_3","Panel1_4","Panel1_7","Panel1_8",
                          "Panel1_9","Panel1_10","Panel1_11")
new_covariates_colname <- c("age_t0","gender_t0","ideology_t0","pid_t0","pol_interest_t0","healthcare_t0",
                          "climate_t0","religion_t0","abortion_t0","immigration_t0")
## change primary column name
names(df_wide) <- setNames(lapply(names(df_wide), function(x) {
  ifelse(x %in% old_primary_colname, new_primary_colname[old_primary_colname == x], x)
}), colnames(df_wide))
## change secondary column name
names(df_wide) <- setNames(lapply(names(df_wide), function(x) {
  ifelse(x %in% old_secondary_colname, new_secondary_colname[old_secondary_colname == x], x)
}), colnames(df_wide))
## change covariates column name
names(df_wide) <- setNames(lapply(names(df_wide), function(x) {
  ifelse(x %in% old_covariates_colname, new_covariates_colname[old_covariates_colname == x], x)
}), colnames(df_wide))


## --------------------------------------- variable coding &  reverse coding ---------------------------------------##
## reverse coding variables:
## 1. gender_norm_moral_t1 Coding: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree
## 2. gender_norm_abnormal_t1: Coding: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree
## 3. gender_norm_trans_moral_wrong_t1: Coding:  -1 = Agree, 0 = Noopinion/don’t know, 1 = Disagree
## 4. trans_teacher_t1: Coding: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree
## 5. trans_bathroom_t1: Coding: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree
## 6. gender_norm_dress_t1: [Coding: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree
df_analysis <- df_wide %>%
  select(sender_id, therm_trans_t1,gender_norm_sexchange_t1,gender_norm_moral_t1,gender_norm_abnormal_t1,gender_norm_trans_moral_wrong_t1,
         trans_teacher_t1,trans_bathroom_t1,gender_norm_dress_t1,florida_trans_policy_t1,florida_trans_policy2_t1,age_t0,gender_t0,ideology_t0,pid_t0,pol_interest_t0,
         healthcare_t0,climate_t0,religion_t0,abortion_t0,immigration_t0,Control_1,Context_1 )
## gender_norm_sexchange_t1
df_analysis$gender_norm_sexchange_t1 <- case_when(df_analysis$gender_norm_sexchange_t1 == "agree" ~ -1,
                                                  df_analysis$gender_norm_sexchange_t1 == "disagree" ~ 1,
                                                  df_analysis$gender_norm_sexchange_t1 == "undecided/don't know" ~ 0,
                                                  df_analysis$gender_norm_sexchange_t1 == "NA" ~ NA_integer_,
                                                  TRUE ~ 999)
## gender_norm_moral_t1
df_analysis$gender_norm_moral_t1 <- case_when(df_analysis$gender_norm_moral_t1 == "agree" ~ -1,
                                              df_analysis$gender_norm_moral_t1 == "disagree" ~ 1,
                                              df_analysis$gender_norm_moral_t1 == "undecided/don't know" ~ 0,
                                              df_analysis$gender_norm_moral_t1 == "NA" ~ NA_integer_,
                                              TRUE ~ 999)
## gender_norm_abnormal_t1
df_analysis$gender_norm_abnormal_t1 <- case_when(df_analysis$gender_norm_abnormal_t1 == "agree" ~ -1,
                                                 df_analysis$gender_norm_abnormal_t1 == "disagree" ~ 1,
                                                 df_analysis$gender_norm_abnormal_t1 == "undecided/don't know" ~ 0,
                                                 df_analysis$gender_norm_abnormal_t1 == "NA" ~ NA_integer_,
                                                 TRUE ~ 999)
## gender_norm_trans_moral_wrong_t1
df_analysis$gender_norm_trans_moral_wrong_t1 <- case_when(df_analysis$gender_norm_trans_moral_wrong_t1 == "agree" ~ -1,
                                                          df_analysis$gender_norm_trans_moral_wrong_t1 == "disagree" ~ 1,
                                                          df_analysis$gender_norm_trans_moral_wrong_t1 == "undecided/don't know" ~ 0,
                                                          df_analysis$gender_norm_trans_moral_wrong_t1 == "NA" ~ NA_integer_,
                                                          TRUE ~ 999)
# trans_teacher_t1
df_analysis$trans_teacher_t1 <- case_when(df_analysis$trans_teacher_t1 == "agree" ~ -1,
                                          df_analysis$trans_teacher_t1 == "disagree" ~ 1,
                                          df_analysis$trans_teacher_t1 == "undecided/don't know" ~ 0,
                                          df_analysis$trans_teacher_t1 == "NA" ~ NA_integer_,
                                          TRUE ~ 999)
# trans_bathroom_t1
df_analysis$trans_bathroom_t1 <- case_when(df_analysis$trans_bathroom_t1 == "agree" ~ -1,
                                           df_analysis$trans_bathroom_t1 == "disagree" ~ 1,
                                           df_analysis$trans_bathroom_t1 == "undecided/don't know" ~ 0,
                                           df_analysis$trans_bathroom_t1 == "NA" ~ NA_integer_,
                                           TRUE ~ 999)
# gender_norm_dress_t1
df_analysis$gender_norm_dress_t1 <- case_when(df_analysis$gender_norm_dress_t1 == "agree" ~ -1,
                                              df_analysis$gender_norm_dress_t1 == "disagree" ~ 1,
                                              df_analysis$gender_norm_dress_t1 == "undecided/don't know" ~ 0,
                                              df_analysis$gender_norm_dress_t1 == "NA" ~ NA_integer_,
                                              TRUE ~ 999)
## florida_trans_policy_t1
df_analysis$florida_trans_policy_t1 <- case_when(df_analysis$florida_trans_policy_t1 == "1" ~ -3,
                                                 df_analysis$florida_trans_policy_t1 == "2" ~ -2,
                                                 df_analysis$florida_trans_policy_t1 == "3" ~ -1,
                                                 df_analysis$florida_trans_policy_t1 == "4" ~ 0,
                                                 df_analysis$florida_trans_policy_t1 == "5" ~ 1,
                                                 df_analysis$florida_trans_policy_t1 == "6" ~ 2,
                                                 df_analysis$florida_trans_policy_t1 == "7" ~ 3,
                                                 df_analysis$florida_trans_policy_t1 == "NA" ~ NA_integer_,
                                                 TRUE ~ 999)
## florida_trans_policy_t2
df_analysis$florida_trans_policy2_t1 <- case_when(df_analysis$florida_trans_policy2_t1 == "1" ~ 3,
                                                  df_analysis$florida_trans_policy2_t1 == "2" ~ 2,
                                                  df_analysis$florida_trans_policy2_t1 == "3" ~ 1,
                                                  df_analysis$florida_trans_policy2_t1 == "4" ~ 0,
                                                  df_analysis$florida_trans_policy2_t1 == "5" ~ -1,
                                                  df_analysis$florida_trans_policy2_t1 == "6" ~ -2,
                                                  df_analysis$florida_trans_policy2_t1 == "7" ~ -3,
                                                  df_analysis$florida_trans_policy2_t1 == "NA" ~ NA_integer_,
                                                  TRUE ~ 999)
## covariates
df_analysis$gender_t0 <- case_when(df_analysis$gender_t0 == "man" ~ 1,
                                   df_analysis$gender_t0 == "woman" ~ 0,
                                   df_analysis$gender_t0 == "something else" ~ 0,
                                   df_analysis$gender_t0 == "NA" ~ NA_integer_,
                                   TRUE ~ 999)
df_analysis$ideology_t0 <- case_when(df_analysis$ideology_t0 == "1" ~ -3,
                                     df_analysis$ideology_t0 == "2" ~ -2,
                                     df_analysis$ideology_t0 == "3" ~ -1,
                                     df_analysis$ideology_t0 == "4" ~ 0,
                                     df_analysis$ideology_t0 == "5" ~ 1,
                                     df_analysis$ideology_t0 == "6" ~ 2,
                                     df_analysis$ideology_t0 == "7" ~ 3,
                                     df_analysis$ideology_t0 == "NA" ~ NA_integer_,
                                     TRUE ~ 999)
df_analysis$pid_t0 <- case_when(df_analysis$pid_t0 == "1" ~ -3,
                                df_analysis$pid_t0 == "2" ~ -2,
                                df_analysis$pid_t0 == "3" ~ -1,
                                df_analysis$pid_t0 == "4" ~ 0,
                                df_analysis$pid_t0 == "5" ~ 1,
                                df_analysis$pid_t0 == "6" ~ 2,
                                df_analysis$pid_t0 == "7" ~ 3,
                                df_analysis$pid_t0 == "NA" ~ NA_integer_,
                                TRUE ~ 999)
df_analysis$pol_interest_t0 <- case_when(df_analysis$pol_interest_t0 == "1" ~ 1,
                                         df_analysis$pol_interest_t0 == "2" ~ -1,
                                         df_analysis$pol_interest_t0 == "3" ~ 0,
                                         df_analysis$pol_interest_t0 == "4" ~ -2,
                                         df_analysis$pol_interest_t0 == "NA" ~ NA_integer_,
                                         TRUE ~ 999)
df_analysis$healthcare_t0 <- case_when(df_analysis$healthcare_t0 == "1" ~ 1,
                                       df_analysis$healthcare_t0 == "2" ~ 2,
                                       df_analysis$healthcare_t0 == "3" ~ 3,
                                       df_analysis$healthcare_t0 == "4" ~ 4,
                                       df_analysis$healthcare_t0 == "NA" ~ NA_integer_,
                                       TRUE ~ 999)
df_analysis$religion_t0 <- case_when(df_analysis$religion_t0 == "1" ~ -1,
                                     df_analysis$religion_t0 == "2" ~ 0,
                                     df_analysis$religion_t0 == "3" ~ 1,
                                     df_analysis$religion_t0 == "NA" ~ NA_integer_,
                                     TRUE ~ 999)
df_analysis$abortion_t0 <- case_when(df_analysis$abortion_t0 == "1" ~ 1,
                                     df_analysis$abortion_t0 == "2" ~ 2,
                                     df_analysis$abortion_t0 == "3" ~ 3,
                                     df_analysis$abortion_t0 == "4" ~ 4,
                                     df_analysis$abortion_t0 == "NA" ~ NA_integer_,
                                     TRUE ~ 999)
df_analysis$immigration_t0 <- case_when(df_analysis$immigration_t0 == "1" ~ -1,
                                        df_analysis$immigration_t0 == "2" ~ 1,
                                        df_analysis$immigration_t0 == "NA" ~ NA_integer_,
                                        TRUE ~ 999)

### NA: user didn't answer; 999: user answered something else other than the provided answer
## fill NA or blankness with NA

df_analysis[] <- lapply(df_analysis, function(x) if(is.list(x)) as.character(unlist(x)) else x)
df_analysis[df_analysis=="NA"] <- NA
df_analysis[df_analysis==""] <- NA


##--------------------------------------- user input cleaning --------------------------------------- ##
## extract as much useful information from input as possible
# Function to extract the age number from the answer
extractAge <- function(answer) { 
  age <- gsub("\\D", "", answer)  # Remove non-digit characters
  age <- as.numeric(age)  # Convert the result to numeric
  return(age)
}

df_analysis$age_t0<- sapply(df_analysis$age_t0, extractAge)

## convert from character to numeric
df_analysis[] <- lapply(df_analysis, function(x) if(is.list(x)) as.character(unlist(x)) else x)
df_analysis[, -c(1,22,23)] <- apply(df_analysis[, -c(1,22,23)], 2, function(x) as.numeric(as.character(x)))


# ## format control -  temporarily treating all answers which don't follow the coding rubric as NA
# df_analysis$therm_trans_t1[!(df_analysis$therm_trans_t1 %in% c(0:10,999))] <- NA
# df_analysis$gender_norm_sexchange_t1[!(df_analysis$gender_norm_sexchange_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$gender_norm_moral_t1[!(df_analysis$gender_norm_moral_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$gender_norm_abnormal_t1[!(df_analysis$gender_norm_abnormal_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$gender_norm_trans_moral_wrong_t1[!(df_analysis$gender_norm_trans_moral_wrong_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$trans_teacher_t1[!(df_analysis$trans_teacher_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$trans_bathroom_t1[!(df_analysis$trans_bathroom_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$gender_norm_dress_t1[!(df_analysis$gender_norm_dress_t1 %in% c(-1,0,1,999))] <- NA
# df_analysis$florida_trans_policy_t1[!(df_analysis$florida_trans_policy_t1 %in% c(-3:3,999))] <- NA
# df_analysis$florida_trans_policy2_t1[!(df_analysis$florida_trans_policy2_t1 %in% c(-3:3,999))] <- NA
# df_analysis$gender_t0[!(df_analysis$gender_t0 %in% c(1,0,999))] <- NA
# df_analysis$ideology_t0[!(df_analysis$ideology_t0 %in% c(-3:3,999))] <- NA
# df_analysis$pid_t0[!(df_analysis$pid_t0 %in% c(-3:3,999))] <- NA
# df_analysis$pol_interest_t0[!(df_analysis$pol_interest_t0 %in% c(-2:1,999))] <- NA
# df_analysis$healthcare_t0[!(df_analysis$healthcare_t0 %in% c(1:3,999))] <- NA
# df_analysis$climate_t0[!(df_analysis$climate_t0 %in% c(0:10,999))] <- NA
# df_analysis$religion_t0[!(df_analysis$religion_t0 %in% c(-1,0,1,999))] <- NA
# df_analysis$abortion_t0[!(df_analysis$abortion_t0 %in% c(1:4,999))] <- NA
# df_analysis$immigration_t0[!(df_analysis$immigration_t0 %in% c(-1,0,1,999))] <- NA

# Get the column names which end with 't0'
varnames_t0 <- grep("t0$", names(df_analysis), value=TRUE)
# Remove 'healthcare_t0' and 'abortion_t0' from the list
varnames_t0 <- varnames_t0[!varnames_t0 %in% c('healthcare_t0', 'abortion_t0')]

# Loop over each variable and replace 999 with mean for t0 variable
for (var in varnames_t0) {
  mean_value <- mean(df_analysis[[var]][df_analysis[[var]] != 999], na.rm = TRUE)
  df_analysis[[var]][df_analysis[[var]] == 999] <- mean_value
  df_analysis[[var]][is.na(df_analysis[[var]])] <- mean_value
}

### replace 999 with NA for categorical variable in t0 variable
df_analysis[c('healthcare_t0', 'abortion_t0')][df_analysis[c('healthcare_t0', 'abortion_t0')]  == 999] <- NA

### raplace 999 with NA for t1 variable
varnames_t1 <- grep("t1$", names(df_analysis), value = TRUE)
df_analysis[varnames_t1][df_analysis[varnames_t1] == 999] <- NA


df_analysis$treat_ind <- NA
df_analysis$treat_ind[!is.na(df_analysis$Context_1)] <- TRUE
df_analysis$treat_ind[!is.na(df_analysis$Control_1)] <- FALSE

## keep observation in either treatment or control
df_analysis <- df_analysis %>%
  filter(!is.na(treat_ind))

## finish all questions; primary: tolerence index; sec: law index; therm 
df_analysis <- df_analysis %>%
  mutate(finished_dv_primary = ifelse(
    rowSums(is.na(df_analysis[, 2:9])) == 0 & rowSums(df_analysis[, 2:9] == 999) == 0,
    TRUE,
    FALSE
  ))

df_analysis <- df_analysis %>%
  mutate(finished_dv_sec = ifelse(
    rowSums(is.na(df_analysis[, 10:11])) == 0 & rowSums(df_analysis[, 10:11] == 999) == 0,
    TRUE,
    FALSE
  ))

df_analysis <- df_analysis %>%
  mutate(finished_dv_therm_thrans = ifelse(
    !is.na(df_analysis$therm_trans_t1) & df_analysis$therm_trans_t1 != 999,
    TRUE,
    FALSE
  ))


## pre-survey finish rate
n_complete_presurvey_rows <- sum(!is.na(df_analysis$treat_ind))
pre_survey_finish_rate<-n_complete_presurvey_rows/nrow(df_wide)*100
pre_survey_n<-nrow(df_wide)


## conditional filtering
## remove age > 99 and < 18, dropped 39 observations
df_analysis <- df_analysis %>%
  filter(age_t0>=18 & age_t0 <= 99)

## remove therm_trans_t1 not in c(0:10, 999), dropped 1 observations 
df_analysis <- df_analysis %>%
  filter(therm_trans_t1 %in% c(0:10,999)| is.na(therm_trans_t1))


##--------------------------------------- calculate survey finish rate---------------------------------------##

## among people who finished the pre-survey question, calculate the randomization possibility
n_control_group <- sum(df_analysis[, 24] != TRUE)
n_treatment_group <- sum(df_analysis[, 24] == TRUE)
percent_control<-n_control_group/(n_treatment_group+n_control_group)*100
percent_treatment<-n_treatment_group/(n_treatment_group+n_control_group)*100
## finish rate for each group
n_complete_control_rows <- sum((df_analysis[, 24] != TRUE) & (df_analysis[, 25] == TRUE & df_analysis[, 26] == TRUE))
finish_rate_control<-n_complete_control_rows/n_control_group*100
n_complete_treatment_rows <-  sum((df_analysis[, 24] == TRUE) & (df_analysis[, 25] == TRUE & df_analysis[, 26] == TRUE))
finish_rate_treatment<-n_complete_treatment_rows/n_treatment_group*100

summary_table <- data.frame(
  Group = c("Pre-survey","Control", "Treatment"),
  N = c(pre_survey_n, n_control_group, n_treatment_group),
  "Percent of Group" = c("NA", percent_control, percent_treatment),
  "Finish rate" = c(pre_survey_finish_rate, finish_rate_control,finish_rate_treatment )
)

summary_table_latex <- xtable(summary_table, caption = "Attrition Summary")
names(summary_table_latex)<-c("Group","N","Percentage", "Finish rate")

# Save the summary table as a LaTeX file
print(summary_table_latex, file = "../tables/attrition_summary.tex")

## ---------------------------------- write-out ---------------------------------##
saveRDS(df_analysis, "../data/df_for_analysis_processed.rds")
