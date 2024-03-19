## load in pacakges 
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(jsonlite)
library(stringi)


###---------------------------description---------------------------------### 
## this R code performs valid flexible output preprocessing tasks on concern_db, including concatenating user input data, 
## joint with sender_db to get sender text and evaluate classification results. 


###-------------------------------concern_db-------------------------------## 
### processing valid flexible output 
## read in data
df_concern<-read.csv("../data/concern_db_may7.csv",colClasses = rep("character", ncol(read.csv("../data/concern_db_may7.csv"))))

## sort by time
df_concern <- df_concern%>% 
  arrange(sender_id, timestamp_utc)

df_concern <- df_concern %>% 
  group_by(sender_id) %>% 
  slice_tail(n = 1)

## get valid flexible input (string length > 40)
flexible_valid<-df %>%
  filter(section == "Flexible") %>%
  filter(section_id %in% c("2")) 

flexible_valid<-flexible_valid %>%
  filter(str_length(sender_text) > 40)

flexible_valid <- flexible_valid %>% 
  arrange(sender_id, desc(timestamp_utc))

## keep the most recent user input
flexible_valid <- flexible_valid %>% 
  group_by(sender_id) %>% 
  slice_tail(n = 1)

# ## if it's from the same user, concatenate the input 
# flexible_valid$merged_sender_text <- ""
# last_sender <- ""
# for (i in 1:nrow(flexible_valid)) {
#   if (flexible_valid[i, "sender_id"] == last_sender) {
#     flexible_valid[i, "merged_text"] <- paste(flexible_valid[i-1, "merged_text"], flexible_valid[i, "sender_text"], sep = ";")
#   } else {
#     flexible_valid[i, "merged_text"] <- flexible_valid[i, "sender_text"]
#     last_sender <- flexible_valid[i, "sender_id"]
#   }
# }


## filter sender_id in valid input
df_concern<-df_concern %>%
  filter(sender_id %in% flexible_valid$sender_id)

df_concern<-df_concern %>%
  distinct(sender_id, concern_response,.keep_all = TRUE)

## join to get sender_text 
df_concern<-left_join(df_concern[],flexible_valid[,-1], by = "sender_id")


## remove ok, no, yes 
df_concern$sender_text <- gsub("^(ok|yes|no)[\\p{P}\\s]*", "", df_concern$sender_text)

## remove leading white space or punctuation 
df_concern$sender_text<- stri_replace_first_regex(df_concern$sender_text, "^\\p{P}*\\s*", "")

## remove tailing ok 
df_concern$sender_text <- gsub("\\bok\\b\\s*$", "", df_concern$sender_text)

# ## remove duplicates 
# df_concern<-df_concern %>%
#   distinct(sender_id, .keep_all = TRUE)

## sort by time 

df_concern$timestamp_utc<-as.POSIXct(df_concern$timestamp_utc)
df_concern<-df_concern %>%
  arrange(timestamp_utc)

## write out 
current_date <- format(Sys.Date(), "%Y-%m-%d")
file_name <- paste0("valid_flexible_input_", current_date, ".csv")
file_path <- "../data"
full_file_path <- file.path(file_path, file_name)
write.csv(df_concern,  file = full_file_path)
