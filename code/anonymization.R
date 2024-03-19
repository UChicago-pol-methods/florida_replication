
## read in data 
df <- read.csv("../data/sender_db_may15.csv",colClasses = rep("character", ncol(read.csv("../data/sender_db_may15.csv"))))

## read in and set seed
seed <- readLines('../data/seed.txt') # keep this private and local
set.seed(seed)

# recode userids
ids <- unique(df$sender_id)
keys <- sample(length(ids))

df$sender_id <- keys[match(df$sender_id, ids)]

## writeout data
saveRDS(df, "../data/df_for_analysis_final.rds")
