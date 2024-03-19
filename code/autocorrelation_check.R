
## autocorrelation check 

df<-read.csv("../data/sender_db_may9.csv",colClasses = rep("character", ncol(read.csv("../data/sender_db_may9.csv"))))

## subset to users who responded after 19/Apr/2023 16:54:48 
df$timestamp_utc<-as.POSIXct(df$timestamp_utc)
df<-df %>%
  filter(timestamp_utc > "2023-04-19 16:55:00")


## ---------------------------------------- concatenation ---------------------------------------## 
datalist = list()
datalist = vector("list", length = n)
for (i in list(c("Video", "4","5","1"),c("Video","6","7","2"),c("Define","4","5","3"),
               c("Define","8","9","4"),c("Flexible","2","3","5"),c("Address","7","8","6")) ){
  subset_df<- subset(df, section == i[1] & section_id %in% c(i[2], i[3]))
  agg_df <- aggregate(sender_text ~ sender_id, subset_df, paste, collapse = ";")
  agg_df$section = i[1]
  agg_df$section_id = i[2]
  datalist[[i[4]]]<-agg_df
}

big_data = do.call(rbind, datalist)

## replace old value with concatenated value 
big_data<-setDT(big_data)
setDT(df)[big_data, sender_text:= i.sender_text,  on = .(sender_id, section,section_id)]


##  -------------------------------- pivot from long to wide  ------------------------------- ## 
df_check<-df %>%
  select(-c(chatbot_text_seq,chatbot_text,redirection_id,redirection_section)) %>%
  group_by(sender_id) %>%
  pivot_wider(names_from = c(section, section_id), values_from = sender_text)

df_check<-df_check %>%
  select(sender_id,timestamp_utc,Control_1,Context_1 )

df_check <- df_check %>%
  mutate(treated = ifelse(!is.na(Context_1), 1, ifelse(!is.na(Control_1), 0, NA)))

df_check<-df_check %>%
  filter(!is.na(treated))

df_check <- df_check[order(df_check$timestamp_utc),]

# Check the autocorrelation of treatment
acf(df_check$treated)
autocorr <- acf(df_check$treated, plot = FALSE)$acf[2]
autocorr
