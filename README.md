# Chabot-replication
Replication files for "Canvassing with chatbots: personalized messaging using natural language processing to change attitudes"

- Pre-registration [here](https://osf.io/yqmjg). 

**Files**
+ `code/` folder
  + Data cleaning script [`pre_processing_senderdb.R`](https://github.com/UChicago-pol-methods/DD-notes-memos/blob/main/chatbot-replication/code/pre_processing_senderdb.R)
  + Analysis scripts [`main_analysis.R`](https://github.com/UChicago-pol-methods/DD-notes-memos/blob/main/chatbot-replication/code/main_analysis.R)
  
**Reproducing results**

1. Save most recent data in the `data/` folder. 
2. To replicate analysis, run
   +  Data cleaning script [`pre_processing_senderdb.R`](https://github.com/UChicago-pol-methods/DD-notes-memos/blob/main/chatbot-replication/code/pre_processing_senderdb.R) then 
   + Analysis scripts [`main_analysis.R`](https://github.com/UChicago-pol-methods/DD-notes-memos/blob/main/chatbot-replication/code/main_analysis.R)
5. Resulting html files (`pre_processing_senderdb.html`, `main_analysis.html`) will be saved in the same folder.
6. Figures and latex tables will be saved in `figures/` and `tables/` folders respectively.
