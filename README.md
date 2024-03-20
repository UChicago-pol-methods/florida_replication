# Florida replication
Repository for Replication files for "Canvassing with chatbots: personalized messaging using natural language processing to change attitudes".

# Instructions for Replication


## Datasets and data cleaning
Before starting the replication, ensure the the below files are contained in a `/data` folder in the parent directory.

- `/data/df_for_analysis_final.rds`

In the `/code` folder, compile `pre_processing_senderdb.R`, which will produce clean working data sets that are also located in `/data`:
   - `/data/df_for_analysis_processed.rds`

## Analysis
Compile the following scripts in sequence:
  - `/code/main_analysisr.R`, for tables and figures in the paper.
  - `/code/estimations_w_bounds.R`, `/code/bootstrapping_manski.R` for tables and figures in the appendix. 
