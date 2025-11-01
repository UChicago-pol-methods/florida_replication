# Florida Replication

Repository for replication files for "Canvassing with chatbots: personalized messaging using natural language processing to change attitudes".


## Replication Instructions

### Step 0: Data Preprocessing
Ensure `/data/df_for_analysis_final.rds` exists, then run:
```bash
cd code
Rscript pre_processing_senderdb.R
```
**Output**: `/data/df_for_analysis_processed.rds` (5,115 observations, 31 variables)

### Step 1: Main Analysis
Run the main analysis script (automatically sources `utils.R`):
```bash
Rscript main_analysis.R
```
**Outputs**:
- `/tables/main_results.tex` - Main treatment effects (all respondents)
- `/tables/dem_results.tex` - Subgroup analysis (Democrats)
- `/tables/rep_results.tex` - Subgroup analysis (Republicans)
- `/tables/ind_results.tex` - Subgroup analysis (Independents)
- `/tables/custom_summary_*.tex` - Custom formatted results tables

### Step 2: Robustness Checks
Run robustness analysis:
```bash
Rscript robustness.R
```
**Outputs**:
- `/tables/balance_table_all.tex` - Covariate balance (all respondents)
- `/tables/balance_table_completers.tex` - Covariate balance (completers)
- `/tables/placebo_leave_one_out.tex` - Leave-one-out covariate placebo test
- `/tables/placebo_outcomes.tex` - Placebo outcomes test

### Step 3: Appendix Analysis (Manski Bounds)
Run 
```bash
Rscript estimations_w_bounds.R    # Manski bounds estimation
Rscript bootstrapping_manski.R     # Bootstrap confidence intervals
```
**Outputs**: Various Manski bounds tables in `/tables/`

### Step 5: Generate Codebook (Optional)
```bash
Rscript datamaid_codebook.R
```
**Output**: `/test-data-codebook.Rmd` - Data dictionary for all variables
