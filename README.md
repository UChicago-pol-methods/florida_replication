# Florida Replication

Repository for replication files for "Deep canvassing with automated conversational agents: personalized messaging to change attitudes".


## Replication Instructions

### Step 0: Data Preprocessing
Ensure `/data/df_for_analysis_final.rds` exists, then run:
```bash
cd code
Rscript pre_processing_senderdb.R
```
**Outputs**: 
- `/data/df_for_analysis_processed.rds` (5,115 observations, 31 variables)
- `/tables/attrition_summary.tex` - Summary of pre-survey, control, and treatment group sizes and completion rates

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
Generate Manski bounds with bootstrapped confidence intervals:
```bash
Rscript bootstrapping_manski.R     # Bootstrap confidence intervals
```
**Note**: If `boot_UB.rds` and `boot_LB.rds` already exist in `/data/`, the script will use them by default to quickly generate the table. To force regeneration of bootstrap results (e.g., after changing parameters), set `USE_EXISTING_BOOTSTRAP <- FALSE` at the top of the script.

**Outputs**: 
- `/data/boot_UB.rds`, `/data/boot_LB.rds` - Bootstrap results (if generated)
- `/tables/manski_bounds_ci.tex` - Manski bounds with 95% confidence intervals

### Step 4: Generate Figures
Run the figures script to create visualizations:
```bash
mkdir -p ../figures  # Create figures directory if it doesn't exist
Rscript figures.R
```
**Outputs**:
- `/figures/treatment_effects.png` - Treatment effects plot for all respondents (3 outcomes × 3 models)
- `/figures/treatment_effects_by_subgroup.png` - Treatment effects by party affiliation (Overall, Democrats, Republicans, Independents)

### Step 5: Generate Codebook (Optional)
```bash
Rscript datamaid_codebook.R
```
**Output**: `/test-data-codebook.Rmd` - Data dictionary for all variables
