# GS_adult_multistate
Example R code to support multi-state modeling for adult Gulf Sturgeon

Gulf Sturgeon CJFAS reproducibility archive
This archive reproduces the statistical analyses, manuscript tables, and manuscript figures for the Gulf Sturgeon multistate analysis. 

What is included
run_all.R runs the complete workflow in order.
code/00_setup.R checks software and required raw inputs.
code/01_build_informed_histories.R applies transmitter-longevity censoring and constructs the adult annual encounter histories.
code/02_fit_seven_river_models.R fits/loads the 25 seven-river candidate multistate models and the profile-CI version of the top model.
code/03_deployment_accounting.R reconstructs adult-tag deployment accounting used in manuscript summaries.
code/04_build_tag_pools.R constructs annual active-transmitter tag pools under first-observed and most-recent-river assignment rules.
code/06_fit_regional_state_model.R fits the separate four-region multistate model used for regional movement/fidelity summaries.
code/07_build_manuscript_tables.R creates manuscript and appendix tables.
code/08_build_manuscript_figures.R creates Figures 2, 3, A1, and A2 plus figure-data CSV files.
code/_project_paths.R contains portable project paths and Program MARK configuration.
Figure 1 (study-area map) is produced separately and is not recreated by this workflow. 

Required raw data
These two source files must go in data/raw/:

MS_TR_1350mmTL_allATc1_MAX_ANNUAL_2010-2022_7basin_20240121.inp
admin_query_export_All_Fish_All_ATags_By_Mark_Status_Min_Date_20240121.csv
The .inp file is the archived annual encounter-history input and already reflects the upstream no-singles filtering and annual encounter-state construction. The CSV contains acoustic-tag deployment/longevity information used to apply transmitter-life censoring.

Because of the ESA listing status of Gulf Sturgeon and absence of data sharing agreements, data requests must be filed with NOAA Fisheries by contacting Nick Farmer (nick.farmer@noaa.gov)

Software requirements
R
Program MARK, called through the R package RMark
R packages: RMark, tidyverse, dplyr, stringr, tidyr, ggplot2, gridExtra, and grid

Running the complete analysis
Start R and run:

source("run_all.R")
The scripts are also designed to be run individually in numerical order.

Workflow and generated files
01_build_informed_histories.R reads the two raw archived inputs and regenerates the key derived objects, including:

data/derived/pre2010.RDS
data/derived/informed.all.stages.ch.RDS
data/derived/informed_at_longevity_all_gs_ms_export_20240121.RDS
data/derived/informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS
02_fit_seven_river_models.R creates model objects in results/models/seven_river/. Existing model objects in that directory are loaded rather than refit, which makes repeated runs faster.

03_deployment_accounting.R and 04_build_tag_pools.R create additional derived RDS objects used by manuscript products.

06_fit_regional_state_model.R writes results/models/regional/regional_state_model.RDS. This model has a different state space/likelihood from the seven-river candidate set; its AICc is not compared with the seven-river model-selection table.

The outputs are written to:

outputs/tables/
outputs/figures/
results/audits/


Reproducibility notes
The workflow uses the archived 2024-01-21 database exports and does not query a live database.
Missing manufacturer tag longevity is assigned five years, preserving the validated analysis rule.
The annual transmitter endpoint convention is preserved from the validated workflow.
freq = -1 in the MARK input represents loss-on-capture/censoring; sum(freq) is not the fish sample size.
Fish first entering the adult analysis in 2022 are excluded because they have no subsequent annual survival interval.
Direct geographically impossible transitions are fixed to zero as specified in the model scripts.
For p(river:time) models, final-occasion detection is constrained equal to the preceding occasion within river to address terminal survival-detection confounding.
The regional-state model is a separate likelihood/state-space model used for regional movement/fidelity summaries and is not part of the seven-river AICc candidate set.
