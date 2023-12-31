########################################################################
# PIPELINE FOR RMST CALCULATION FOR MULTIPLE EVENT TYPES

## RUN FOR RESULTS ## 

# Prof. LJ Wei, Dr. Zack McCaw, Stephanie Armbruster 
########################################################################


### add visualsation and test for survival profile difference between treatment 
# and control group 

# Installation - first time only ------------------------------------------
# only run once at beginning to install packages
# once installed, skip this step 
# install.packages(c("tidyverse", 
#                    "survRM2", 
#                    "plyr", 
#                    "cowplot"))
# devtools::install_github(repo = "zrmacc/MRMST", 
#                          force = TRUE)


# Load - every time -------------------------------------------------------
# load libraries
# required every time code is run 
library(tidyverse)
library(survRM2)
library(MRMST)
library(survival)
library(reconstructKM)

# set seed to ensure reproducibility 
set.seed(2023)

# Input -------------------------------------------------------------------

# !!!This section requires input depending on data set!!!

# set working directory to where Pipeline folder is saved
setwd("~/.../Pipeline")
source("functions_pipeline.R")

# input path to file 
data_raw <- read_csv("path_to_file.csv", 
                     show_col_types = FALSE)

# insert column name for multiple events of interest
times_to_events <- c("t2hosp", 
                     "t2mi", 
                     "t2trans")

# insert column names for status variables for multiple events
# i.e. 1 if event observed, 0 otherwise 
statuses_for_events <- c("c4hosp", 
                         "c4mi", 
                         "c4trans")

# insert column name for composite event and its status
time_to_death <- "t2dth"
status_for_death <- "c4dth"

# insert column name for treatment indicator 
arm <- "trt"

# insert desired confidence interval coverage
alpha <- 0.05

# insert desired grid density for number at risk output
nar_grid <- 10

# insert rounding precision for number at risk output
round_to <- 50

# Workflow ----------------------------------------------------------------
# 1. Construct composite endpoints 
data_comp <- construct_composite_endpoints(data = data_raw, 
                                           t2event_names = times_to_events,
                                           c4event_names = statuses_for_events,
                                           t2death_name = time_to_death, 
                                           c4death_name = status_for_death,
                                           
                                           composite_name = "Death")

# rename treatment column into "arm"
# data_comp <- rename_arm(data = data_comp, 
#                         arm = arm)
# for PARADISE_MI
data$arm <- data$treatment

# generate column name for composite endpoints 
times_for_composite_events <- paste(times_to_events, 
                                 "Death", 
                                 sep = "")
statuses_for_composite_events <- paste(statuses_for_events, 
                                       "Death", 
                                       sep = "")

  
# 2. Find maximal truncation time 
global_tau <- find_truncation_time(data = data_comp, 
                                   t2event_names = times_for_composite_events, 
                                   c4event_names = statuses_for_composite_events)
# save global truncation time in Data folder 
save(global_tau, 
     file = "Data/global_tau.RData")

# 3. Survival estimate
# compute Kaplan-Meier survival estimate and save plot in Figures folder 
KM_S <- KM_survival(data = data_comp,
                    t2event_names = times_for_composite_events, 
                    c4event_names = statuses_for_composite_events, 
                    
                    alpha = alpha, 
                    nar_grid = nar_grid, 
                    round_to = round_to)

# show KM survival curves, the pointwise confidence interval 
# and the number at risk
# showing plots requires hitting <RETURN> in the console 
par(ask = TRUE)
KM_S$curves
par(ask = FALSE)

# compute Cox proportional hazard model 
Cox <- HR(data = data_comp,
          t2event_names = times_for_composite_events, 
          c4event_names = statuses_for_composite_events)


# save KM survival estimates in Data folder
save(KM_S, 
     file = "Data/KM_survival_estimates.RData")
# save Cox model estimates in Data folder
save(Cox, 
     file = "Data/COX_estimates.RData")

# !! Choose truncation time !!
chosen_tau <- global_tau

# 4. RMST calculation
RMST_results <- MRMST::TwoSample(
  arm = data_comp$arm,
  status = data_comp %>% 
    dplyr::select(all_of(statuses_for_composite_events)), 
  time = data_comp %>% 
    dplyr::select(all_of(times_for_composite_events)),
  tau = chosen_tau, 
  alpha = alpha
)
# save RMST results in Data folder 
save(RMST_results, 
     file = "Data/RMST_estimates.RData")

# 5. WLW procedure
WLW_results <- WLW_TwoSamples(data = data_comp, 
                              time = times_for_composite_events, 
                              status = statuses_for_composite_events)
save(WLW_results, 
     file = "Data/WLW_estimates.RData")

# 6. Stratification 
