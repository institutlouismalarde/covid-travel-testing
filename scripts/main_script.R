# - - - - - - - - - - - - - - - - - - - - - - - 
# Main code for FP traveller
# Author: Adam Kucharski
# - - - - - - - - - - - - - - - - - - - - - - -

# Load libraries
library(tidyverse)
library(zoo)
library(mgcv)
library(readxl)
library(lubridate)
library(incidence)
library(gridExtra)
library(DescTools)
library(gratia)
library(MASS)

# Load and clean data and functions
source("R/data_load.R")

# Load model and plotting functions
source("R/model_functions.R")
source("R/plotting_functions.R")

# Figure 1
figure_schematic_testing(after_arrival_test=4)

# Figure 2
figure_phase_bias(btt=2,att=4) # Defined with before and after travel test days

# Figure 3
figure_basic_data()

# Figure 4
figure_reconstruct_epidemics()

# Figure S1 - antigen test 1 day before travel (+ 1 day in transit)
figure_reconstruct_epidemics(test_type="Antigen",btt=2)

# Figure S2
figure_explore_cumulative_incidence()

