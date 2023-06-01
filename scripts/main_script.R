# - - - - - - - - - - - - - - - - - - - - - - - 
# Main code for FP traveller
# Author: Adam Kucharski
# - - - - - - - - - - - - - - - - - - - - - - -

# Load libraries
library(dplyr)
library(magrittr)
library(readr)
library(zoo)
library(mgcv)
library(readxl)
library(lubridate)
library(incidence)
library(gridExtra)
library(DescTools)
library(devtools)
library(MASS)
library(data.table)

# Load library for highest posterior density interval calculation
# devtools::install_github("rmcelreath/rethinking@slim")
library(rethinking)

# Set working directory
# setwd("~/Documents/GitHub/covid-travel-testing/")

# Load and clean data and functions
source("R/data_load.R")

# Load model and plotting functions
source("R/model_functions.R")
source("R/plotting_functions.R")

# Figure 1
figure_schematic_testing(after_arrival_test=4)

# Figure 2 + Figure S1
figure_phase_bias(btt=3,att=4) # Defined with before and after travel test days

# Figure 3
figure_basic_data()

# Figure 4
figure_reconstruct_epidemics()

# Figure S2 - antigen test 1 day before travel (+ 1 day in transit)
figure_reconstruct_epidemics(test_type2="Antigen",btt=2)

# Figure S3
figure_explore_cumulative_incidence()

