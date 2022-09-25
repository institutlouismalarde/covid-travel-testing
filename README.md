# covid-travel-testing

Code for analysis of arrival travel testing data in French Polynesia.

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/covid-travel-testing/")
`
The main script to reproduce the analysis is in `scripts/main_script.r`. This calls the following R files:

> `R/data_load.R` - Script to load and format incidence data for travellers arriving in French Polynesia, cases in US/France, and published antigen and PCR test positivity estimates.

> `R/model_functions.R` - Functions to simulate travel testing protocols under different scenarios, estimate depature prevalence, and run statistical models on prevalence data

> `R/plotting_functions.R` - Functions to plot travel testing data, simulation recovery analysis, and departure prevalence estimates.

### Citation

Adam J Kucharski, Kiyojiken Chung, Maite Aubry, Iotefa Teiti, Anita Teissier, Vaea Richard, Timothy W Russell, Raphaëlle Bos, Sophie Olivier, Van-Mai Cao-Lormeau. Real-time surveillance of international SARS-CoV-2 prevalence using systematic traveller arrival screening, 2022
