# covid-travel-testing

Code for analysis of arrival travel testing data in French Polynesia.

_Note: this is currently a working repository, so non-archived code and data are likely to change over time_

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/covid-travel-testing/")
`
The main script to reproduce the analysis is in `scripts/main_script.r`. This calls the following R files:

> `R/data_load.R` - Script to load and format incidence data for travellers arriving in French Polynesia, cases in US/France, and published antigen and PCR test positivity estimates.

> `R/model_functions.R` - Functions to simulate travel testing protocols under different scenarios, estimate depature prevalence, and run statistical models on prevalence data

> `R/plotting_functions.R` - Functions to plot travel testing data, simulation recovery analysis, and departure prevalence estimates.

### Archived versions

• Code to accompany the medRxiv V1 pre-print is in `V1_code`, with file paths as above.

### Citation

[Adam J Kucharski, Kiyojiken Chung, Maite Aubry, Iotefa Teiti, Anita Teissier, Vaea Richard, Timothy W Russell, Raphaëlle Bos, Sophie Olivier, Van-Mai Cao-Lormeau. Real-time surveillance of international SARS-CoV-2 prevalence using systematic traveller arrival screening: An observational study. PLOS Medicine, 2023]([https://www.medrxiv.org/content/10.1101/2022.10.12.22280928v1](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1004283))

The paper above also uses epidemiolgical and population data from additional sources for context (details in manuscript and in `R/data_load.R` file on this repository). Please ensure you cite all relevant data sources if subsequently building on our analysis.
