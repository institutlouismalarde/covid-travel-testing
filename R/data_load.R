# - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and format
# - - - - - - - - - - - - - - - - - - - - - - -

# Load and clean Our World In Data France case data -------------------------------------------------

owid_france_clean <- function(fr_cases){
  date_list <- which(is.na(fr_cases$new_cases))
  yy <- fr_cases
  
  for(ii in 1:length(date_list)){
    jj <- date_list[ii]
    yy[jj,]$new_cases <- mean(c(yy[jj-1,]$new_cases,yy[jj+1,]$new_cases))
  }
  
  # Fix remaining two entries
  list_rem <- which(is.na(yy$new_cases))
  yy[list_rem,]$new_cases <- mean(c(yy[min(list_rem)-1,]$new_cases,yy[max(list_rem)+1,]$new_cases))
  
  yy
  
}


# Load US infection-induced seroprevalence data ---------------------------------------------------
# Source: https://covid.cdc.gov/covid-data-tracker/#national-lab

# Load data and calculate midpoint of serosurvey
us_national_sero <- read_csv("data/us_antibody_data.csv") |> 
                    mutate(date_mid = date_start + round(as.numeric(date_end-date_start)/2))


# Load wastewater data ----------------------------------------------------

wastewater_nat <- fread("https://raw.githubusercontent.com/biobotanalytics/covid19-wastewater-data/master/wastewater_by_region.csv")



# Load tourism data -------------------------------------------------------
# Source: https://data.ispf.pf/themes/SystemeProductif/Tourisme/Details.aspx
# tourism_data <- read_csv("data/Tourism_data_ISPF.csv")

# Define protocol periods --------------------------------------------------

# COV-CHECK dates
covd1 <- as.Date(c("2020-07-15","2021-04-30"))
covd2 <- as.Date(c("2021-08-12","2022-03-27"))
covd1_q <- as.Date(c("2021-02-20","2021-04-30"))

# LABM dates
labd <- as.Date(c("2021-05-01","2021-08-11"))

# Import incidence values --------------------------------------------------

# Weekly incidence for all tests - weeks start on Mondays - - - 
travel_incidence_n <- read_csv("data/inc_travel_incidence_n.csv")

# COV-CHECK tests by origin - - - 
passenger_tests_1 <- read_csv("data/inc_passenger_tests_1.csv")
passenger_tests_2 <- read_csv("data/inc_passenger_tests_2.csv")

# Lab flight origins by day - - - 
traveller_labm <- read_csv("data/inc_traveller_labm.csv")

# COV-CHECK positives by country of origin - - - 
country_incidence <- read_csv("data/inc_country_incidence.csv")

# Ct values - - - 
travel_test_ct <- read_csv("data/inc_travel_test_ct.csv")


# Local reported cases - - - 
fp_cases <- read_csv("data/inc_fp_cases.csv")


# - - - 
# Calculate test origins

# Calculate test origins using season 1 COV-CHECK exact data:
match_date <- match(travel_incidence_n$dates,passenger_tests_1$dates)

tests_fr_s1 <- passenger_tests_1$France[match_date]
tests_us_s1 <- passenger_tests_1$USA[match_date]

# Calculate test origins using season 2 COV-CHECK exact data:
match_date <- match(travel_incidence_n$dates,passenger_tests_2$dates)

tests_fr_s2 <- passenger_tests_2$France[match_date]
tests_us_s2 <- passenger_tests_2$USA[match_date]

# Split into three fits for COV-CHECK and lab testing
range1 <- which(!is.na(tests_fr_s1))
range2 <- which(!is.na(tests_us_s2))
range_lab <- (max(range1)+1):(min(range2)-1) # Extra week of overlap with S2

# - - -
# Estimate LABM travel tests by origin

# Convert to weekly
traveller_labm2 <- traveller_labm %>% mutate(week = ceiling_date(date, "week",week_start=1))
traveller_wk <- traveller_labm2 %>% group_by(week) %>% summarize(fr_tot=sum(fr_in),us_tot=sum(us_in),tot_tot=sum(total_in))

# Match and estimate travellers
match_date <- match(travel_incidence_n$dates[range_lab],traveller_wk$week) # Only set for lab dates
tests_during_period <- rowSums(travel_incidence_n[,c("NON","OUI")])[range_lab]

tests_fr_lab <- (tests_during_period*(traveller_wk$fr_tot/traveller_wk$tot_tot)[match_date]) %>% round()
tests_us_lab <- (tests_during_period*(traveller_wk$us_tot/traveller_wk$tot_tot)[match_date]) %>% round()

# Set to zero if no tests
tests_fr_lab[is.na(tests_fr_lab)] <- 0; tests_us_lab[is.na(tests_us_lab)] <- 0

# Positive counts estimated
pos_counts_fr_lab <- (travel_incidence_n$OUI[range_lab]*(traveller_wk$fr_tot/traveller_wk$tot_tot)[match_date])
pos_counts_us_lab <- (travel_incidence_n$OUI[range_lab]*(traveller_wk$us_tot/traveller_wk$tot_tot)[match_date])


# Load wider data --------------------------------------------------

# PCR positivity and Ct curves from Hellewell et al, BMC Med, 2021
pcr_curve <- read_csv("data/PCR_curve_summary.csv") # PCR curve
lft_curve <- read_csv("data/LFT_curve_summary.csv") # Antigen test curve
hcw_ct <- read_csv("data/hellewell_ct_values_since_exposure.csv") # HCW Ct

# Our World In Data surveillance data for France and US
fr_cases <- read_csv("data/inc_fr_cases.csv")
us_cases <- read_csv("data/inc_us_cases.csv")

# Convert PCR curves to daily probabilities
pcr_curve <- pcr_curve %>% mutate(p_neg = 1-median)
lft_curve <- lft_curve %>% mutate(p_neg = 1-median)

day_list <- c(0:30)
p_by_day <- pcr_curve[match(day_list,pcr_curve$days_since_infection),]
l_by_day <- lft_curve[match(day_list,lft_curve$days_since_infection),]
