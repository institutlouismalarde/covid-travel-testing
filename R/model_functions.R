# - - - - - - - - - - - - - - - - - - - - - - - 
# Functions for FP traveller analysis
# - - - - - - - - - - - - - - - - - - - - - - -


# Plot data and binom CI --------------------------------------------------

plot_CI <- function(dates,xx,nn,colA="black") {

  # Remove NA entries in denominators if needed
  dates_p <- ymd(dates[!is.na(nn) & nn>0])
  xx_p <- xx[!is.na(nn) & nn>0] %>% round()
  nn_p <- nn[!is.na(nn) & nn>0] %>% round()
  
  for(ii in 1:length(nn_p)){
    test_r <- binom.test(xx_p[ii],nn_p[ii])
    CI1 <- as.numeric(test_r$conf.int)[1]
    CI2 <- as.numeric(test_r$conf.int)[2]
    points(dates_p[ii],100*xx_p[ii]/nn_p[ii],col=colA,pch=19);
    lines(c(dates_p[ii],dates_p[ii]),100*c(CI1,CI2),col=colA)
  }
  
}

# Output mid and 95% interval --------------------------------------------------

c.text<-function(x1,x2,x3,sigF=2){
  bp1=signif(c(x1,x2,x3),sigF)
  paste(bp1[1],"% (",bp1[2],"-",bp1[3],"%)",sep="")
}

# Calculate CI for cumulative incidence --------------------------------------------------

cumul_CI <- function(sims,cumulative=T) {
  
  # Cumulative sums
  if(cumulative==T){
    sim_cumul <- apply(sims,2,cumsum)
  }else{
    sim_cumul <- sims
  }
  
  # Calculate 95% interval
  sim_95 <- apply(sim_cumul,1,function(x){quantile(x,c(0.025,0.5,0.975))}) %>% t()
  
  # sim_95 <- apply(sims,1,function(x){quantile(x,c(0.025,0.5,0.975))})
  # 
  # # Cumulative sums
  # if(cumulative==T){
  #   sim_cumul <- apply(sim_95,1,cumsum)
  # }else{
  #   sim_cumul <- apply(sim_95,1,cumsum)
  # }
  
  cumul_out <- data.frame(lower=sim_95[,1],median=sim_95[,2],upper=sim_95[,3])
  cumul_out
  
}



# Plot GAM for dates ------------------------------------------------------

plot_GAM <- function(dates,x_val,n_val=NULL,return_vals=T,kk=NULL,family_f="binomial",col1="blue",colf = rgb(0,0,1,0.2)){
  # DEBUG   kk <- NULL; dates <- travel_incidence_n$dates[range1];family_f="binomial"; col1="blue"; colf = rgb(0,0,1,0.2); x_val <- pos_counts_fr[range1]; n_val <- tests_fr_s1[range1]

  if(family_f=="binomial"){
    input <- data.frame(date = dates,vals = x_val,n = round(n_val) )
    input$date <- ymd(input$dat)
    input$date <- as.numeric(input$dat)
    input <- input[!is.na(input$n) & input$n>0,] # Remove NA entries
    
    x_date <- input$date
    y_binom <- cbind(input$vals,input$n)
    
    if(is.null(kk)){
      model1 <- gam(y_binom~s(x_date),family = "binomial")
    }else{
      model1 <- gam(y_binom~s(x_date,k=kk),family = "binomial")
    }
    mult_v <- 100 # Percentage multiplier
  }
  
  if(family_f=="gaussian"){
    input <- data.frame(date = dates,vals = x_val)
    input$date <- ymd(input$dat)
    input$date <- as.numeric(input$dat)
    x_date <- input$date
    y_val <- input$vals
    
    if(is.null(kk)){
      model1 <- gam(y_val~s(x_date),family = "gaussian")
    }else{
      model1 <- gam(y_val~s(x_date,k=kk),family = "gaussian")
    }
    mult_v <- 1 # Set multiplier to 1
  }

  x_date_in <- seq( min(input$date),max(input$date),1)
  preds <- predict(model1, newdata = list(x_date=x_date_in), type = "link", se.fit = TRUE)
  
  # yy <- simulate(model1,nsim=10)
  
  critval <- 1.96; upperCI <- preds$fit + (critval * preds$se.fit); lowerCI <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  fitPlotF <- model1$family$linkinv(fit); CI1plotF <- model1$family$linkinv(upperCI);  CI2plotF <- model1$family$linkinv(lowerCI)
  
  x_date2 <- as.Date(x_date_in)
  polygon(c(x_date2,rev(x_date2)),mult_v*(c(CI1plotF,rev(CI2plotF))),col=colf,lty=0)
  lines(x_date2, mult_v*(fitPlotF) ,col=col1,lwd=2)
  
  # Simulate incidence from GAM prevalence estimates
  sims <- simulate(model1, nsim = 1e3, seed = 42) 
  sims_inc <- apply(sims,1,prev_to_inc) %>% t()
  
  
  if(return_vals==T){
    list(pred_date=x_date2,pred_med=100*fitPlotF,
         pred_CI1=mult_v*CI1plotF,pred_CI2=mult_v*CI2plotF,
         sim_date = as.Date(x_date),
         sim_out = mult_v*sims,
         sim_inc = mult_v*sims_inc)
  }

}


# Plot polygon for GAM outputs ------------------------------------------------------

plot_polygon <- function(x_date_in,fitPlotF,CI1plotF,CI2plotF,scaleF,col1="blue",colf = rgb(0,0,1,0.2)){
 
  x_date2 <- x_date_in
  polygon(c(x_date2,rev(x_date2)),scaleF*(c(CI1plotF,rev(CI2plotF))),col=colf,lty=0)
  lines(x_date2, scaleF*(fitPlotF) ,col=col1,lwd=2)

}


# Estimate growth from data (deprecated) -----------------------------------------------

growth_calc <- function(x,dd=7,type="back",fill_val=0){
  
  d0 <- dd-1
  y <- log(tail(x,-d0)/head(x,-d0))/dd
  if(type=="back"){y <- c(rep(fill_val,d0),y)  }
  if(type=="mid"){y <- c(rep(fill_val,floor(d0/2)),y,rep(fill_val,ceiling(d0/2)))   }
  
  y
  
}

# Growth phase function ----------------------------------------------------

phase_p <- function(x_day,growth=0,maxd){ exp(-growth*x_day)/sum(exp(-growth*(1:maxd))) }


# Forward calculation function -----------------------------------------------

# Generate distribution for missed at departure

# Test T days later
delay_test <- function(n_depart_x = 310, # Number of daily departing infections
                       before_travel_test = 3, # Test 1: days pre-departure
                       after_arrival_test = 4, # Test 1: days pre-departure
                       daily_growth=0, # Daily growth rate at departure
                       dep_test="PCR"
                       ){

  # n_depart_x = 31; before_travel_test = 3; after_arrival_test = 0; daily_growth=0; dep_test="PCR"
  
  # Set up vector of days
  max_days <- nrow(p_by_day)
  x_days <- 1:max_days
  
  # Define departure distribution
  if(dep_test=="PCR"){
    if(before_travel_test==0){
      neg_depart <- p_by_day$p_neg
    }else{
      neg_depart <- c(rep(1,before_travel_test),head(p_by_day$p_neg,-before_travel_test)) # Probability missed at departure given test before
    }
      neg_dist <- neg_depart # Distribution over days since infection in testing
  }
  
  if(dep_test=="Antigen"){
    if(before_travel_test==0){
      neg_depart <- l_by_day$p_neg
    }else{
      neg_depart <- c(rep(1,before_travel_test),head(l_by_day$p_neg,-before_travel_test)) # Probability missed at departure given test before
    }
    neg_dist <- neg_depart # Distribution over days since infection in testing
  }
  
  # Calculate proportion arriving, including epidemic phase shift
  phase_shift <- phase_p(1:max_days,growth=daily_growth,maxd=max_days)
  n_arrive <- n_depart_x*phase_shift*neg_dist # Number who arrive having been missed
  
  if(after_arrival_test==0){
    out_arrive_delay <- n_arrive # Probability missed at arrival on day X
  }else{
    out_arrive_delay <- c(rep(0,after_arrival_test),head(n_arrive,-after_arrival_test)) # Probability missed at arrival on day X
  }
  out_arrive_detect <- out_arrive_delay*p_by_day$median # Number detected
  
  list(n_total_dep = n_depart_x*phase_shift,
       n_total = n_arrive,
       n_arrive_delay=out_arrive_delay, # Number at arrival
       n_arrive_detect=out_arrive_detect) # Number detected
}



# Forward incidence simulation --------------------------------------------

simulate_prev <- function(n_incidence, # Time series of daily departing infections
                         btt = 3, # Test 1: days pre-departure
                         aat = 4, # Test 1: days pre-departure
                         n_max, # Number of time points
                         dep_test="PCR"
                         ){

  max_d <- 30 # PCR positivity max days in calc

  prop_incidence <- n_incidence # Proportion infected per day
  
  n_travellers <- rep(500,n_max)

  # Departure prevalence:
  prev_depart <- rep(0,n_max)
  
  # Rolling window tallying up incidence
  for(ii in 1:n_max){
    calc_window <- ii:min(n_max,ii+max_d) # Window to calculate prevalence
    prev_ii <- prop_incidence[ii]*p_by_day$median # Use median
    prev_ii_cut <- prev_ii[1:min(max_d+1,(n_max-ii+1))] # Avoid overrunning end of time series
    prev_depart[calc_window] <- prev_depart[calc_window] + prev_ii_cut # Match to length

  }
  
  # Define departure positives:
  prev_negative_dep <- rep(0,n_max)
  prev_pos_dep <- rep(0,n_max)
  prev_positive_arr <- rep(0,n_max)
  prev_tested_arr <- rep(0,n_max)
  
  # Define departure distribution
  if(dep_test=="PCR"){
    neg_depart <- c(rep(1,btt),head(p_by_day$p_neg,-btt)) # Probability missed at departure given test before
    neg_dist <- neg_depart # Distribution over days since infection in testing
  }

  # Calculate proportion arriving based on incidence
  for(ii in 1:n_max){ # Iterate over day of arrival
    pick_window <- max(ii-max_d,1):ii # Incidence values to use
    calc_window <- 1:(min(max_d+1,ii)) # Window of test values to use
    negative_dep_ii <- rev((prop_incidence)[pick_window])*neg_dist[calc_window] # Negative at departure by time of incidence
    
    prev_pos_dep[ii] <- sum(rev((prop_incidence)[pick_window])*(1-neg_dist)[calc_window])
    prev_negative_dep[ii] <- sum(negative_dep_ii)
    
    # Arrival detect:
    if(aat==0){
      out_arrive_delay <- negative_dep_ii # Probability missed at arrival on day X
    }else{
      out_arrive_delay <- c(rep(0,aat),head(negative_dep_ii,-aat)) # Probability missed at arrival on day X
    }

    arrive_cut <- length(out_arrive_delay) # Adjust for censoring at start of time series
    
    # XX REFACTOR TO INCLUDE DENOMINATORS IN BINOMIAL
    # Shift to match day of departure and arrival test prevalence in timeseries
    if(ii<(n_max-aat+1)){
      prev_tested_arr[ii+aat] <- sum(negative_dep_ii)/sum((prop_incidence)[pick_window]) # Proportion tested at arrival
      prev_positive_arr[ii+aat] <- sum(out_arrive_delay*p_by_day$median[1:arrive_cut]) # Detected at arrival
    }
    
  }
  
  list(in_incidence=prop_incidence,
       out_pos_depart=prev_depart,
       out_pos_arrive=prev_positive_arr
       )
  
}

# Back calculation function -----------------------------------------------

estimate_vals <- function(n_detect,
                          before_travel_test=3,
                          after_arrival_test=4,
                          daily_growth=0,dep_test="PCR"){ 
  
  # DEBUG: n_detect=1; before_travel_test=3; after_arrival_test=4; daily_growth=0; dep_test="PCR"
  
  max_days <- nrow(p_by_day)
  n_depart <- 1

  # Define departure distribution
  if(dep_test=="PCR"){
    neg_depart <- c(rep(1,before_travel_test),head(p_by_day$p_neg,-before_travel_test)) # Probability missed at departure given test before
    neg_dist <- neg_depart
  }
  
  if(dep_test=="Antigen"){
    neg_depart <- c(rep(1,before_travel_test),head(l_by_day$p_neg,-before_travel_test)) # Probability missed at departure given test before
    neg_dist <- neg_depart
  }
  
  # HAVE DEPRECATED BELOW GROWTH ADJUSTMENTS
  # Add daily growth calculation depending on whether vector - rev data to look backwards
  if(length(daily_growth)==1){ phase_shift <- phase_p(1:max_days,growth=daily_growth,maxd=max_days) }
  if(length(daily_growth)>=max_days){phase_shift <- rev(tail(daily_growth,max_days))}
  if(length(daily_growth)>1 & length(daily_growth)<max_days){ phase_shift <- rev(c(rep(0,max_days-length(daily_growth)),daily_growth)) }

  phase_shift <- phase_shift/sum(phase_shift)
  
  # This is u_i in manuscript
  n_arrive <- n_depart*phase_shift*neg_dist # Number who arrive having been missed, days since infection
 
  if(after_arrival_test==0){
    n_arrive_delay <- n_arrive # Arrival on day X
  }else{
    n_arrive_delay <- c(rep(0,after_arrival_test),head(n_arrive,-after_arrival_test)) # Arrival on day X before test
  }

  # This is a_i in manuscript
  n_arrive_detect <- n_arrive_delay*p_by_day$median # Distribution of timings of those detected
  
  # - - -
  # Start with those detected
  
  # Normalise detected numbers at arrival to get distribution of times
  norm_detect <- n_detect*n_arrive_detect/sum(n_arrive_detect)
  
  # Estimate scaling of under-detection at point of test based on different points in time
  n_arrive_delay_est <- norm_detect/p_by_day$median
  
  # This is r_arrive = sum(n_arrive_delay_est) in manuscript
  
  # Scale number who would have tested positive on day 0 by estimate of total who could have been detected
  if(after_arrival_test==0){
    n_arrive_est <- n_arrive*sum(n_arrive_delay_est)/sum(n_arrive)  # removed 
  }else{
    # Adjust for time shift - as this would reduce detection independent of test performance
    n_arrive_est <- n_arrive*sum(n_arrive_delay_est)/sum(head(n_arrive,-after_arrival_test)) #c(tail(n_arrive_delay_est,-after_arrival_test))
  }
  
  # Estimate total overall at point of departure test
  n_depart_est <- n_arrive_est/neg_dist
  
  # Estimate total that would have been positive in departure survey
  n_depart_survey <- n_depart_est*p_by_day$median
  
  list(#n_depart_est = sum(n_depart_est), 
       #n_arrive_est = sum(n_arrive_est), # Not used
       #ratio_depart = sum(n_depart_est)/n_detect,
       #ratio_arrive = sum(n_arrive_est)/n_detect,
       ratio_total = sum(n_depart_survey)/n_detect # ratio of those that would have been detected at departure
       )
  
}

# Convert prevalence to incidence using weekly values -----------------------------------------

prev_to_inc <- function(x,d_growth=0,weekly_adjust=7){
  
  # Define maximum PCR positivity considered
  max_d <- 31

  # Calculate mean duration of positivity and overall probability density of positivity
  mean_pcr <- round(sum((0:30)*p_by_day$median/sum(p_by_day$median)))
  pcr_density <- sum(p_by_day$median) # Duration of positivity
  
  # Estimate incidence (note need to then shift so that incidence estimates are mean_pcr ahead of prevalence)
  # And scale to weekly value
  weekly_adjust*x/pcr_density

}


# Convert incidence to prevalence using daily values via deconvolution ----------------------

prev_to_inc_daily <- function(prevalence){
  
  # Define transition matrix - 
  f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
  n_pcr_days <- p_by_day$median
  
  for(ii in 1:n_inf){
    i_max <- min(ii+n_pcr_days,n_inf)
    j_max <- min(n_inf-ii+1,n_pcr_days)
    
    f_matrix[ii:i_max,ii] <- p_by_day$median[1:j_max]
    
  }
  
  invert_f <- ginv(f_matrix)
  estimated_incidence <- invert_f %*% prevalence
  
}



