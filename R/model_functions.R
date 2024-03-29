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


# Grid lines by year ------------------------------------------------------

grid_year <- function(){
  y1 <- as.Date("2021-01-01")
  y2 <- as.Date("2022-01-01")
  
  lines(c(y1,y1),c(0,1e6),lty=3,col="light gray")
  lines(c(y2,y2),c(0,1e6),lty=3,col="light gray")
  
}

# Plot data and binom CI --------------------------------------------------

plot_CI_def <- function(dates,data_in,colA="black") {
  
  xx <- data_in[,1]
  c_lower <- data_in[,2]
  c_upper <- data_in[,3]
  
  for(ii in 1:length(dates)){
    points(dates[ii],xx[ii],col=colA,pch=19);
    lines(c(dates[ii],dates[ii]),c(c_lower[ii],c_upper[ii]),col=colA)
  }
  
}

# Output mid and 95% interval --------------------------------------------------

# c.text<-function(x1,x2,x3,sigF=2){
#   bp1=signif(c(x1,x2,x3),sigF)
#   paste(bp1[1],"% (",bp1[2],"-",bp1[3],"%)",sep="")
# }

c.text<-function(x1,x2,x3,sigF=2){
  bp1=signif(c(x1,x2,x3),sigF)
  paste(bp1[1],"% (95% CrI: ",bp1[2],", ",bp1[3],"%)",sep="")
}



# Estimate MLE and 95% profile likelihood ---------------------------------

prev_estimate <- function(data_in, 
                          before_travel_test=3, # 
                          after_arrival_test=4,
                          test_type="PCR",
                          sim_out_n = 1e3 # 
                          ){
  
  # before_travel_test=3; after_arrival_test=4; sim_out_n=100
  
  y <- as.numeric(data_in[1]) # Number positive at arrival
  x <- as.numeric(data_in[2]) # Number tested at arrival

  # Calculate detection proportion at arrival
  n_depart_N <- nrow(p_by_day)
  output_n_det <- delay_test(n_depart_N,before_travel_test,after_arrival_test,dep_test=test_type)$n_arrive_detect; output_n_det <- output_n_det[output_n_det>0]
  
  # Define probabilities
  p_n1_p2 <- sum(output_n_det)/n_depart_N # P(negative at departure & positive at arrival | infected)
  
  # Check whether use PCR or LFT distribution for departure
  if(test_type=="PCR"){
    p_n1 <- sum(1-p_by_day$median)/n_depart_N   # P(negative at departure | infected)
  }else{
    p_n1 <- sum(1-l_by_day$median)/n_depart_N   # P(negative at departure | infected)
  }
  p_n1_pcr <- sum(1-p_by_day$median)/n_depart_N   # Define PCR for prevalence

  # Calculate profile likelihood for w (i.e. probability of infection at departure)
  xx_search <- seq(0,1,0.0001)
  likelihood_f <- function(w){dbinom(y,x,w*p_n1_p2/(w*p_n1+(1-w)),log=T)}
  pp_out <- sapply(xx_search,likelihood_f)
  
  # Extract MLE and 95% CI:
  mle_est <- y/(x*p_n1_p2 + y - y*p_n1) #xx_search[which.max(yy_out)]
  prof_lik <- xx_search[(max(pp_out)-pp_out)<qchisq(0.95,1)/2]
  
  prev_est <- 100*c(mle_est,min(prof_lik),max(prof_lik))*(1-p_n1_pcr) # convert back into probability test positive)

  posterior_p <- exp(pp_out)/sum(exp(pp_out)) # Define posterior of p (assuming uniform prior)
  
  # Calculate variance of p
  mean_p <- sum(xx_search*posterior_p)
  variance_p <- sum((xx_search-mean_p)^2*posterior_p) # Define variance of p
                                 
  sim_out <- sample(xx_search,sim_out_n,prob = posterior_p,replace=T)*(1-p_n1_pcr) # convert back into probability test positive)

  list(prev_est=prev_est,
       map_est=mle_est,
       var_est=variance_p,
       sim_out=sim_out
  )
}

# Plot bootstrap for prevalence ------------------------------------------------------

bootstrap_est <- function(dates,
                          data_input, # First column positives at arrival, second column tests
                          return_vals=T,
                          bootstrap_n=1e3,
                          before_travel_test=3,
                          after_arrival_test=4,
                          test_type="PCR",
                          kk=NULL,
                          gam_add=T,
                          col1="blue",
                          colf = rgb(0,0,1,0.2)){
  
  # DEBUG before_travel_test=3; after_arrival_test=4; bootstrap_n=1e2; kk=NULL;data_input=data_fr1; dates=travel_incidence_n$dates[range1];col1="blue"; colf = rgb(0,0,1,0.2);

  # Format dates and remove data points with no tests
  n_test <- data_input[,2]
  x_date <- dates
  if(class(x_date)=="Date"){x_date <- ymd(dates)}
  x_date <- as.numeric(x_date)
  
  valid_x <- (!is.na(n_test) & n_test>1)  # Remove NA entries
  x_date <- x_date[valid_x]
  data_in <- data_input[valid_x,]

  # Set up bootstrap matrix
  store_sim <- matrix(NA,nrow=length(x_date),ncol=bootstrap_n)
  prev_store <- matrix(NA,nrow=length(x_date),ncol=4) # Final entry is variance
  
  # Simulate prevalence from posterior for each time poinnt
  for(ii in 1:length(x_date)){
    est_vals <- prev_estimate(data_in[ii,],before_travel_test,after_arrival_test,test_type,sim_out_n=bootstrap_n)
    sim_ii <- est_vals$sim_out # simulated values
    store_sim[ii,] <- sim_ii
    prev_store[ii,] <- c(est_vals$prev_est,est_vals$var_est)
  }
  
  # - - -
  # Fit weighted GAMs for estimate data point, weighted by inverse variance of each estimate
  x_date_range <- seq( min(x_date),max(x_date),1)

  mult_v <- 100 # Set multiplier to 100 (e.g. percent vs proportion)

  # Adjust for lognormal fitting
  y_val <- prev_store[,1]
  # y_val[y_val==0] <- 1e-2
  # y_val <- log(y_val)
  #   
  weight_var <- 1/prev_store[,4]; weight_var <- weight_var/mean(weight_var)

  model1 <- gam(y_val~s(x_date),family = "gaussian",weights=weight_var)

  # Predict curve from data
  preds <- predict(model1, newdata = list(x_date=x_date_range), type = "link", se.fit = TRUE)

  critval <- 1.96; upperCI <- preds$fit + (critval * preds$se.fit); lowerCI <- preds$fit - (critval * preds$se.fit)
  store_gam <- preds$fit
  CI1plotF <- model1$family$linkinv(lowerCI)
  CI2plotF <- model1$family$linkinv(upperCI)
  CI1plotF <- pmax(0,CI1plotF) # constrain to be positive

  # Calculate 95% HPDI for parameter
  sim_95 <- apply(mult_v*store_sim,1,function(x){rethinking::HPDI(x,prob=0.95)}) %>% t()
  CI1plotP <- sim_95[,1]
  #fitPlotP <- sim_95[,2]
  CI2plotP <- sim_95[,2]
  

  # Plot curves and output estimates from posterior
  x_date2 <- as.Date(x_date_range)
  
  #plot(x_date,prev_store[,1]) # XX DEBUG - REMOVE LATER XX
  plot_CI_def(x_date,cbind(prev_store[,1],CI1plotP,CI2plotP)) # Use MAP and HDPI
  
  # Overlay GAM of fitted points
  if(gam_add==T){
    polygon(c(x_date2,rev(x_date2)),(c(CI1plotF,rev(CI2plotF))),col=colf,lty=0)
    lines(x_date2, (store_gam) ,col=col1,lwd=2)
  }
  
  # Simulate incidence from bootstrap prevalence estimates
  sims <- simulate(model1, nsim = 1e3, seed = 42); sims[sims<0] <- 0 # truncate
  sims_inc <- apply(sims,1,prev_to_inc) %>% t()
  
  #sims_inc <- apply(mult_v*store_sim,1,prev_to_inc) %>% t()

  
  if(return_vals==T){
    list(pred_date=x_date,
         pred_med=store_gam,
         pred_CI1=CI1plotF,
         pred_CI2=CI2plotF,
         out_MAP = prev_store[,1],
         HDPI_1 = CI1plotP,
         HDPI_2 = CI2plotP,
         sim_date = as.Date(x_date),
         #sim_out = sims,
         sim_inc = sims_inc)
  }
  
  # CUT CODE
  
  # Plot posterior draws
  # for(kk in 1:100){
  #   points(x_date,mult_v*store_sim[,kk],pch=19,cex=0.5,col=rgb(0.5,0.5,0.5,0.1))
  # }
  
  
}

# Plot GAM for dates ------------------------------------------------------
# DEPRECATED?

plot_GAM <- function(dates,x_val,n_val=NULL,return_vals=T,kk=NULL,family_f="binomial",col1="blue",colf = rgb(0,0,1,0.2)){
  # DEBUG   kk <- NULL; dates <- travel_incidence_n$dates[range1];family_f="binomial"; col1="blue"; colf = rgb(0,0,1,0.2); x_val <- pos_counts_fr[range1]; n_val <- tests_fr_s1[range1]
  # DEBUG dates = x_weeks; x_val=prev_pos_dep[w_s]; n_val=rep(n_travel_daily,length(w_s))
  
  
  if(family_f=="binomial"){
    input <- data.frame(date = dates,vals = x_val,n = round(n_val) )
    if(class(input$date)=="Date"){input$date <- ymd(input$dat)}
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
  # By definition simulation uses original timescale of data
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

    # Shift to match day of departure and arrival test prevalence in timeseries
    if(ii<(n_max-aat+1)){
      prev_tested_arr[ii+aat] <- sum(negative_dep_ii)/sum((prop_incidence)[pick_window]) # Proportion tested at arrival
      prev_positive_arr[ii+aat] <- sum(out_arrive_delay*p_by_day$median[1:arrive_cut]) # Detected at arrival
    }
    
  }
  
  # CHECKS
  prev_pos_dep[50]/(prev_positive_arr[50]/prev_tested_arr[50])
  
  list(in_incidence=prop_incidence,
       out_pos_depart=prev_depart,
       out_pos_arrive=prev_positive_arr
       )
  
}


# Simulate and recover arrival test data ----------------------------------------------

# Predict PCR positivity distribution on arrival
stochastic_arrival <- function( n_incidence, # Time series of daily departing infections
                                btt = 3, # Test 1: days pre-departure
                                att = 4, # Test 1: days pre-departure
                                n_max, # Number of time points
                                n_travel_vol = 2000 ,# Number of weekly travellers
                                test_type="PCR",
                                test_type_infer="PCR",
                                col_list_in = col_list2[[2]],
                                add_gam = F,
                                x_weeks,
                                w_s
                                ){
                                
  # DEBUG n_incidence = 0.2*dnorm(x_days,mean=50,sd=22) ; btt=3; att=4; n_max=150; test_type="PCR"; test_type_infer="PCR"

  # Set seed for simulations
  set.seed(10)
  
  # Set up variables and lengths
  max_d <- 30 # PCR positivity max days in calc (includes 0)
  n_travellers <- rep(n_travel_vol,n_max)
  
  prop_incidence <- n_incidence # Proportion infected on given day
  
  n_depart_N <- nrow(p_by_day)   # Total days in PCR curve

  #prev_depart <- rep(0,n_max) 

  # XX NEED TO ACCOUNT FOR EPIDEMIC PHASE HERE...
  
  # Rolling window tallying up incidence
  # for(ii in 1:n_max){
  #   calc_window <- ii:min(n_max,ii+max_d) # Window to calculate prevalence
  #   prev_ii <- prop_incidence[ii]*p_by_day$median # Use median
  #   prev_ii_cut <- prev_ii[1:min(max_d+1,(n_max-ii+1))] # Avoid overrunning end of time series
  #   prev_depart[calc_window] <- prev_depart[calc_window] + prev_ii_cut # Match to length
  #   
  # }
  
  # Define variables for positives and number tested:
  prev_negative_dep <- rep(0,n_max)
  prev_pos_dep <- rep(0,n_max)
  prev_positive_arr <- rep(0,n_max)
  prev_tested_arr <- rep(0,n_max)
  e_pos_dep <- rep(0,n_max)
  
  # Define departure distribution
  if(test_type=="PCR"){
    neg_test <- p_by_day$p_neg # Probability missed at departure given test before
    neg_depart <- c(rep(1,btt),head(neg_test,-btt)) # Probability missed at departure given test before
    neg_dist <- neg_depart # Distribution over days since infection in testing
  }
  if(test_type!="PCR"){
    neg_test <- 1-test_type*(1-p_by_day$p_neg) # 80% sensitivity
    neg_depart <- c(rep(1,btt),head(neg_test,-btt)) # Probability missed at departure given test before
    neg_dist <- neg_depart # Distribution over days since infection in testing
  }
  
  # Calculate proportion arriving based on incidence
  for(ii in 1:n_max){ # Iterate over day of arrival
    pick_window <- max(ii-max_d,1):ii # Incidence values to use
    calc_window <- 1:(min(max_d+1,ii)) # Window of test values to use
    n_window <- length(calc_window)
    
    # Tally up those who would test negative at departure
    #n_travel_d <- round(n_travel_vol/(max_d+1)) # Distribute over possible infection days
    
    n_inf <- rbinom(n_window,n_travel_vol,rev((prop_incidence)[pick_window])) # Number infected at departure based on recent incidence
    
    negative_dep_ii <- rbinom(n_window,n_inf,neg_dist[calc_window]) # Negative at departure by time of incidence
    positive_dep_ii <- sum(n_inf - negative_dep_ii)
    
    # Expected positives if tested at departure - use 'neg_test' here as no shift
    departure_e_pos <- n_travel_vol*rev((prop_incidence)[pick_window])*(1-neg_test[calc_window]) # Prevalence at departure test (rather than at departure with btt lag)
    e_pos_dep[ii] <- sum(departure_e_pos)
    
    # Tally positives and negatives at departure
    prev_negative_dep[ii] <- sum(negative_dep_ii)
    prev_pos_dep[ii] <- sum(positive_dep_ii)
    
    # Arrival detect:
    if(att==0){
      out_arrive_delay <- negative_dep_ii # Probability missed at arrival on day X
    }else{
      out_arrive_delay <- c(rep(0,att),head(negative_dep_ii,-att)) # Probability missed at arrival on day X
    }
    
    arrive_cut <- length(out_arrive_delay) # Adjust for censoring at start of time series
    
    # XX REFACTOR TO INCLUDE DENOMINATORS IN BINOMIAL
    # Shift to match day of departure and arrival test prevalence in timeseries
    # Note burn in if testing after day 0 arrival
    if(ii<(n_max-att+1)){
      prev_tested_arr[ii+att] <- n_travel_vol - sum(positive_dep_ii) # Number tested at arrival
      prev_positive_arr[ii+att] <- sum(rbinom(arrive_cut,out_arrive_delay,p_by_day$median[1:arrive_cut])) # Detected at arrival
    }
    
  }
  
  # XX DEBUG
  # Plot by week
  plot(x_weeks,100*(e_pos_dep/n_travel_vol)[w_s],type="l",ylab="%",col="white",ylim=c(0,5),xaxs="i",yaxs="i",xlab="days",main="")

  data_week <- data.frame(prev_positive_arr,prev_tested_arr); data_week <- data_week[w_s,]
  
  # Estimate prevalence (assume baseline 3 day pre arrival and 4 day post)
  est_prev <- bootstrap_est(dates=x_weeks, data_input=data_week,gam_add=add_gam,
                            before_travel_test=3,after_arrival_test=4,test_type=test_type_infer)

  points(x_weeks,100*(e_pos_dep/n_travel_vol)[w_s],lwd=2,col=col_list_in) # Expected prevalence
  
  points(x_weeks,100*prev_positive_arr[w_s]/prev_tested_arr[w_s],col=col_list_in,pch=0) # Prevalence at arrival
  
  # Plot GAM comparison on original data with larger sample size:
  if(add_gam==T){
    plot_GAM(x_weeks,round(e_pos_dep[w_s]),rep(n_travel_vol,length(w_s)),
             return_vals=T,kk=NULL,family_f="binomial",col1="orange",colf = rgb(1,0.5,0,0.2))
  }

  
  # correlation
  # break_n <- seq(-5.25,5.25,0.5)
  # hist(100*tail(prev_pos_dep,-4)/n_travel_vol-est_prev$out_MAP,breaks=break_n)


  list(in_incidence=prop_incidence,
       out_test_depart = n_travel_vol,
       out_pos_depart = prev_pos_dep,
       out_test_arrive = prev_tested_arr,
       out_pos_arrive=prev_positive_arr
  )
  
  
}

# Back calculation function (deprecated) -----------------------------------------------

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
  pcr_density <- sum(p_by_day$median) # Density of positivity
  
  # Estimate incidence (note need to then shift so that incidence estimates are mean_pcr ahead of prevalence)
  # And scale to weekly incidence value
  weekly_adjust*x/pcr_density

}


# Convert prevalence to incidence using daily values via deconvolution ----------------------

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



