# - - - - - - - - - - - - - - - - - - - - - - - 
# Plotting functions for FP traveller analysis
# - - - - - - - - - - - - - - - - - - - - - - -

# Figure 1: Schematic with delays and detection at arrival -----------------------------------------------

figure_schematic_testing <- function(n_depart_N=31, 
                                     after_arrival_test=0,
                                     before_travel_test=3) {
  
  # DEBUG  n_depart_N=31; after_arrival_test=4; before_travel_test=3
  
  par(mfrow=c(2,3),mgp=c(2,0.7,0),mar = c(4,3,1,1))
  col_red <- rgb(1,0.5,0.5)
  col_green <- rgb(0.6,1,0.6)
  x_shift <- 6
  x_labels <- 0:30
  
  letter_x <- 2
  
  # PCR detection
  
  plot(pcr_curve$days_since_infection,pcr_curve$median,xlab="days since infection",
       ylab="probability",xaxs="i",col="white",ylim=c(0,1),type="l",yaxs="i",main="Test positivity")
  
  polygon(c(pcr_curve$days_since_infection,rev(pcr_curve$days_since_infection)),c(pcr_curve$lower_95,rev(pcr_curve$upper_95)),col=rgb(0,0,1,0.2),lty=0)
  lines(pcr_curve$days_since_infection,pcr_curve$median,col="blue",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # At departure
  out_depart <- delay_test(n_depart_N,before_travel_test,after_arrival_test)$n_total
  out_arrive <- delay_test(n_depart_N,before_travel_test,after_arrival_test)$n_total
  
  plot(x_labels,x_labels,xlab="days since infection",ylab="probability",
       ylim=c(0,1),type="l",yaxs="i",xaxs="i",col="white",main="Departure testing")
  
  polygon(c(x_labels,rev(x_labels)),c(0*out_arrive,rev(rep(n_depart_N,length(x_labels)))),col=col_green,lty=0)
  polygon(c(x_labels,rev(x_labels)),c(0*out_arrive,rev(out_arrive)),col=col_red,lty=0)
  lines(c(before_travel_test,before_travel_test),c(-10,1e3),lty=2) 
  text(x=x_shift,y=0.05,labels="missed",adj=0,col="dark red")
  text(x=x_shift,y=0.8,labels="detected",adj=0,col="dark green")
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # At arrival
  output_n_all <- delay_test(n_depart_N,before_travel_test,after_arrival_test)$n_arrive_delay
  x_labels2 <- x_labels[output_n_all>0]; output_n_all <- output_n_all[output_n_all>0]
  output_n_det <- delay_test(n_depart_N,before_travel_test,after_arrival_test)$n_arrive_detect; output_n_det <- output_n_det[output_n_det>0]
  
  plot(x_labels,x_labels,xlab="days since infection",ylab="probability",
       ylim=c(0,1),type="l",yaxs="i",xaxs="i",col="white",main="Arrival testing")
  polygon(c(x_labels2,rev(x_labels2)),c(0*output_n_all,rev(output_n_all)),col=col_red,lty=0)
  polygon(c(x_labels2,rev(x_labels2)),c(output_n_all-output_n_det,rev(output_n_all)),col=col_green,lty=0)
  lines(c(after_arrival_test,after_arrival_test),c(-10,1e3),lty=2)
  text(x=x_shift,y=0.8,labels="detected",adj=0,col="dark green")
  text(x=x_shift,y=0.05,labels="missed",adj=0,col="dark red")
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  

  # - - - 
  # Plot incidence simulation
  
  n_max <- 150
  peak_sim <- 1e-04; cut_sim <- 50
  prop_incidence1 <- 0.2*dnorm(1:n_max,mean=50,sd=22) # Proportion infected per day
  prop_incidence2 <-  0.15*dnorm(1:n_max,mean=120,sd=30) # Proportion infected per day
  
  out1 <- simulate_prev(prop_incidence1,btt=1,aat=0,n_max)
  out2 <- simulate_prev(prop_incidence2,btt=1,aat=0,n_max)
  out1_4 <- simulate_prev(prop_incidence1,btt=3,aat=4,n_max)
  out2_4 <- simulate_prev(prop_incidence2,btt=3,aat=4,n_max)
  
  # Plot incidence
  ymax <- 3.5
  plot(100*out1$in_incidence,type="l",ylab="%",col="white",ylim=c(0,0.4),xaxs="i",yaxs="i",xlab="days",main="Incidence")
  grid(nx=NULL,ny=NA,col="light gray")
  lines(100*out1$in_incidence,col="dark orange",lwd=2)
  lines(100*out2$in_incidence,col="dark blue",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

  plot(100*out1$out_pos_depart,type="l",ylab="%",ylim=c(0,ymax),xaxs="i",yaxs="i",col="white",xlab="days",main="Departure prevalence")
  grid(nx=NULL,ny=NA,col="light gray")
  lines(100*out1$out_pos_depart,col="dark orange",lwd=2)
  lines(100*out2$out_pos_depart,col="dark blue",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  plot(100*out1$out_pos_arrive,type="l",ylab="%",ylim=c(0,ymax),col="white",xaxs="i",yaxs="i",xlab="days",main="Arrival prevalence")
  grid(nx=NULL,ny=NA,col="light gray")
  lines(100*out1_4$out_pos_arrive+100*out2_4$out_pos_arrive,col="grey",lwd=2)
  lines(100*out1_4$out_pos_arrive,col="orange",lwd=2)
  lines(100*out2_4$out_pos_arrive,col="blue",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # Output figure
  dev.copy(pdf,paste("plots/Fig1_schematic.pdf",sep=""),width=7,height=5)
  dev.off()
  
  
}


# Figure 2: Growth phase bias -----------------------------------------------

figure_phase_bias <- function(btt=3,att=4,add_gam_in=F){
  
  # btt=3; att=4
  
  # Sample on departure
  n_total <- 1e3
  x_labels <- 0:30
  max_days <- nrow(p_by_day)
  
  # Predict PCR positivity distribution on arrival
  # Shift distribution departure by phase and testing
  pred_pcr <- function(n_travel=1,
                       before_travel_test=btt,
                       after_arrival_test=att,
                       daily_growth=0){

    pd_depart <- c(rep(1,before_travel_test),head(p_by_day$p_neg,-before_travel_test)) # Probability missed at departure given test before
    phase_shift <- phase_p(1:max_days,growth=daily_growth,maxd=max_days)
    pd_depart_both <- pd_depart*phase_shift
    
    lft_depart <- c(rep(1,before_travel_test),head(l_by_day$p_neg,-before_travel_test))
    lft_depart_both <- lft_depart*phase_shift
    
    # Shift distribution for arrivals
    if(after_arrival_test==0){
      pd_arrive <- p_by_day$median
    }else{
      pd_arrive <- c(tail(p_by_day$median,-after_arrival_test),rep(0,after_arrival_test))
    }

    p_detect_all_raw <- n_travel*pd_depart_both*pd_arrive # Cases detected - relative to time of departure
    p_detect_all <- p_detect_all_raw/sum(p_detect_all_raw)
    
    lft_detect_all <- lft_depart_both*pd_arrive # Cases detected - relative to time of departure
    lft_detect_all <- lft_detect_all/sum(lft_detect_all)
    
    # Normalise PCR detection function at departure
    p_all <- p_by_day$median/sum(p_by_day$median)
    p_dep_raw <- n_travel*(1-pd_depart)*phase_shift
    p_dep <- p_dep_raw/sum(p_dep_raw)
    

    # Output results
    list(p_shift = n_travel*phase_shift,
         out_base = p_all,
         out_raw_depart = pd_depart_both,
         out_raw_arrive = p_detect_all_raw,
         out_depart = p_dep,
         out_arrive = p_detect_all)
  
  }
  
  # - - -
  # Plot curves for growing and declining epidemics
  
  col_red <- rgb(1,0.5,0.5)
  col_green <- rgb(0.8,1,0.8)
  col_green2 <- rgb(0.6,1,0.6)
  
  # Iterate over epidemic types
  growth_label <- c("10% daily decline","Stable epidemic","10% daily growth")
  
  par(mfrow=c(2,3),mgp=c(2,0.7,0),mar = c(3,3,1,1))
  
  grow_seq <- c(-0.1,0,0.1)
  letter_x <- 1
  
  for(kk in 1:length(grow_seq)){
  
    p_all <- pred_pcr(before_travel_test=btt,after_arrival_test=att,daily_growth = grow_seq[kk])
  
    n_travellers <- 100
    plot(-x_labels,n_travellers*p_all$p_shift,type="l",xlim=c(min(-x_labels),att+1),ylim=c(0,n_travellers/10),xaxs="i",yaxs="i",
         col="white",xlab="day of infection relative to arrival",ylab="infected travellers",main=growth_label[kk])
    # lines(-x_labels,n_travellers*p_all$out_raw_depart,col="red")
    # lines(-x_labels,n_travellers*p_all$out_raw_arrive,col="orange")
    # 
    polygon(c(-x_labels,rev(-x_labels)),n_travellers*c(0*p_all$p_shift,rev(p_all$p_shift)),col=col_green,lty=0)
    polygon(c(-x_labels,rev(-x_labels)),n_travellers*c(0*p_all$out_raw_depart,rev(p_all$out_raw_depart)),col=col_red,lty=0)
    polygon(c(-x_labels,rev(-x_labels)),n_travellers*c(0*p_all$out_raw_arrive,rev(p_all$out_raw_arrive)),col=col_green2,lty=0)

    
    lines(-x_labels,n_travellers*p_all$p_shift,col="black",lwd=2)
    lines(-c(btt,btt),c(-10,1e3),lty=2) # minus 1 to set range = 0:2
    lines(c(att,att),c(-10,1e3),lty=2) # minus 1 to set range = 0:2
    lines(c(0,0),c(-10,1e3),col="light grey",lwd=2) # Arrival time
    
    if(kk==2){
      text(x=-28,y=1,labels="missed",adj=0,col="dark red")
      text(x=-17,y=2.6,labels="departure test",adj=0,col="dark green")
      text(x=-5,y=0.7,labels="arrival test",adj=0,col="dark green")
    }

    
    title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
    
  }

  # - - -
  # Reconstruction of curves
  
  n_max <- 150
  x_days <- 1:n_max
  w_s <- seq(1,n_max,7) # Weekly scale
  x_weeks <- x_days[w_s]
  
  col_list1 <- list("dark orange","purple","gold")
  col_list2 <- list("dark orange","dark cyan","purple")
  
  for(ii in 1:3){
    # Define incidence curves
    if(ii==1){
      prop_incidence1 <- 0.2*dnorm(x_days,mean=50,sd=22) # Proportion infected per day
    }
    
    if(ii==2){
      a1 <- 0.1
      peak_sim <- 0.005; peak_sim2 <- 0.002; cut_sim <- 50; cut_sim2 <- 120
      prop_incidence1<- rep(0,n_max)
      prop_incidence1[1:cut_sim] <- peak_sim*exp(a1*(1:cut_sim))/max(exp(a1*(1:cut_sim)))
      prop_incidence1[(cut_sim+1):n_max] <- peak_sim*exp(-a1*((cut_sim+1):n_max))/max(exp(-a1*((cut_sim+1):n_max)))
    }
    
    if(ii==3){
      peak_sim <- 0.003; peak_sim2 <- 0.001; cut_sim <- 50; cut_sim2 <- 120
      prop_incidence1<- rep(0,n_max)
      prop_incidence1[(cut_sim+1):n_max] <- peak_sim
    }
  
    
    btt1 <- 1; aat1 <- 0; # Scenario 1
    btt2 <- 3; aat2 <- 4; # Scenario 2
    
    # Simulate and recover dynamics
    stochastic_arrival(prop_incidence1,btt = btt2,att = aat2,n_max,n_travel_vol = 5000,test_type="PCR",
                       test_type_infer="PCR",col_list_in = col_list2[[ii]], add_gam = add_gam_in,x_weeks,w_s)
    
    title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  }
  
  dev.copy(pdf,paste("plots/Fig2_phase_bias.pdf",sep=""),width=8,height=5)
  dev.off()
  
  # - - - 
  # Plot Figure S2 with more scenarios:
  
  par(mfcol=c(3,3),mgp=c(2,0.7,0),mar = c(3,3,1,1))
  letter_x <- 1
  
  title_label <- c("2d pre + 4d post","2d + 0d","2d + 4d (80% sensitivity)")
  
  for(ii in 1:3){
    for(kk in 1:2){
    # Define incidence curves
    prop_incidence1 <- 0.2*dnorm(x_days,mean=50,sd=22) # Proportion infected per day
    
    if(ii==1 | 3){btt2 <- 3; aat2 <- 4} # Scenario 2
    if(ii==2){btt2 <- 1; aat2 <-0}; # Scenario 1
    if(ii==3){test_inf <- 0.8}else{test_inf <- "PCR"} # Test sensitivity used at arrival/departure
      
    if(kk==1){travel_n=1000}
    if(kk==2){travel_n=5000}
    
    # Simulate and recover dynamics
    stochastic_arrival(prop_incidence1,btt = btt2,att = aat2,n_max,n_travel_vol = travel_n,test_type=test_inf,
                       test_type_infer="PCR",col_list_in = col_list2[[ii]], add_gam = add_gam_in,x_weeks,w_s)
    text(x=100,y=3.5,adj=0,labels=paste0("n=",travel_n))
    
    title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
    
    if(kk==1){title(main=title_label[ii],adj=0.5)}
    
    }
  }
  
  dev.copy(pdf,paste("plots/FigS1_phase_bias.pdf",sep=""),width=8,height=4)
  dev.off()
  
  
}



# Figure 3: Summary incidence --------------------------------------------------

figure_basic_data <- function() {
  
  # Set up values
  ylab_v <- 6e3; y_arrow <- 3e3
  xlim_v <- c(as.Date("2020-07-01"),as.Date("2022-04-01")) #c(min(travel_incidence$dates),max(travel_incidence$dates))
  grid_v <- 5 #seq(as.Date("2020-01-01"),as.Date("2023-01-01"),"quarter")
  letter_x <- 1
  
  layout(matrix(c(1,2,3,1,2,3,4,5,6), ncol=3))
  par(mgp=c(2,0.7,0),mar = c(2,3,1,3))
  
  # - - -
  # Testing data
  plot(travel_incidence_n$dates,rowSums(travel_incidence_n[,c("NON","OUI")]),
       xlim=xlim_v,ylim=c(0,7000),xaxs="i",yaxs="i",col="white",xlab="2020/21",ylab="weekly tests",type="l",lwd=2)
  
  grid_year()

  polygon(c(covd1[1],covd1_q[2],covd1_q[2],covd1[1]),c(0,0,1e5,1e5),col=rgb(0,0,1,0.1),lty=0)
  polygon(c(labd[1],covd2[2],covd2[2],labd[1]),c(0,0,1e5,1e5),col=rgb(0,1,1,0.1),lty=0)
  
  text(x=covd1[1],y=ylab_v,labels="Day 4 test",adj=0,col="dark blue")
  text(x=labd[1],y=ylab_v,labels="Day 0 test",adj=0,col="dark cyan")
  
  lines(travel_incidence_n$dates,rowSums(travel_incidence_n[,c("NON","OUI")]),lwd=2) #COV-CHECK
  
  par(new=TRUE)
  plot(travel_incidence_n$dates,travel_incidence_n$OUI,col="dark orange",xaxt="n",
       yaxt="n",ylab="",ylim=c(0,250),type="l",lwd=2,xlim=xlim_v,xaxs="i",yaxs="i")
  
  axis(4,col="dark orange",col.axis="dark orange")
  mtext("positive tests", side=4, line=1.5,col="dark orange",cex=0.7) # Label for 2nd axis
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - -
  # Plot positives
  plot(ymd(travel_incidence_n$dates),100*travel_incidence_n$OUI/rowSums(travel_incidence_n[,c("OUI","NON")]),
       xlim=ymd(xlim_v),ylim=c(0,5),col="white",xaxs="i",yaxs="i",xlab="2020/21",ylab="% positive")
  grid_year()
  
  plot_CI(travel_incidence_n$dates,travel_incidence_n$OUI,rowSums(travel_incidence_n[,c("NON","OUI")]))
  
  # Define ranges for plotting
  range1a <- 1:42 # Protocol with day 4 test
  range2a <- 43:91  # Protocol with day 0 test
  out_fp1 <- plot_GAM(travel_incidence_n$dates[range1a],travel_incidence_n$OUI[range1a],rowSums(travel_incidence_n[range1a,c("NON","OUI")]))
  out_fp2 <- plot_GAM(travel_incidence_n$dates[range2a],travel_incidence_n$OUI[range2a],rowSums(travel_incidence_n[range2a,c("NON","OUI")]))
  
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - -
  # Case data
  plot(fp_cases$date,fp_cases$new_cases,col="white",xaxs="i",yaxs="i",type="l",ylim=c(0,8000),ylab="weekly cases in FP", xlim=xlim_v,lwd=2)
  grid_year()
  lines(fp_cases$date,fp_cases$new_cases,lwd=2)
  
  case_1_wt <- as.Date("2020-03-09");   case_1_cc <- as.Date("2020-07-27")
  case_1_alpha <- as.Date("2020-12-28")
  case_1_delta <- as.Date("2021-06-08")
  case_1_ba1 <- as.Date("2021-12-12")
  case_1_ba2 <- as.Date("2022-01-14")
  
  arrows(case_1_cc,y_arrow,case_1_cc,0,angle=30,length=0.07,col="red"); text(x=case_1_cc,y=1.2*y_arrow,labels="Wildtype",col="red")
  arrows(case_1_alpha,y_arrow,case_1_alpha,0,angle=30,length=0.07,col="red"); text(x=case_1_alpha,y=1.2*y_arrow,labels="Alpha",col="red")
  arrows(case_1_delta,y_arrow,case_1_delta,0,angle=30,length=0.07,col="red"); text(x=case_1_delta,y=1.2*y_arrow,labels="Delta",col="red")
  arrows(case_1_ba1,y_arrow,case_1_ba1,0,angle=30,length=0.07,col="red"); text(x=case_1_ba1,y=1.2*y_arrow,labels="Omicron",col="red")
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - - 
  # Ct values
  par(mgp=c(2,0.7,0),mar = c(3,3,1,1))
  break_x <- seq(10,40,5)
  hcw_ct_p <- hcw_ct %>% filter(ct < 37 & x_date>=0) # Select Ct positives
  
  hist(travel_test_ct$ct_value,breaks=break_x,ylab="probability",prob=T,ylim=c(0,0.08),
       xlab="Ct value",main="",col=NULL,xaxs="i",yaxs="i",lwd=2,border="orange") # rgb(0,0.8,1,0.6)
  
  hist(hcw_ct_p$ct,breaks=break_x,add=T,prob=T,col=NULL,border="dark cyan",lwd=2) # rgb(1,0.5,0,0.5)
  
  text(x=12,y=0.06,labels="HCW",adj=0,col="dark cyan")
  text(x=12,y=0.07,labels="Travellers",adj=0,col="orange")
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - - 
  # Detection probability
  growth_vals <- c(-0.1,0.1)
  growth_label <- c("10% daily decline","10% daily growth")
  test_days <- c(0:8)
  
  for(ii in 2:1){
    
    # Calculate probability of arrival detectionn
    y_vals <- sapply(test_days,function(x){
      out_ii <- delay_test(n_depart_x=31,before_travel_test=3,after_arrival_test=x,
                           daily_growth = growth_vals[ii])
      out_ii_ant <- delay_test(n_depart_x=31,before_travel_test=2,after_arrival_test=x,
                               daily_growth = growth_vals[ii],
                               dep_test="Antigen")
      
      early_cut1 <- 10 + x
      early_cut2 <- 5 + x
      
      n2 <- out_ii$n_arrive_detect;  n2_ant <- out_ii_ant$n_arrive_detect
      ntot_arr <- out_ii$n_total; ntot_arr_ant <- out_ii_ant$n_total
      ntot_dep <- out_ii$n_total_dep; ntot_dep_ant <- out_ii_ant$n_total_dep
      c(100*sum(n2[1:early_cut1])/sum(ntot_arr[1:early_cut1]),
        100*sum(n2[1:early_cut2])/sum(ntot_arr[1:early_cut2]),
        100*sum(n2_ant[1:early_cut1])/sum(ntot_arr_ant[1:early_cut1]),
        100*sum(n2_ant[1:early_cut2])/sum(ntot_arr_ant[1:early_cut2])
      )
    })
    
    plot(test_days,y_vals[1,],col="dark green",lty=2,
         xlab="day of post-arrival test",ylab="% infections detected",ylim=c(0,80),
         main = "")  # Number detected within 5d (arrival)
    points(test_days,y_vals[2,],col="dark green",pch=3) # Number detected within 10d (arrival)
    points(test_days,y_vals[3,],col="green",pch=1) # Number detected within 5d (depart)
    points(test_days,y_vals[4,],col="green",pch=3) # Number detected within 10d (depart)
    title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
    text(x=4,y=77,labels=growth_label[ii],adj=0.5,font=2)
    
    if(ii==2){
      sze <- 0.8
      text(x=2,y=28,labels="48hr dep. PCR",adj=0,col="dark green",cex=sze)
      text(x=2,y=20,labels="24hr dep. antigen",adj=0,col="green",cex=sze)
      text(x=2,y=12,labels="+ within 5d of infection",adj=0,col="dark grey",cex=sze)
      text(x=2,y=4,labels="o within 10d of infection",adj=0,col="dark grey",cex=sze)
    }
    
  }
  
  
  dev.copy(pdf,paste("plots/Fig3_data.pdf",sep=""),width=7,height=6)
  dev.off()
  
}



# Figure 4: Reconstruct epidemics -----------------------------------------------

figure_reconstruct_epidemics <- function(test_type1="PCR",test_type2="PCR",btt=3){
  
  # DEBUG test_type="PCR"; btt=3

  # Pop size
  us_pop <- 329.5e6
  fr_pop <- 67.39e6

  # Calculate positives by country of origin
  match_date <- match(travel_incidence_n$dates,country_incidence$dates)
  adj_case_counts <- country_incidence[match_date,]
  
  pos_counts_fr <- adj_case_counts$France
  pos_counts_us <- adj_case_counts$USA

  # Estimate recent growth rates in prevalence data
  
  # Calculate moving average for France and US
  fr_cases_ma <- rollmean(fr_cases$new_cases,k=7,fill=c(NA,NA,NA))
  us_cases_ma <- rollmean(us_cases$new_cases,k=7,fill=c(NA,NA,NA))

  # First phase - testing on day 4
  # scale_depart_fr1 <- estimate_vals(10,daily_growth=0,before_travel_test=3,dep_test = "PCR")$ratio_total
  # scale_depart_us1 <- estimate_vals(10,daily_growth=0,before_travel_test=3,dep_test = "PCR")$ratio_total
  # 
  # # Second phrase - testing on day 0
  # scale_depart_fr2 <- estimate_vals(10,daily_growth=0,before_travel_test=btt,after_arrival_test=0,dep_test = test_type)$ratio_total
  # scale_depart_us2 <- estimate_vals(10,daily_growth=0,before_travel_test=btt,after_arrival_test=0,dep_test = test_type)$ratio_total
  # 
  # # Combine the two different scalings
  # scale_depart_fr <- scale_depart_fr2; scale_depart_fr[range1] <- scale_depart_fr1[range1]
  # scale_depart_us <- scale_depart_us2; scale_depart_fr[range1] <- scale_depart_us1[range1]
  
  # Calculate positive values by test phase
  # fr_pos_vals1 <- scale_depart_fr1*pos_counts_fr
  # us_pos_vals1 <- scale_depart_us1*pos_counts_us
  # 
  # fr_pos_vals2 <- scale_depart_fr2*pos_counts_fr
  # us_pos_vals2 <- scale_depart_us2*pos_counts_us
  # 
  # fr_pos_vals_lab <- scale_depart_fr2*pos_counts_fr_lab
  # us_pos_vals_lab <- scale_depart_fr2*pos_counts_us_lab
  # 
  # - - -
  # Plot incidence
  xlim_v <- c(min(travel_incidence_n$dates),max(travel_incidence_n$dates))
  xlim_v2 <- ymd(xlim_v)     
  

  par(mfrow=c(4,2),mgp=c(2,0.7,0),mar = c(2,3,1,3))
  
  letter_x <- 1
  
  # France arrivals plot
  plot(travel_incidence_n$dates[range1],tests_fr_s1[range1],main="France",
       col="dark orange",ylab="FP arrivals tested weekly",ylim=c(0,4000),type="l",lwd=2,xlim=xlim_v2,xaxs="i",yaxs="i")
  lines(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2]),c(tests_fr_lab,tests_fr_s2[range2]),col="gold",lwd=2)

  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # US arrivals plot
  plot(travel_incidence_n$dates[range1],tests_us_s1[range1],main="USA",
       col="dark orange",ylab="FP arrivals tested weekly",ylim=c(0,4000),type="l",lwd=2,xlim=xlim_v2,xaxs="i",yaxs="i")
  lines(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2]),c(tests_us_lab,tests_us_s2[range2]),col="gold",lwd=2)

  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - -
  # France values
  pos_counts_fr[is.na(pos_counts_fr)] <- 0
  pos_counts_fr[is.na(pos_counts_fr)] <- 0
  
  plot(travel_incidence_n$dates,round(pos_counts_fr),col="white",ylim=c(0,8),xaxs="i",yaxs="i",ylab="departure prevalence (%)")
  grid_year()
  
  data_fr1 <- data.frame(pos=round(pos_counts_fr)[range1],num=round(tests_fr_s1)[range1])
  data_fr2 <- data.frame(pos=round(c(pos_counts_fr_lab,pos_counts_fr[range2])),num=c(tests_fr_lab,tests_fr_s2[range2]))
  
  # Include shift to get estimate at departure
  out_fr1r <- bootstrap_est(travel_incidence_n$dates[range1]-7,data_fr1,test_type=test_type1)
  out_fr2r <- bootstrap_est(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2])-3,data_fr2,before_travel_test=btt,after_arrival_test=0,test_type=test_type2)
  
  # Observed prevalence
  # plot_CI(travel_incidence_n$dates[range1],round(pos_counts_fr)[range1], round(tests_fr_s1)[range1])
  # plot_CI(travel_incidence_n$dates[range2],round(pos_counts_fr)[range2], round(tests_fr_s2)[range2])
  # plot_CI(travel_incidence_n$dates[range_lab],round(pos_counts_fr_lab), tests_fr_lab)
  # 
  #plot_CI_def(,  apply(data_fr1,1,function(x){prev_estimate(x,3,4)$prev_est}) )



  # # 
  # # Smoothing observed prevalence -  merge lab and COV-CHECK 2 data (both day 0)
  # out_fr1r <- plot_GAM(travel_incidence_n$dates[range1],pos_counts_fr[range1],tests_fr_s1[range1])
  # out_fr2r <- plot_GAM(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2]),
  #                      c(pos_counts_fr_lab,pos_counts_fr[range2]),
  #                      c(tests_fr_lab,tests_fr_s2[range2]))
  # #out_fr_labr <- plot_GAM(,,)

  
  # Estimated prevalence
  # 
  # plot_polygon(out_fr1r$pred_date,out_fr1r$pred_med,out_fr1r$pred_CI1,out_fr1r$pred_CI2,scale_depart_fr1[1],col1=colA,colf=colB)
  # plot_polygon(out_fr2r$pred_date,out_fr2r$pred_med,out_fr2r$pred_CI1,out_fr2r$pred_CI2,scale_depart_fr2[1],col1=colA,colf=colB)

  colA <- "red"; colB <-rgb(0,0,1,0.2) # rgb(1,0,0,0.1) colf = 
  
  x_text <- min(out_fr1r$pred_date)+20
  text(x=x_text+5,y=7,labels="  Estimated prevalence",col="black",adj=0)
  text(x=x_text,y=6,labels="- Smoothed GAM estimate",col="blue",adj=0)
  points(x=x_text,y=7,pch=19,cex=0.8)
  
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - -
  # US values
  pos_counts_us[is.na(pos_counts_us)] <- 0
  pos_counts_us[is.na(pos_counts_us)] <- 0
  
  plot(travel_incidence_n$dates,round(pos_counts_fr),col="white",ylim=c(0,8),xaxs="i",yaxs="i",ylab="departure prevalence (%)")
  grid_year()
  
  # Observed prevalence
  
  data_us1 <- data.frame(pos=round(pos_counts_us)[range1],num=round(tests_us_s1)[range1])
  data_us2 <- data.frame(pos=round(c(pos_counts_us_lab,pos_counts_us[range2])),num=c(tests_us_lab,tests_us_s2[range2]))
  
  out_us1r <- bootstrap_est(travel_incidence_n$dates[range1]-7,data_us1,test_type=test_type1)
  out_us2r <- bootstrap_est(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2])-3,data_us2,before_travel_test=btt,after_arrival_test=0,test_type=test_type2)
  
  
  # plot_CI(travel_incidence_n$dates[range2],round(pos_counts_us)[range2], round(tests_us_s2)[range2])
  # plot_CI(travel_incidence_n$dates[range1],round(pos_counts_us)[range1], round(tests_us_s1)[range1])
  # plot_CI(travel_incidence_n$dates[range_lab],round(pos_counts_us_lab), tests_us_lab)
  # 
  # # Smoothing observed prevalence -  merge lab and COV-CHECK 2 data (both day 0)
  # out_us1r <- plot_GAM(travel_incidence_n$dates[range1],pos_counts_us[range1],tests_us_s1[range1])
  # out_us2r <- plot_GAM(c(travel_incidence_n$dates[range_lab],travel_incidence_n$dates[range2]),
  #                      c(round(pos_counts_us_lab),pos_counts_us[range2]),
  #                      c(tests_us_lab,tests_us_s2[range2]))
  # #out_us2r_only <- plot_GAM(c(travel_incidence_n$dates[range2]),c(pos_counts_us[range2]),c(tests_us_s2[range2])) # COV-CHECK 2 for antibody comparison
  # 
  # # Estimated prevalence
  # plot_polygon(out_us1r$pred_date,out_us1r$pred_med,out_us1r$pred_CI1,out_us1r$pred_CI2,scale_depart_us1[1],col1=colA,colf=colB)
  # plot_polygon(out_us2r$pred_date,out_us2r$pred_med,out_us2r$pred_CI1,out_us2r$pred_CI2,scale_depart_us2[1],col1=colA,colf=colB)

  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # - - -
  # Compare to France/US case data
  
  # France cases plot
  plot(fr_cases$date,1e5*fr_cases_ma/fr_pop,ylab="daily domestic cases (per 100k)",type="l",
       lwd=2,xlim=xlim_v2,ylim=c(0,600),xaxs="i",yaxs="i")
  grid_year()

  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # US cases plot
  ymax <- 300
  plot(us_cases$date,1e5*us_cases_ma/us_pop,ylab="daily domestic cases (per 100k)",
       type="l",lwd=2,xlim=xlim_v2,ylim=c(0,300),xaxs="i",yaxs="i")
  grid_year()
  
  par(new=TRUE)
  waste_y <- max(wastewater_nat$effective_concentration_rolling_average)/(max(1e5*us_cases_ma/us_pop,na.rm=T)/ymax)
  
  plot(wastewater_nat$sampling_week,wastewater_nat$effective_concentration_rolling_average,col="dark orange",xaxt="n",
       yaxt="n",ylab="",ylim=c(0,waste_y),type="l",lwd=2,xlim=xlim_v,xaxs="i",yaxs="i")
  
  axis(4,col="dark orange",col.axis="dark orange")
  mtext("SARS-CoV-2 concentration", side=4, line=1.5,col="dark orange",cex=0.7) # Label for 2nd axis


  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  #par(las=1)
  
  # - - -
  # Match to France seroprevalence data 1

  shift_s <- 17 # number of days to seroconvert
  mean_pcr <- round(sum((0:30)*p_by_day$median/sum(p_by_day$median))) # mean PCR duration

  y_cinc <- cumul_CI(out_fr1r$sim_inc,cumulative=T)
  x_date_p <- out_fr1r$sim_date - mean_pcr + shift_s  # Shift to get incidence (mean PCR duration = 9 days)
  
  # Match to seroprevalence timings
  date1 <- as.Date("2020-10-11")
  date2 <- as.Date("2021-02-14")
  
  d1 <- which.min(abs(x_date_p - as.Date("2020-10-11")))
  d2 <- which.min(abs(x_date_p - as.Date("2021-02-14")))
  
  sero1 <- 9; sero2 <- 13.2
  
  d3 <- as.Date("2020-11-01"); sero3a <- 8.5
  d4 <- as.Date("2021-01-31"); sero4a <- 15
  
  plot(x_date_p,y_cinc$median,
       xlim=xlim_v2,type="l",xaxs="i",yaxs="i",ylim=c(0,62),col="white",ylab="cumulative % infected")
  grid_year()
  points(c(date1,date2),c(sero1,sero2),pch=19)
  lines(c(date1,date1),c(6.9,11)); lines(c(date2,date2),c(10.8,15.6))
  points(c(d3),sero3a,pch=17)
  points(c(d4),sero4a,pch=15)

  # Plot estimates
  polygon(c(x_date_p[d1:d2],rev(x_date_p[d1:d2])),
          c((y_cinc$lower[d1:d2]-min(y_cinc$lower[d1:d2])+sero1),
            rev((y_cinc$upper[d1:d2]-min(y_cinc$upper[d1:d2])+sero1)) ),
          col=colB,lty=0)
  lines(x_date_p[d1:d2],(y_cinc$median[d1:d2]-min(y_cinc$median[d1:d2])+sero1),col="blue")
  
  # Plot cases and comparison
  cumsum_fr <- cumsum(coalesce(fr_cases_ma, 0)) + fr_cases_ma*0
  lines(fr_cases$date,100*cumsum_fr/fr_pop,col="grey")

  x_text <- min(out_fr1r$pred_date)+20
  text(x=x_text+5,y=50,labels="  Domestic seroprevalence surveys",col="black",adj=0)
  text(x=x_text,y=40,labels="- Estimated from FP arrivals",col="blue",adj=0)
  text(x=x_text,y=30,labels="- Cumulative domestic cases",col="dark gray",adj=0)
  points(x=x_text,y=50,pch=19,cex=0.8)
  
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # Print peak values for France and US in each wave:

  c.text(max(out_fr1r$pred_med),
         max(out_fr1r$pred_CI1),
         max(out_fr1r$pred_CI2)) %>% print()
  
  c.text(max(out_fr2r$pred_med),
         max(out_fr2r$pred_CI1),
         max(out_fr2r$pred_CI2)) %>% print()
  
  c.text(max(out_us1r$pred_med),
         max(out_us1r$pred_CI1),
         max(out_us1r$pred_CI2)) %>% print()
  
  c.text(max(out_us2r$pred_med),
         max(out_us2r$pred_CI1),
         max(out_us2r$pred_CI2)) %>% print()
  
  # - - -
  # Match to US seroprevalence data 1
  y_cinc <- cumul_CI(out_us1r$sim_inc)
  x_date_p <- out_us1r$sim_date - mean_pcr + shift_s # Shift to get incidence (mean PCR duration = 9 days)

  d1 <- which.min(abs(x_date_p - as.Date("2020-07-27")))
  d2 <- length(x_date_p)
  
  sero1 <- 100*us_national_sero$mid[1] # Initial US national serosurvey estimate
  
  # Match to US seroprevalence data 2
  y_cinc2 <- cumul_CI(out_us2r$sim_inc)
  x_date_p2 <- out_us2r$sim_date - mean_pcr # Shift to get incidence (mean PCR duration = 9 days)

  d1a <- which.min(abs(x_date_p2 - as.Date("2021-05-06")))
  d2a <- which.min(abs(x_date_p2 - as.Date("2022-02-11")))
  
  sero1a <- us_national_sero |> filter(date_mid=="2021-05-06")
  sero1a <- 100*sero1a$mid[1] # US national serosurvey estimate at start of round 2 travel testing

  
  plot(x_date_p[d1:d2],(y_cinc$median[d1:d2]-min(y_cinc$median[d1:d2]))+sero1,
       xlim=xlim_v2,type="l",xaxs="i",yaxs="i",ylim=c(0,62),col="white",ylab="cumulative % infected")
  grid_year()

  # Plot seroprevalence data over time in US (national labs)
  
  for(ii in 1:nrow(us_national_sero)){
    s_pick <- us_national_sero[ii,]
    points(s_pick$date_mid,100*s_pick$mid,pch=19,cex=0.7)
    lines(c(s_pick$date_mid,s_pick$date_mid),100*c(s_pick$lower,s_pick$upper))
  }

  # Plot estimates for first wave
  polygon(c(x_date_p[d1:d2],rev(x_date_p[d1:d2])),
          c((y_cinc$lower[d1:d2]-min(y_cinc$lower[d1:d2])+sero1),
            rev((y_cinc$upper[d1:d2]-min(y_cinc$upper[d1:d2])+sero1)) ),
          col=colB,lty=0)
  lines(x_date_p[d1:d2],(y_cinc$median[d1:d2]-min(y_cinc$median[d1:d2])+sero1),col="blue")
  
  # Plot estimates for second wave
  polygon(c(x_date_p2[d1a:d2a],rev(x_date_p2[d1a:d2a])),
          c((y_cinc2$lower[d1a:d2a]-min(y_cinc2$lower[d1a:d2a])+sero1a),
            rev((y_cinc2$upper[d1a:d2a]-min(y_cinc2$upper[d1a:d2a]))+sero1a) ),
          col=colB,lty=0)
  lines(x_date_p2[d1a:d2a],(y_cinc2$median[d1a:d2a]-min(y_cinc2$median[d1a:d2a])+sero1a),col="blue")

  
  # Plot cases and comparison
  cumsum_us <- cumsum(coalesce(us_cases_ma, 0)) + us_cases_ma*0
  lines(us_cases$date,100*cumsum_us/us_pop,col="gray")
  
  
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  
  
  dev.copy(pdf,paste("plots/Fig4_estimates_",test_type1,"_",test_type2,"_",btt,".pdf",sep=""),width=6,height=7)
  dev.off()

}

# References for seroprevalence:

# US national frequent data given in data_load.R

# Prevalence France 1
# https://www.santepubliquefrance.fr/etudes-et-enquetes/covid-19-etudes-pour-suivre-la-part-de-la-population-infectee-par-le-sars-cov-2-en-france
# COVID-19 une étude pour connaître la part de la population infectée par le coronavirus en France
# weeks 41 (05-11 October 2020) and 6 (8-14 February 2021)
# La séroprévalence nationale est passée de 9,0 % [6,9 -11,0] 
# 13,2 % [10,8-15,6] au cours de la seconde vague

# Prevalence France 2
# 4 % de la population a développé des anticorps contre le SARS-CoV-2 entre mai et novembre 2020
# 8.5% in November
# ER1202_0

# Prevalence France 3
# https://www.medrxiv.org/content/10.1101/2022.07.29.22278190v1.full
# 15% end of Jan 2021


# Figure S2: Compare simulated incidence and prevalence estimated -----------------------------------------------

figure_explore_cumulative_incidence <- function(){
  
  # Simulate infection dynamics
  
  xx <- 0:150
  data_infections <- 0.1*(2+sin(8*pi*(xx-20)/365))
  n_inf <- length(data_infections)
  plot(data_infections)
  
  # Define transition matrix to construct prevalence
  f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
  n_pcr_days <- length(p_by_day$median)
  
  for(ii in 1:n_inf){
    i_max <- min(ii+n_pcr_days-1,n_inf)
    j_max <- min(n_inf-ii+1,n_pcr_days)
    
    f_matrix[ii:i_max,ii] <- p_by_day$median[1:j_max]
    
  }
  
  invert_f <- ginv(f_matrix)
  
  estimate_incidence <- function(prevalence){
    
    # Define transition matrix - 
    n_inf <- length(prevalence)
    f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
    n_pcr_days <- length(p_by_day$median)
    
    for(ii in 1:n_inf){
      i_max <- min(ii+n_pcr_days-1,n_inf)
      j_max <- min(n_inf-ii+1,n_pcr_days)
      
      f_matrix[ii:i_max,ii] <- p_by_day$median[1:j_max]
      
    }
    
    invert_f <- ginv(f_matrix)
    output <- invert_f %*% prevalence
    output
    
  }
  
  data_prev <- f_matrix %*% data_infections
  
  par(mfrow=c(2,2),mgp=c(2,0.7,0),mar = c(3,3,1,1))
  
  # Incidence
  plot(data_infections,xaxs="i",yaxs="i",ylab="daily incidence (%)",ylim=c(0,0.35),xlab="days")
  
  letter_x <- 1
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # Prevalence
  plot(data_prev,type="l",lwd=2,xaxs="i",yaxs="i",ylab="daily prevalence (%)",ylim=c(0,3),xlab="days")
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # Reconstruction of incidence
  plot(data_infections,xaxs="i",yaxs="i",ylab="daily incidence (%)",ylim=c(0,0.35),xlab="days")
  
  inc1 <- estimate_incidence(data_prev)
  inc2 <- prev_to_inc(data_prev,weekly_adjust=1) %>% tail(.,-9)
  lines(inc1,col="blue",lwd=2)
  lines(inc2,col="red",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  # Reconstruction of cumulative incidence
  plot(cumsum(data_infections),xaxs="i",yaxs="i",ylab="cumulative incidence (%)",ylim=c(0,50),xlab="days")
  
  lines(cumsum(inc1),col="blue",lwd=2)
  lines(cumsum(inc2),col="red",lwd=2)
  title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1
  
  dev.copy(pdf,paste("plots/FigS3_incidence_estimates.pdf",sep=""),width=8,height=6)
  dev.off()
  
}


#  Plot pooled testing scenarios ------------------------------------------

figure_pool_test <- function(){
  
  # Set up example numbers
  n_plane <- 2 # Number of planes
  t_passengers <- 200 # Passengers tested per plane
  pool_n_plane <- 10 # Pools per plane
  pool_daily <- pool_n_plane*n_plane # Pools per day
  
  pool_size <- t_passengers/pool_n_plane # Pool size
  
  
  # Plot data and CI
  plot_CI <- function(xx,med,CI1,CI2,colA="black") {
    
    for(ii in 1:length(med)){
      points(xx[ii],med[ii],col=colA,pch=19); lines(c(xx[ii],xx[ii]),c(CI1[ii],CI2[ii]),col=colA)
    }
    
  }
  
  # Calculate MLE and profile likelihood
  calculate_MLE <- function(x_pool, # Positive pools
                            n_pool, # Total pools
                            size # Size of a pool
  ){
    
    # DEBUG: x_pool = 20; n_pool = pool_daily; size = pool_size
    
    xx <- seq(0,1,0.0001)
    
    # Probability at least one positive in pool
    p_pool <- 1-(1-xx)^size
    
    # Binomial likelihood
    b_lik <- dbinom(x_pool,n_pool,prob = p_pool,log=T)
    
    mle_x <- max(100*xx[which(b_lik==max(b_lik))])
    ci1_x <- 100*xx[which(b_lik>(max(b_lik)-1.92))] %>% min()
    ci2_x <- 100*xx[which(b_lik>(max(b_lik)-1.92))] %>% max()
    
    c_vec <- c(mle_x,ci1_x,ci2_x)
    c_vec
    #ci_text <- paste0(mle_x,"% (",ci1_x,"-",ci2_x,"%)")
    
  }
  
  
  
  # Plot estimates
  
  par(mfcol=c(2,2),mar=c(3,3,1,1),las=0,mgp=c(2,0.7,0))
  passenger_list <- c(100,500,100,500)
  pool_list <- c(5,20)
  
  col_list <- c("orange","blue","purple")
  
  for(kk in 1:4){
    
    t_passengers <- passenger_list[kk] # Passengers tested
    
    for(jj in 1:2){
      
      pool_size <- pool_list[jj] # Pool size
      pool_daily <- t_passengers/pool_size # Pools tested
      total_screened <- pool_daily*pool_size # Total screened
      
      pool_positives <- c(0:ceiling(1*pool_daily))
      
      pool_positives_est <- sapply(pool_positives,function(x){calculate_MLE(x,pool_daily,pool_size)})
      pool_out <- data.frame(positive_pools=pool_positives,estimated_prevalence=t(pool_positives_est))
      
      if(kk==1 | kk==2){xmax=100;ymax=100;headerK=paste0(total_screened," passengers")}
      if(kk==3 | kk==4){xmax=40;ymax=20;headerK=""}
      
      if(jj==1){
        plot(100*pool_out$positive_pools/pool_daily,pool_out$estimated_prevalence.1,
             ylab="estimated prevalence (%)",
             xlab="% pools positive",
             col="white",ylim=c(0,ymax),xlim=c(0,xmax),
             main=headerK)
      }
      if(kk==1 | kk==2){
        polygon(c(0,40,40,0),c(0,0,20,20),border="light gray")
        text(x=0,y=80,labels="Pool size:",adj=0)
        text(x=5,y=80-10*jj,labels=pool_size,adj=0,col=col_list[[jj]])
        text(x=0,y=25,labels=LETTERS[kk+2],adj=0,col="gray",cex=0.7)
      }
      
      
      title(main=LETTERS[kk],adj=0)
      plot_CI(100*pool_out$positive_pools/pool_daily,pool_out$estimated_prevalence.1,
              pool_out$estimated_prevalence.2,pool_out$estimated_prevalence.3,col=col_list[[jj]])
      
      
    }
    
    lines(c(0,100),c(0,100),lty=2)
    
  }
  
  dev.copy(pdf,paste("plots/Fig5.pdf",sep=""),width=8,height=6)
  dev.off()
  
}
  

