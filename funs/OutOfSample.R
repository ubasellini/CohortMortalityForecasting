

################################################################
## Functions to compute out-of-sample outcomes
################################################################

OutOfSample_function <- function(ages,cohorts,h.out,
                                 cMX.obs,cMX.sim,
                                 alpha=80){
  
  ## dimensions
  m <- length(ages)
  nc <- length(cohorts)
  
  ## cohorts and rates for out-of-sample exercise
  cohorts.out <- cohorts[1]:cohorts[nc-h.out]
  nc.out <- length(cohorts.out)
  cMX.out <- cMX[,cohorts%in%cohorts.out]
  
  ## years of out-of-sample region
  last.t.out <- max(cohorts)-h.out
  t.out <- (last.t.out+1):max(cohorts)
  
  ## find out-of-sample region
  df.data <- tibble(ages=rep(ages,nc),cohorts=rep(cohorts,each=m),
                    rates=c(cMX.obs),exposures=c(cE)) %>% 
    mutate(years=cohorts+ages,
           oos=case_when(
             years %in% t.out ~ 1,
             TRUE ~ 0),
           oos=ifelse(cohorts>max(cohorts) - h.out,0,oos))
  
  ## alpha level 
  alpha.lev <- alpha/100
  
  ## median, lower and upper ETA
  ETA.sim <- log(cMX.sim)
  ETA.sim.upp <- ETA.sim.low <- ETA.sim.med <- matrix(NA,m,nc.out)
  i <- 1
  for (i in 1:m){
    ETA.sim.med[i,] <- apply(ETA.sim[i,,], 1, median, na.rm = T)
    ETA.sim.upp[i,] <- apply(ETA.sim[i,,], 1, quantile, prob=1- (1-alpha.lev)/2, na.rm = T)
    ETA.sim.low[i,] <- apply(ETA.sim[i,,], 1, quantile, prob=(1-alpha.lev)/2, na.rm = T)
  }
  
  ## variance of forecast (for DSS)
  VAR.sim <- matrix(NA,m,nc.out)
  i <- j <- 1
  for (i in 1:m){
    for (j in 1:nc.out){
      VAR.sim[i,j] <- var(ETA.sim[i,j,])
    }
  }
  
  ## data-frame with simulation results
  df.out <- tibble(ages=rep(ages,nc.out),cohorts=rep(cohorts.out,each=m),
                   cMX_med=exp(c(ETA.sim.med)),
                   cMX_low=exp(c(ETA.sim.low)),
                   cMX_upp=exp(c(ETA.sim.upp)),
                   var.eta=c(VAR.sim))
  
  ## combine data with forecasts
  df.fin <- df.data %>% 
    left_join(df.out,by = join_by(ages, cohorts))
  
  ## compute RMSE log-rates
  df.rmse <- df.fin %>% 
    filter(oos==1) %>%
    filter(rates!=0) %>%  ## remove rates == 0
    mutate(sq.diff=(log(rates)-log(cMX_med))^2)
  rmse <- sqrt(sum(df.rmse$sq.diff)/nrow(df.rmse))
  
  ## compute CPD log-rates
  df.cpd <- df.fin %>% 
    filter(oos==1) %>% 
    filter(rates!=0) %>%  ## remove rates == 0
    mutate(coverage=case_when(
      log(rates) > log(cMX_low) & log(rates) < log(cMX_upp) ~ 1,
      TRUE  ~ 0))
  cpd <- 100*abs(alpha.lev - sum(df.cpd$coverage)/nrow(df.cpd))
  
  ## compute DSS log-rates
  df.dss <- df.fin %>% 
    filter(oos==1) %>%
    filter(rates!=0) %>%  ## remove rates == 0
    mutate(dss=((log(rates)-log(cMX_med))^2)/var.eta + 2*log(sqrt(var.eta)))
  dss <- sum(df.dss$dss)/nrow(df.dss)
  
  ## output
  out <- list(rmse=rmse,df.rmse=df.rmse,
              cpd=cpd,df.cpd=df.cpd,
              dss=dss,df.dss=df.dss)
  return(out)
}
