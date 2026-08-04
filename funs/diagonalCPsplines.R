


## --------------------------------------------------------- ##
##
## Functions for fitting the diagonal CP-splines
##
## --------------------------------------------------------- ##

# ## cleaning the workspace
# rm(list=ls(all=TRUE))
# 
# ## set up the directory where .R is saved (R-studio command)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("~/WORK/CohortMortalityForecasting")
# 
# ## load packages
# library(tidyverse)
# library(forecast)
# library(svcm)
# library(MortalitySmooth)
# library(MASS)
# 
# ## load useful functions
# source("funs/diagonalLC.R")
# source("funs/cohortLC.R")
# source("funs/cohortLLC.R")
# source("funs/CCPsplines.R")
# source("funs/OutOfSample.R")
# source("funs/LifetableMX.R")
# ## all countries and sexes of interest
# all.cou <- c("AUS","FRATNP","SWE","USA")
# all.sex <- c("F","M")
# 
# ## length of out-of-sample exercise
# h.out <- 10
# 
# ##----- looping over different populations -----
# 
# ## set seed for reproducibility
# set.seed(1)
# 
# ## number of simulations
# n.sim <- 250
# 
# ## width PIs
# alpha <- 80
# 
# ## for loop
# i <- j <- 1
# 
# ## select country
# cou <- all.cou[i]
# 
# ## select sex
# sex <- all.sex[j]
# 
# ## print current exercise
# cat("analysing",cou,sex,"\n")
# 
# ## loading data
# load(file=paste0("data/input/",cou,"_",sex,".Rdata"))
# 
# ## select out-of-sample data
# last.t <- max(c)-h.out
# t.out <- (last.t+1):max(c)
# c1 <- c[1]:(last.t)
# nc1 <- length(c1)
# cZ1 <- cZ[,c%in%c1]
# cE1 <- cE[,c%in%c1]
# t1 <- t[1]:last.t
# nt1 <- length(t1)
# Z1 <- Z[,t%in%t1]
# E1 <- E[,t%in%t1]
# 
# ## matrices containing index of APC
# AGES <- matrix(x,nrow=m,ncol=nc1)
# COHORTS <- matrix(rep(c1,each=m),nrow=m)
# PERIODS <- COHORTS + AGES
# 
# ## set to NAs data after last observed year 
# cZ1[PERIODS>last.t] <- NA
# cE1[PERIODS>last.t] <- NA
# cMX1 <- cZ1/cE1
# cLMX1 <- log(cMX1)

## fitting the diagonal CP-splines

# ## arguments
# ages=x
# years=t1
# cohorts=c1
# Z=Z1
# E=E1
# cZ=cZ1
# cE=cE1
# sex=sex
# n.sim=n.sim
# yrs.exclude=NULL
# decrease=TRUE
# levels=c(95,50)
# print.last=TRUE
# print.all=TRUE
# print.sim=TRUE

diagonal_CPS_fun <- function(ages, years, cohorts, 
                             Z, E, cZ, cE, sex,
                             n.sim=100,
                             yrs.exclude=NULL,
                             decrease=TRUE,
                             levels=c(95,50),
                             print.last=TRUE,
                             print.all=FALSE,
                             print.sim=TRUE){
  ## printing the model being fitted
  cat("fitting the diagonal-CP-splines model", '\n')
  ## dimensions
  m <- length(ages)
  n <- length(years)
  nc <- length(cohorts)
  ## find number of forecast years needed to complete cohorts
  t.fore <- (years[n]+1):(cohorts[nc]+m+10)
  n.fore <- length(t.fore)
  
  E2 <- matrix(999, m, n+n.fore)
  E2[,1:n] <- E
  Z2 <- matrix(99, m, n+n.fore)
  Z2[,1:n] <- Z
  WEI <- matrix(0, m, n+n.fore)
  WEI[,1:n] <- 1
  WEI[E==0] <- 0
  WEI[is.na(Z)] <- 0
  WEI1 <- WEI[,1:n]
  ## eventually add weights=0 for user-defined years (e.g. wars)
  if(!is.null(yrs.exclude)){
    # WEI1 <- WEI
    ## matrix with years
    yrs.mat <- row(Z)+col(Z)+c[1]-2
    ## where to exclude
    whi.exclude <- which(yrs.mat%in%yrs.exclude)
    ## include in the weights
    WEI[whi.exclude] <- 0
  }
  ## to get good starting values for the lambdas-search
  if (print.all) cat("Getting starting values ...", "\n")
  FIT0 <- PSinfant(Y=Z, E=E, WEI=WEI1, lambdas=c(1,100), verbose=FALSE)
  
  ## function to extract the BIC for a given lambda
  BIC1 <- function(par){
    FIT1 <- PSinfant(Y=Z, E=E, lambdas=par, WEI=WEI1, ALPHAS.st=FIT0$ALPHAS)
    FIT1$bic
  }
  ## optimizing lambdas using greedy grid search
  if (print.all) cat("Optimizing smoothing parameters ...", "\n")
  OPT1 <- cleversearch(BIC1, lower=c(-4, 1), upper=c(0, 5),
                       ngrid=5, logscale=TRUE, verbose=FALSE)
  ## optimal smoothing parameters lambdas
  lambdas.hat <- OPT1$par
  if (print.all) cat("Optimal smoothing parameters:", lambdas.hat, "\n")
  ## estimating mortality with optimal lambdas
  if (print.all) cat("Fitting observed data and computing derivatives ...", "\n")
  FIT1 <- PSinfant(Y=Z, E=E, WEI=WEI1, lambdas=lambdas.hat, 
                   ALPHAS.st=FIT0$ALPHAS, verbose=FALSE)
  ## extract deltas
  ETA1a <- FIT1$ETA1a
  ETA1t <- FIT1$ETA1t
  ## compute levels and deltas over ages
  p.a.up <- (100 - (100-levels[1])/2)/100
  p.a.low <- ((100-levels[1])/2)/100
  delta.a.up <- apply(ETA1a, 1, quantile,
                      probs=p.a.up, na.rm=TRUE)
  delta.a.low <- apply(ETA1a, 1, quantile,
                       probs=p.a.low, na.rm=TRUE)
  ## compute levels and deltas over years
  p.t.up <- (100 - (100-levels[2])/2)/100
  p.t.low <- ((100-levels[2])/2)/100
  delta.t.up <- apply(ETA1t, 1, quantile,
                      probs=p.t.up, na.rm=TRUE)
  delta.t.low <- apply(ETA1t, 1, quantile,
                       probs=p.t.low, na.rm=TRUE)
  if(decrease){
    delta.t.up[delta.t.up>0] <- 0
    delta.t.low[delta.t.low>0] <- 0
  }
  deltas <- list(delta.a.up=delta.a.up,delta.a.low=delta.a.low,
                 delta.t.up=delta.t.up,delta.t.low=delta.t.low)
  
  S <- WEI
  S <- 1-S
  if (print.all) cat("Fitting and forecasting data ...", "\n")
  FIT <- CPSfunction(Y=Z2, E=E2, WEI=WEI, lambdas=lambdas.hat,
                     deltas=deltas, S=S, verbose=FALSE)
  
  if (print.last) cat("converged at iter.", FIT$it, ", conv. criteria", FIT$d.dev, '\n')
  ## compute e0 and e-dagger for the point estimates
  if (print.all) cat("Computing cohort estimated log-mortality ...", "\n")
  
  ## extract cohort forecasts
  # library(fields)
  ## observed rates
  MX <- Z/E
  cMX <- cZ/cE
  
  ## fitted rates (only observed period)
  MX.fit <- exp(FIT$ETA[,1:n])
  # library(fields)
  # image.plot(t1,x,t(log(MX)))
  # image.plot(t1,x,t(log(MX.fit)))
  cMX.fit <- matrix(NA,nrow(cMX),ncol(cMX))
  ## for the last observed year, cohort data can't be computed as it needs 
  ## the first forecast. Hence include the fitted forecast as obtained from CP-splines
  MX.fore1 <- exp(FIT$ETA[,n+1])
  MX.fit.ext <- cbind(MX.fit,MX.fore1)
  # image(c(t1,2010),x,t(log(MX.fit.ext)))
  AGES <- matrix(ages,nrow=m,ncol=n+1)
  PERIODS <- matrix(rep(c(years,years[n]+1),each=m),nrow=m)
  COHORTS  <-  PERIODS - AGES
  c.temp <- COHORTS[1,1]:max(cohorts)
  nc.temp <- length(c.temp)
  j <- i <- 1
  for (i in 1:nc.temp){
    my.coh <- which(c==c.temp[i])
    for (j in 1:m){
      whi <- which(COHORTS[j,]==c.temp[i])
      whi1 <- which(COHORTS[j,]==c.temp[i]+1) 
      if (length(whi)!=0 & length(whi1)!=0){
        mx.temp <- MX.fit.ext[j,whi]
        mx.temp1 <- MX.fit.ext[j,whi1]
        ## using Carl's adjustment (AP2 estimator)
        cMX.fit[j,my.coh] <- (mx.temp + mx.temp1) / 2
      }
    }
  }
  
  if (print.all) cat("Simulating coefficients ...", "\n")
  require(MASS)
  ALPHAS.sim <- mvrnorm(n=n.sim, 
                        mu=c(FIT$ALPHAS),
                        Sigma=FIT$Valphas)
  Ba <- FIT$Ba
  Bt <- FIT$Bt
  ## empty matrices/arrays to store bootstrap results
  ETA.sim <- array(NA,dim=c(nrow(Z),ncol(Z2),n.sim))
  for(sim in 1:n.sim){
    alphas.s <- ALPHAS.sim[sim,]
    ALPHAS.s <- matrix(alphas.s, ncol(Ba), ncol(Bt))
    ETA.s <- MortSmooth_BcoefB(Ba, Bt, ALPHAS.s)
    ## saving outcomes
    ETA.sim[,,sim] <- ETA.s
  }
  
  if (print.all) cat("Computing cohort simulated log-mortality, ex and e-dagger ...", "\n")
  ## all time periods, and APC matrices
  t <- years
  t.all <- c(t,t.fore)
  nt.all <- length(t.all)
  AGES <- matrix(ages,nrow=m,ncol=nt.all)
  PERIODS <- matrix(rep(t.all,each=m),nrow=m)
  COHORTS  <-  PERIODS - AGES
  
  ## objects to store results
  ETA <- LE <- ED <- array(NA,dim=c(m,nc,n.sim))
  sim <- 1
  for (sim in 1:n.sim){
    if (print.sim & sim %in% round(seq(0,n.sim,length.out = 11))) cat("simulation",sim,"/",n.sim,"\n")
    MX.temp <- exp(ETA.sim[,,sim])
    cMX.temp <- matrix(NA,m,nc)
    i <- 1
    for (i in 1:nc){
      whi <- which(COHORTS==c[i])
      whi1 <- which(COHORTS==c[i]+1) 
      mx.temp <- MX.temp[whi]
      mx.temp1 <- MX.temp[whi1]
      if (length(mx.temp)<m){
        mx.temp <- c(rep(NA,m-length(mx.temp)),mx.temp)
        mx.temp1 <- c(rep(NA,m-length(mx.temp1)),mx.temp1)
      } 
      ## using Carl's adjustment (AP2 estimator)
      cMX.temp[,i] <- (mx.temp + mx.temp1) / 2
    }
    
    ## create dataframe to store forecast and select only forecast
    ## (i.e. keep observed data otherwise)
    df.temp <- tibble(ages=rep(ages,nc),cohorts=rep(cohorts,each=m),
                      observed=c(cMX.fit),forecast=c(cMX.temp)) %>% 
      mutate(rates=case_when(
        is.na(observed) ~ forecast,
        TRUE            ~ observed
      ))
    
    ETA[,,sim] <- log(df.temp$rates)
    
    ## life-table
    for (i in 1:nc){
      if (!is.na(ETA[1,i,sim])){
        clt <- lifetable.mx(x=x,mx=exp(ETA[,i,sim]),sex=sex)
        LE[,i,sim] <- clt$ex
        ED[,i,sim] <- eDagger(clt)
      }
    }
    
  }
  ## output
  out <- list(LE.sim=LE,ED.sim=ED,ETA.sim=ETA)
}


# bla <- diagonal_CPS_fun(ages=x,years=t1,cohorts=c1,
#                         Z=Z1,E=E1,cZ=cZ1,cE=cE1,
#                         sex=sex)
