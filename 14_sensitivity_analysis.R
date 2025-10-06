## --------------------------------------------------------- ##
##
##  FILE 14: making sensitivity analysis to be plotted in Fig A1
##
##  sessionInfo() details:
##  
##  R version 4.4.0 (2024-04-24 ucrt)
##  Platform: x86_64-w64-mingw32/x64
##  Running under: Windows 11 x64 (build 26100)
##  
##  locale: LC_COLLATE=English_United Kingdom.utf8
##
##  attached base packages:
##  splines stats graphics grDevices utils datasets methods base 
## 
##  other attached packages:
##  lubridate_1.9.3   forcats_1.0.0    
##  stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.5      
##  tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0  
##
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading packages
library(tidyverse)

## load useful functions
source("funs/cohortLC.R")

## ---- two cases: SWE F and FRATNP M-----
my.cases <- c("SWE_F","FRATNP_M")

## loop over cases
case <- "FRATNP_M"
for (case in my.cases){
  
  ## load data 
  load(file=paste0("data/input/",case,".Rdata"))
  
  ## remove all data prior to a given year, and save results
  
  ## matrices containing index of APC
  AGES <- matrix(x,nrow=m,ncol=nc)
  COHORTS <- matrix(rep(c,each=m),nrow=m)
  PERIODS <- COHORTS + AGES
  
  ## years to loop over
  my.years <- seq(1880,1950)
  n <- length(my.years)
  
  ## objects to store results
  ALPHAS <- BETAS <- matrix(NA,m,n)
  KAPPAS <- matrix(NA,nc,n)
  colnames(ALPHAS) <- colnames(BETAS) <- colnames(KAPPAS) <- my.years
  
  ## starting lambda.hat for all populations
  lambda.hat <- 0
  
  ## for loop
  i <- 1
  for (i in 1:n){
    
    ## start year
    cat("fitting with data from", my.years[i],"\n")
    
    ## set to NAs data before specific year 
    first.t <- my.years[i]
    cZ1 <- cZ
    cE1 <- cE
    cZ1[PERIODS<first.t] <- NA
    cE1[PERIODS<first.t] <- NA
    cLMX1 <- log(cZ1/cE1)
    
    ## fit model
    fit_cohort_LC <- cohort_LC_fun(ages=x,cohorts=c,Z=cZ1,E=cE1,sex=sex,simulate=F)
    
    if (!fit_cohort_LC$conv){
      ## if no convergence, loop over higher lambdas until convergence
      j <- 1
      for (j in 1:6){
        ## increase lambda, starting from lambda of previous fitting period
        lambda.hat <- max(10^(j+1),lambda.hat)
        fit_cohort_LC <- cohort_LC_fun(ages=x,cohorts=c,Z=cZ1,E=cE1,sex=sex,simulate=F,
                                       lambdaA = lambda.hat, lambdaB = lambda.hat)
        ## break if convergence and update lambda hat
        if (fit_cohort_LC$conv){
          lambda.hat <- fit_cohort_LC$lambda
          break
        } 
      }
    } 
    
    ## save 
    if (fit_cohort_LC$conv){
      ALPHAS[,i] <- fit_cohort_LC$Alpha
      BETAS[,i] <- fit_cohort_LC$Beta
      KAPPAS[,i] <- fit_cohort_LC$Kappa
    }
  }
  
  ## saving results
  df.alphas.temp <- tibble(par="Alpha",cou=case,scale=x) %>% 
    bind_cols(ALPHAS) 
  df.betas.temp <- tibble(par="Beta",cou=case,scale=x) %>% 
    bind_cols(BETAS) %>% mutate()
  df.kappas.temp <- tibble(par="Kappa",cou=case,scale=c) %>% 
    bind_cols(KAPPAS)
  df.res.temp <- df.alphas.temp %>% 
    bind_rows(df.betas.temp) %>% 
    bind_rows(df.kappas.temp)
  
  if (which(case==my.cases)==1){
    df.res <- df.res.temp
  }else{
    df.res <- df.res %>% 
      bind_rows(df.res.temp)
  }
}

## saving
save(df.res,file = "results/14_sensitivity.Rdata")

## END
