## --------------------------------------------------------- ##
##
##  FILE 12: 20-year out-of-sample validation 
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
##  MASS_7.3-60.2         MortalitySmooth_2.3.4 lattice_0.22-6       
##  svcm_0.1.2            Matrix_1.7-0          forecast_8.22.0      
##  lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1        
##  dplyr_1.1.4           purrr_1.0.2           readr_2.1.5          
##  tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1        
##  tidyverse_2.0.0    
##
## --------------------------------------------------------- ##


##----- setting up the basic steps -----

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load packages
library(tidyverse)
library(forecast)
library(svcm)
library(MortalitySmooth)
library(MASS)

## load useful functions
source("funs/diagonalLC.R")
source("funs/cohortLC.R")
source("funs/cohortLLC.R")
source("funs/CCPsplines.R")
source("funs/OutOfSample.R")

## all countries and sexes of interest
all.cou <- c("AUS","FRATNP","SWE","USA")
all.sex <- c("F","M")

## length of out-of-sample exercise
h.out <- 20

##----- looping over different populations -----

## set seed for reproducibility
set.seed(1)

## number of simulations
n.sim <- 250

## width PIs
alpha <- 80

## for loop
i <- j <- 1
for (i in 1:length(all.cou)){
  ## select country
  cou <- all.cou[i]
  for (j in 1:length(all.sex)){
    ## select sex
    sex <- all.sex[j]
  
    ## print current exercise
    cat("analysing",cou,sex,"\n")
    
    ## loading data
    load(file=paste0("data/input/",cou,"_",sex,".Rdata"))
    
    ## select out-of-sample data
    last.t <- max(c)-h.out
    t.out <- (last.t+1):max(c)
    c1 <- c[1]:(last.t)
    nc1 <- length(c1)
    cZ1 <- cZ[,c%in%c1]
    cE1 <- cE[,c%in%c1]
    t1 <- t[1]:last.t
    nt1 <- length(t1)
    Z1 <- Z[,t%in%t1]
    E1 <- E[,t%in%t1]
    
    ## matrices containing index of APC
    AGES <- matrix(x,nrow=m,ncol=nc1)
    COHORTS <- matrix(rep(c1,each=m),nrow=m)
    PERIODS <- COHORTS + AGES
    
    ## set to NAs data after last observed year 
    cZ1[PERIODS>last.t] <- NA
    cE1[PERIODS>last.t] <- NA
    cMX1 <- cZ1/cE1
    cLMX1 <- log(cMX1)
    
    ## ---- fitting the four cohort forecasting models -----

    ## fitting the diagonal LC model
    diagonal_LC <- diagonal_LC_fun(ages=x,years=t1,cohorts=c1,
                                   Z=Z1,E=E1,cZ=cZ1,cE=cE1,sex=sex,
                                   n.sim=n.sim)
    ## fitting the cohort LC model (only for FRATNP and SWE)
    if (cou == "FRATNP"| cou == "SWE"){
      cohort_LC <- cohort_LC_fun(ages=x,cohorts=c1,Z=cZ1,E=cE1,sex=sex,
                                 simulate = T,n.sim = n.sim)
    }
    ## fitting the cohort LLC model
    cohort_LLC <- cohort_LLC_fun(ages=x,cohorts=c1,Z=cZ1,E=cE1,sex=sex,
                                 simulate = T,n.sim = n.sim)
    ## fitting the CCP-splines model
    CCPsplines <- CCPsplines_fun(ages=x,cohorts=c1,Z=cZ1,E=cE1,sex=sex,
                                 simulate = T,n.sim = n.sim)
    
    
    ## compute out-of-sample statistics for different models
    oos_dLC <- OutOfSample_function(ages=x,cohorts = c,h.out = h.out,cMX.obs = cMX,
                                    cMX.sim=exp(diagonal_LC$ETA.sim),alpha=alpha)
    if (cou == "FRATNP"| cou == "SWE"){
      oos_LC <- OutOfSample_function(ages=x,cohorts = c,h.out = h.out,cMX.obs = cMX,
                                     cMX.sim=exp(cohort_LC$ETA.sim),alpha=alpha)
    }else{
      oos_LC <- list(rmse.log=NA,cpd.log=NA,dss.log=NA)
    }
    oos_LLC <- OutOfSample_function(ages=x,cohorts = c,h.out = h.out,cMX.obs = cMX,
                                    cMX.sim=exp(cohort_LLC$ETA.sim),alpha=alpha)
    oos_CCP <- OutOfSample_function(ages=x,cohorts = c,h.out = h.out,cMX.obs = cMX,
                                    cMX.sim=exp(CCPsplines$ETA.sim),alpha=alpha)
    
    ## saving results of interest
    
    ## RMSE
    rmse.temp <- round(c(oos_dLC$rmse,oos_LC$rmse,oos_LLC$rmse,oos_CCP$rmse),2)
    
    ## CPD
    cpd.temp <- round(c(oos_dLC$cpd,oos_LC$cpd,oos_LLC$cpd,oos_CCP$cpd),2)
    
    ## DSS
    dss.temp <- round(c(oos_dLC$dss,oos_LC$dss,oos_LLC$dss,oos_CCP$dss),2)
    
    ## all results
    res.temp <- c(rmse.temp,cpd.temp,dss.temp)
    df.res.temp <- tibble(cou=cou,sex=sex,model=c("dLC","LC","LLC","CCP"),
                          rmse=rmse.temp,cpd=cpd.temp,dss=dss.temp)
    
    ## saving results
    if (i==1 & j==1){
      res <- res.temp
      df.res <- df.res.temp
    }else{
      res <- rbind(res,res.temp) 
      df.res <- df.res %>% 
        bind_rows(df.res.temp)
    }
    
    
  }
  
  
}

## printing values for manuscript table
for (i in 1:nrow(res)){
  print(paste(res[i,], collapse = " & "))
}

## checking
df.res

## saving
save(res,df.res,file = "results/12_out20y.Rdata")

## finding the minimum
df.res %>% 
  pivot_longer(-c(cou,sex,model)) %>% 
  group_by(cou,sex,name) %>% 
  summarise(min   = min(value,na.rm = T),
            model = model[which.min(value)]) %>% print(n=100)

## END
