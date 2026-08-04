## --------------------------------------------------------- ##
##
##  FILE 16: seed sensitivity of CCP-splines forecast accuracy
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
source("funs/CCPsplines.R")
source("funs/OutOfSample.R")

## analyse two cases: French males and Swedish females

## all countries and sexes of interest
all.cou <- c("FRATNP","SWE")
all.sex <- c("M","F")

## length of out-of-sample exercise
h.out <- 10

##----- looping over different populations -----

## number of simulations
n.sim <- 250

## width PIs
alpha <- 80

## number of seeds
my.seeds <- 1:50
n.seeds <- length(my.seeds)

## for loop
i <- 1
for (i in 1:length(all.cou)){
  ## select country
  cou <- all.cou[i]
  ## select sex
  sex <- all.sex[i]
  
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
  
  ## fit and forecast with CCP each time with a different seed
  my.s <- 2
  for (my.s in 1:n.seeds){
    ## set seed
    my.seed <- my.seeds[my.s]
    set.seed(my.seed)
    cat("seed",my.seed , "/", n.seeds,"\n")
    ## fitting the CCP-splines model
    CCPsplines <- CCPsplines_fun(ages=x,cohorts=c1,Z=cZ1,E=cE1,sex=sex,
                                 simulate = T,n.sim = n.sim)
    
    ## compute out-of-sample statistics 
    oos_CCP <- OutOfSample_function(ages=x,cohorts = c,h.out = h.out,cMX.obs = cMX,
                                    cMX.sim=exp(CCPsplines$ETA.sim),alpha=alpha)
    
    ## saving results of interest
    df.res.temp <- tibble(cou=cou,sex=sex,model=c("CCP"),
                          seed = my.seed,
                          rmse=round(oos_CCP$rmse,2),
                          cpd=round(oos_CCP$cpd,2),
                          dss=round(oos_CCP$dss,2))
    
    ## saving results
    if (i==1 & my.s==1){
      df.res <- df.res.temp
    }else{
      df.res <- df.res %>% 
        bind_rows(df.res.temp)
    }
    
  }

}

## saving
save(df.res,file = "results/16_out10y_CCP_seed.Rdata")










## END
