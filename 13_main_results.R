## --------------------------------------------------------- ##
##
##  FILE 13: generating the main results of the paper 
##           (Section 3.2)
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
source("funs/CCPsplines.R")
source("funs/CIs.R")

## all countries and sexes of interest
all.cou <- c("AUS","FRATNP","SWE","USA")
all.sex <- c("F","M")


##----- looping over different populations -----

## set seed for reproducibility
set.seed(1)

## number of simulations
n.sim <- 250

## width PIs
alpha <- 80

## for loop
i <- j <- 2
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
    
    ## ---- fitting two cohort forecasting models -----

    ## fitting the diagonal LC model and extract median and CIs
    diagonal_LC <- diagonal_LC_fun(ages=x,years=t,cohorts=c,
                                   Z=Z,E=E,cZ=cZ,cE=cE,sex=sex,
                                   n.sim=n.sim)
    dLC_out <- sim_to_ci(diagonal_LC$ETA.sim,diagonal_LC$LE.sim,alpha = alpha)
    ## fitting the CCP-splines model and extract median and CIs
    CCPsplines <- CCPsplines_fun(ages=x,cohorts=c,Z=cZ,E=cE,sex=sex,
                                 simulate = T,n.sim = n.sim)
    CCP_out <- sim_to_ci(CCPsplines$ETA.sim,CCPsplines$LE.sim,alpha = alpha)
    
    ## saving temporary results for eta and e0
    df.lmx.temp <- tibble(cou=cou,sex=sex,age=rep(x,nc),cohort=rep(c,each=m),
                              obs=c(cLMX),
                              dLC_med=c(dLC_out$ETA.sim.med),
                              dLC_upp=c(dLC_out$ETA.sim.upp),
                              dLC_low=c(dLC_out$ETA.sim.low),
                              CCP_med=c(CCP_out$ETA.sim.med),
                              CCP_upp=c(CCP_out$ETA.sim.upp),
                              CCP_low=c(CCP_out$ETA.sim.low))
    ## plot rates over age 
    # whi.coh <- which(c%in%c(1950,2019))
    # matplot(x,cLMX[,whi.coh],t="p",pch=1,ylim=range(CCP_out$ETA.sim.med[,whi.coh]))
    # matlines(x,CCP_out$ETA.sim.med[,whi.coh],t="l",col=2:3)
    ## plot rates over time 
    # whi.age <- which(x%in%c(30,40))
    # matplot(c,t(cLMX[whi.age,]),t="p",pch=1,ylim=range(CCP_out$ETA.sim.med[whi.age,]))
    # matlines(c,t(CCP_out$ETA.sim.med[whi.age,]),t="l",lty=1)
    # legend("bottomleft",legend=x[whi.age],pch=1,col=1:2)
    # 
    # whi.age <- which(x==40)
    # pdf("figs/USAf_age40.pdf",width=10,height=6)
    # plot(c,t(cLMX[whi.age,]),t="p",pch=1,ylim=range(CCP_out$ETA.sim.med[whi.age,]))
    # lines(c,t(CCP_out$ETA.sim.med[whi.age,]),t="l",col=2:3)
    # dev.off()
    
    
    if (cou=="SWE"|cou=="FRATNP"){
      c.full <- c[1]:1919
      nc.full <- length(c.full)
      e0.obs <- apply(cMX[,1:nc.full],2,e0.mx,x=x,sex=sex)
      e0.obs <- c(e0.obs,rep(NA,nc-nc.full))
    }else{
      e0.obs <- rep(NA,nc)
    }
    df.e0.temp <- tibble(cou=cou,sex=sex,cohort=c,
                         obs=e0.obs,
                         dLC_med=dLC_out$e0.sim.med,
                         dLC_upp=dLC_out$e0.sim.upp,
                         dLC_low=dLC_out$e0.sim.low,
                         CCP_med=CCP_out$e0.sim.med,
                         CCP_upp=CCP_out$e0.sim.upp,
                         CCP_low=CCP_out$e0.sim.low)
    
    
    ## saving results
    if (i==1 & j==1){
      df.lmx <- df.lmx.temp
      df.e0 <- df.e0.temp
    }else{
      df.lmx <- df.lmx %>% 
        bind_rows(df.lmx.temp) 
      df.e0 <- df.e0 %>% 
        bind_rows(df.e0.temp)
    }
  }
}


## saving
save(df.lmx,df.e0,file = "results/13_main_results.Rdata")

## END
