## --------------------------------------------------------- ##
##
##  FILE 01: download period and cohort mortality data from HMD 
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
##  stats graphics grDevices utils datasets methods base    
## 
##  other attached packages:
##  lubridate_1.9.3  forcats_1.0.0    stringr_1.5.1    dplyr_1.1.4     
##  purrr_1.0.2      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1    
##  ggplot2_3.5.1    tidyverse_2.0.0  HMDHFDplus_2.0.8
##
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading packages
library(HMDHFDplus)
library(tidyverse)

## set your HMD username and password
username <- "XXXXXXXXXXX"
password <- "XXXXXXXXXXX"

## all countries of interest and two sexes
countries <- c("SWE","FRATNP","USA","AUS")
sexes <- c("M","F")
n.cou <- length(countries)
n.sex <- length(sexes)

## common cohorts, time period and age group for all countries analysed
c <-  1850:2019
t <- 1950:2019
x <- 0:100
m <- length(x)
nt <- length(t)
nc <- length(c)

## looping over countries
i <- j <- 1
for (i in 1:n.cou){
  cou <- countries[i]
  ## looping over sexes
  for (j in 1:n.sex){
    sex <- sexes[j]
    
    ## which country looking at
    cat(paste("downloading data for",cou,"-",sex),"\n")
    
    ## ------ PERIOD DATA -------------
    ## download deaths and exposures by period
    Z.temp <- readHMDweb(CNTRY = cou, item = "Deaths_1x1", 
                         fixup = TRUE,
                         username = username, password = password)
    E.temp <- readHMDweb(CNTRY = cou, item = "Exposures_1x1", 
                         fixup = TRUE,
                         username = username, password = password)
    
    ## find first observed year in period data
    first.t <- min(Z.temp$Year)
    
    ## select correct column for females/males
    if (sex=="F"){
      Z.temp <- Z.temp %>% 
        select(Year,Age,Deaths=Female)
      E.temp <- E.temp %>% 
        select(Year,Age,Exposures=Female)
    }else{
      Z.temp <- Z.temp %>% 
        select(Year,Age,Deaths=Male)
      E.temp <- E.temp %>% 
        select(Year,Age,Exposures=Male)
    }
    
    ## deaths
    Z.full <- expand.grid(Age=x,Year=t)
    Z.temp2 <- Z.temp %>% 
      filter(Year %in% t) %>% 
      mutate(AgeGroup=case_when(
        Age>x[m] ~ x[m],
        TRUE    ~ Age
      )) %>% 
      group_by(Year,AgeGroup) %>% 
      summarise(Deaths=sum(Deaths,na.rm = T)) %>% rename(Age=AgeGroup)
    Z.full <- Z.full %>% 
      left_join(Z.temp2,by = join_by(Age, Year)) %>% mutate(Country=cou,Sex=sex)
    Z <- Z.full %>% select(Deaths) %>% 
      pull() %>% matrix(.,nrow = m,ncol=nt)
    ## exposures
    E.full <- expand.grid(Age=x,Year=t)
    E.temp2 <- E.temp %>% 
      filter(Year %in% t) %>% 
      mutate(AgeGroup=case_when(
        Age>x[m] ~ x[m],
        TRUE    ~ Age
      )) %>% 
      group_by(Year,AgeGroup) %>% 
      summarise(Exposures=sum(Exposures,na.rm = T)) %>% rename(Age=AgeGroup)
    E.full <- E.full %>% 
      left_join(E.temp2,by = join_by(Age, Year)) %>% mutate(Country=cou,Sex=sex)
    E <- E.full %>% select(Exposures) %>% 
      pull() %>% matrix(.,nrow = m,ncol=nt)
    ## rates and log-rates
    MX <- Z/E
    LMX <- log(MX)
    
    ## single period dataset
    Period.all.temp <- Z.full %>% 
      left_join(E.full)
    if (j==1 & i==1){
      Period.all <- Period.all.temp
    }else{
      Period.all <- Period.all %>% 
        bind_rows(Period.all.temp)
    }

    ## ------ COHORT DATA -------------
    ## download rates and exposures by cohort
    cMX.temp <- readHMDweb(CNTRY = cou, item = "cMx_1x1", 
                           fixup = TRUE,
                           username = username, password = password)
    
    cE.temp <- readHMDweb(CNTRY = cou, item = "cExposures_1x1", 
                           fixup = TRUE,
                           username = username, password = password)

    ## select correct column for females/males
    if (sex=="F"){
      cMX.temp <- cMX.temp %>% 
        select(Year,Age,Rate=Female)
      cE.temp <- cE.temp %>% 
        select(Year,Age,Exposures=Female)
    }else{
      cMX.temp <- cMX.temp %>% 
        select(Year,Age,Rate=Male)
      cE.temp <- cE.temp %>% 
        select(Year,Age,Exposures=Male)
    }
    
    ## single tibble
    df.cohort <- cMX.temp %>% 
      left_join(cE.temp,by = join_by(Year, Age)) %>% 
      mutate(Deaths=Rate*Exposures,
             AgeGroup=case_when(
               Age>x[m] ~ x[m],
               TRUE    ~ Age
             )) %>% 
      group_by(Year,AgeGroup) %>% 
      summarise(Exposures=sum(Exposures,na.rm = T),
                Deaths=sum(Deaths,na.rm = T),
                Rate=Deaths/Exposures) %>% rename(Age=AgeGroup) %>% 
      filter(Year%in%c) %>% mutate(Country=cou,Sex=sex)
    
    ## compute cohort deaths and exposures for most recent data
    ## using Lexis triangles and HMD protocol
    
    ## download deaths and exposures by Lexis triangles from HMD
    Z.lexis <- readHMDweb(CNTRY = cou, item = "Deaths_lexis", 
                          fixup = TRUE,
                          username = username, password = password)
    POP <- readHMDweb(CNTRY = cou, item = "Population", 
                      fixup = TRUE,
                      username = username, password = password)
    
    ## select correct column for females/males
    if (sex=="F"){
      Z.lexis <- Z.lexis %>% 
        select(Year,Age,Cohort,Deaths=Female)
      POP <- POP %>% 
        select(Year,Age,Pop=Female1)
    }else{
      Z.lexis <- Z.lexis %>% 
        select(Year,Age,Cohort,Deaths=Male)
      POP <- POP %>% 
        select(Year,Age,Pop=Male1)
    }
    
    ## subset of cohorts
    c.sub <- c[!c%in%df.cohort$Year]
    nc.sub <- length(c.sub)
    
    ## most recent cohorts
    c.sub.r <- c.sub[c.sub>1950]
    nc.sub.r <- length(c.sub)
    
    ## cohort deaths from Lexis
    df.deaths.coh <- Z.lexis %>% 
      filter(Cohort%in%c.sub) %>% 
      group_by(Cohort,Age) %>% 
      summarise(Deaths=sum(Deaths,na.rm = T))
    
    ## cohort exposures from Lexis (assuming uniform distribution of births)
    cE.sub <- matrix(NA,m,nc.sub)
    k <- l <- 1
    for (k in 1:nc.sub){
      my.cohort <- c.sub[k]
      l <- m
      for (l in 1:m){
        my.age <- x[l]
        my.year <- my.cohort+my.age
        if (my.year<max(Z.lexis$Year)){
          t1 <- POP %>% 
            filter(Age==my.age,Year==(my.year+1)) %>% 
            select(Pop) %>% pull()
          t2 <- Z.lexis %>% 
            filter(Age==my.age,Year==(my.year),Cohort==my.cohort) %>%
            select(Deaths) %>% pull()
          if (length(t1)>0 & length(t2)>0){
            t3 <- Z.lexis %>% 
              filter(Age==my.age,Year==(my.year+1),Cohort==my.cohort) %>%
              select(Deaths) %>% pull()
            cE.sub[l,k] <- t1 + 1/3*(t2-t3)
          }else{
            cE.sub[l,k] <- NA
          }
        }else{
          break
        }
      }
    }
    
    df.expo.coh <- expand.grid(Age=x,Cohort=c.sub) %>% 
      mutate(Exposures=c(cE.sub))
    
    df.cohort.sub <- df.expo.coh %>% 
      left_join(df.deaths.coh) %>% 
      rename(Year=Cohort) %>% 
      mutate(Rate=Deaths/Exposures,Country=cou,Sex=sex)
    
    df.cohort.fin <- df.cohort %>% 
      bind_rows(df.cohort.sub) %>% 
      mutate(Period=Age+Year,
             Exposures=case_when(
               Period > 2019 ~ NA,
               Period < first.t ~ NA,
               TRUE ~ Exposures),
             Deaths=case_when(
               Period > 2019 ~ NA,
               Period < first.t ~ NA,
               TRUE ~ Deaths),
             Rate=case_when(
               Period > 2019 ~ NA,
               Period < first.t ~ NA,
               TRUE ~ Rate)) %>% 
      arrange(Year,Age)
    

    ## matrices
    cMX <- df.cohort.fin %>% ungroup() %>% select(Rate) %>% 
      pull() %>% matrix(.,nrow = m,ncol=nc)
    cE <- df.cohort.fin %>% ungroup() %>% select(Exposures) %>% 
      pull() %>% matrix(.,nrow = m,ncol=nc)
    cZ <- df.cohort.fin %>% ungroup() %>% select(Deaths) %>% 
      pull() %>% matrix(.,nrow = m,ncol=nc)
    cLMX <- log(cMX)
    
    ## single cohort dataset
    if (j==1 & i==1){
      Cohort.all <- df.cohort.fin
    }else{
      Cohort.all <- Cohort.all %>% 
        bind_rows(df.cohort.fin)
    }
    
    ## saving individual-country data
    save(cE,cLMX,cMX,cZ,E,LMX,MX,Z,
         c,cou,m,nc,nt,sex,t,x,
         file=paste0("data/input/",cou,"_",sex,".Rdata"))

  }
}

## saving overall data
save(Period.all,Cohort.all,file=paste0("data/input/all_countries",".Rdata"))

## END