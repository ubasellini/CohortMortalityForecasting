## --------------------------------------------------------- ##
##
##  FILE 23: analysis and plot for Figure 3
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
##  forecast_8.22.0   patchwork_1.2.0   viridis_0.6.5     viridisLite_0.4.2
##  lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4      
##  purrr_1.0.2       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1     
##  ggplot2_3.5.1     tidyverse_2.0.0 
##
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading packages
library(tidyverse)
library(viridis)
library(patchwork)
library(forecast)

## loading cohort LC functions
source("funs/cohortLC.R")

## loading overall data
load(file=paste0("data/input/all_countries",".Rdata"))

## simulations for CIs
n.sim <- 250

## width CIs
alpha.lev <- 95/100

## ---- Swedish females  -----
## load data 
load(file=paste0("data/input/SWE_F",".Rdata"))

## fitting the cohort LC model
fit_cohort_LC <- cohort_LC_fun(ages=x,cohorts=c,Z=cZ,E=cE,
                               sex=sex,simulate=T,n.sim=n.sim)

## parameters plots
par(mfrow=c(1,3))
plot(x,fit_cohort_LC$Alpha)
plot(x,fit_cohort_LC$Beta);abline(h=0)
plot(c,fit_cohort_LC$Kappa)
par(mfrow=c(1,1))

## ggplot for parameters(in inset)
df.pars.swe <- tibble(Scale=c(x,x,c),
                  Value=c(fit_cohort_LC$Alpha,fit_cohort_LC$Beta,fit_cohort_LC$Kappa),
                  Index=c(rep("alpha[x]",m),rep("beta[x]",m),rep("kappa[c]",nc)),
                  Country="SWE")

fig_pars_swe <- df.pars.swe %>% 
  ggplot(aes(x=Scale,y=Value))+
  geom_line(linewidth=1)+
  facet_wrap(.~Index,scales = "free",labeller = label_parsed)+
  theme_bw(base_size = 16)+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=8))+
  labs(x=NULL,y=NULL)
fig_pars_swe

## median, lower and upper ETA
ETA.sim.upp <- ETA.sim.low <- ETA.sim.med <- matrix(NA,m,nc)
i <- 1
for (i in 1:m){
  ETA.sim.med[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, median, na.rm = T)
  ETA.sim.upp[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, quantile, prob=1- (1-alpha.lev)/2, na.rm = T)
  ETA.sim.low[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, quantile, prob=(1-alpha.lev)/2, na.rm = T)
}

## storing predictions
df.fit <- tibble(Country="SWE",Sex="F",
                 Year=rep(c,each=m),Age=rep(x,nc),
                 LC.SWE=c(exp(fit_cohort_LC$ETA)),
                 LC.SWE_med=c(exp(ETA.sim.med)),
                 LC.SWE_upp=c(exp(ETA.sim.upp)),
                 LC.SWE_low=c(exp(ETA.sim.low)))
Cohort.all <- Cohort.all %>% 
  left_join(df.fit)

## plot
Cohort.all %>% 
  filter(Sex=="F") %>% 
  filter(Country=="SWE") %>% 
  filter(Age%in%c(40,60,80)) %>% 
  mutate(Age=factor(Age,levels=c("40","60","80"),
                    labels=c("Age 40","Age 60", "Age 80"))) %>%
  ggplot(aes(x=Year, group=Age))+
  geom_ribbon(aes(ymin=LC.SWE_low, ymax=LC.SWE_upp),
              alpha=0.45, fill="lightpink") +
  geom_line(aes(y=LC.SWE),color=2,linewidth=0.95) +
  geom_line(aes(y=LC.SWE_med),color=3,linewidth=0.95,linetype="dashed") +
  geom_point(aes(y=Rate,shape = Age,color=Age),size=2.5)+
  scale_y_log10()

## ---- USA females  -----
## load data 
load(file=paste0("data/input/USA_F",".Rdata"))

## try fitting the cohort LC model
fit_cohort_LC <- cohort_LC_fun(ages=x,cohorts=c,Z=cZ,E=cE,sex=sex)

## since the model does not converge, we include smoothing in the 
## estimation procedure, and select smoothing parameters by minimising BIC

## series of lambdas
lambdasA <- lambdasB <- 10^seq(1,7)
nla <- nlb <- length(lambdasA)
## array to store BICs
BICs <- array(NA, dim=c(nla, nlb),
              dimnames=list(lambdasA, lambdasB))
## looping over lambas
for (la in 1:nla){
  for (lb in 1:nlb){
    ## print smoothin parameters
    cat("la =",lambdasA[la],"lb =",lambdasA[lb],"\n")
    ## fitting the smooth cohort LC model
    fit_cohort_LC <- fit_cohort_LC_fun(ages=x,cohorts=c,Z=cZ,E=cE,sex=sex,
                                       lambdaA = lambdasA[la],
                                       lambdaB = lambdasA[lb],
                                       max.iter = 750)
    ## save BIC if converged
    if (fit_cohort_LC$conv) BICs[la,lb] <- fit_cohort_LC$BIC
    
  }
}
## choose optimal BIC 
minBIC <- which(BICs==min(BICs, na.rm=TRUE), arr.ind=TRUE)
## optimal lambdas 
lambdaA.hat <- lambdasA[minBIC[1]]
lambdaB.hat <- lambdasB[minBIC[2]]

## fitting the cohort LC model (computing CIs)
fit_cohort_LC <- cohort_LC_fun(ages=x,cohorts=c,Z=cZ,E=cE,sex=sex,
                               lambdaA = lambdaA.hat,lambdaB = lambdaB.hat,
                               simulate=T,n.sim=n.sim)

## parameters plots
par(mfrow=c(1,3))
plot(x,fit_cohort_LC$Alpha)
plot(x,fit_cohort_LC$Beta);abline(h=0)
plot(c,fit_cohort_LC$Kappa);abline(v=1930)
par(mfrow=c(1,1))

## ggplot for parameters(in inset)
df.pars.usa <- tibble(Scale=c(x,x,c),
                  Value=c(fit_cohort_LC$Alpha,fit_cohort_LC$Beta,fit_cohort_LC$Kappa),
                  Index=c(rep("alpha[x]",m),rep("beta[x]",m),rep("kappa[c]",nc)),
                  Country="USA")

fig_pars_usa <- df.pars.usa %>% 
  ggplot(aes(x=Scale,y=Value))+
  geom_line(linewidth=1)+
  facet_wrap(.~Index,scales = "free",labeller = label_parsed)+
  theme_bw(base_size = 20)+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=8))+
  labs(x=NULL,y=NULL)
fig_pars_usa

## single plot for LC parameters
df.pars <- df.pars.swe %>% 
  bind_rows(df.pars.usa)


fig_pars <- df.pars %>% 
  ggplot(aes(x=Scale,y=Value,color=Country))+
  geom_line(linewidth=1.25)+
  facet_wrap(.~Index,scales = "free",labeller = label_parsed)+
  scale_color_manual("Country", name="",values=c("cyan4","darkgoldenrod3"))+
  geom_text(aes(x = 1945,y=100),data = data.frame(Index = "kappa[c]"),
            label = c("Sweden"),color="cyan4",size=7) +
  geom_text(aes(x = 1875,y=-20),data = data.frame(Index = "kappa[c]"),
            label = c("USA"),color="darkgoldenrod3",size=7) +
  theme_bw(base_size = 22)+
  theme(legend.position="none")+
  # theme(axis.text.y=element_text(size=10),
  #       axis.text.x=element_text(size=8),
  #       legend.position="none")+
  labs(x=NULL,y="Cohort LC parameters")
fig_pars

## median, lower and upper ETA
ETA.sim.upp <- ETA.sim.low <- ETA.sim.med <- matrix(NA,m,nc)
i <- 1
for (i in 1:m){
  ETA.sim.med[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, median, na.rm = T)
  ETA.sim.upp[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, quantile, prob=1- (1-alpha.lev)/2, na.rm = T)
  ETA.sim.low[i,] <- apply(fit_cohort_LC$ETA.sim[i,,], 1, quantile, prob=(1-alpha.lev)/2, na.rm = T)
}

## storing predictions
df.fit <- tibble(Country="USA",Sex="F",
                 Year=rep(c,each=m),Age=rep(x,nc),
                 LC.USA=c(exp(fit_cohort_LC$ETA)),
                 LC.USA_med=c(exp(ETA.sim.med)),
                 LC.USA_upp=c(exp(ETA.sim.upp)),
                 LC.USA_low=c(exp(ETA.sim.low)))
Cohort.all <- Cohort.all %>% 
  left_join(df.fit)

## plot
Cohort.all %>% 
  filter(Sex=="F") %>% 
  filter(Country=="USA") %>% 
  filter(Age%in%c(10,60,80)) %>% 
  mutate(Age=factor(Age,levels=c("10","60","80"),
                    labels=c("Age 10","Age 60", "Age 80"))) %>%
  ggplot(aes(x=Year, group=Age))+
  geom_ribbon(aes(ymin=LC.USA_low, ymax=LC.USA_upp),
              alpha=0.45, fill="lightpink") +
  # geom_line(aes(y=LC.USA),color=2,linewidth=0.95) +
  geom_point(aes(y=Rate,shape = Age,color=Age),size=2.5)+
  scale_y_log10()+
  theme_bw()


##---- final plot ------

## text for age groups
ann_text <- data.frame(Country = factor("Sweden",levels = c("Sweden","USA")),
                       Sex="F",
                       Age=factor("10",levels = c("10","60","80")))

## plotting
fig_rates <- Cohort.all %>% 
  mutate(LC_med=case_when(
    !is.na(LC.SWE_med)~LC.SWE_med,
    !is.na(LC.USA_med)~LC.USA_med,
    TRUE~NA),
    LC_upp=case_when(
      !is.na(LC.SWE_upp)~LC.SWE_upp,
      !is.na(LC.USA_upp)~LC.USA_upp,
      TRUE~NA),
    LC_low=case_when(
      !is.na(LC.SWE_low)~LC.SWE_low,
      !is.na(LC.USA_low)~LC.USA_low,
      TRUE~NA)) %>% 
  filter(Sex=="F") %>% 
  filter(Country=="USA"|Country=="SWE") %>% 
  filter(Age%in%c(10,60,80)) %>% 
  mutate(Age=factor(Age,levels=c("10","60","80"),
                    labels=c("Age 10","Age 60", "Age 80")),
         Country=case_when(
           Country=="SWE"~"Sweden",
           TRUE ~ Country),
         Country=factor(Country)) %>%
  ggplot(aes(x=Year, group=Age))+
  geom_point(aes(y=Rate,shape = Age,color=Age),size=2.5)+
  scale_shape_manual("Age", name="",values=c(16,18,17))+
  scale_color_manual("Age", name="",values=c(1,"grey30","grey60"))+
  geom_ribbon(aes(ymin=LC_low, ymax=LC_upp),
              alpha=0.45, fill="lightpink") +
  geom_line(aes(y=LC_med),color=2,linewidth=0.9) +
  facet_grid(.~Country)+
  scale_y_log10()+
  theme_bw(base_size = 22)+
  geom_text(aes(x = 1925,y=0.0885),data = ann_text,
            label = c("age 80"),color="grey60",size=7) +
  geom_text(aes(x = 1925,y=0.013),data = ann_text,
            label = c("age 60"),color="grey30",size=7) +
  geom_text(aes(x = 1925,y=0.00017),data = ann_text,
            label = c("age 10"),color="grey0",size=7) +
  geom_text(aes(x = 1993,y=0.0205),data = ann_text,
            label = c("cohort LC"),color=2,size=7)+
  theme(legend.position="none") +
  labs(x="Cohort",y="Mortality rate (log scale)")
fig_rates

## include inset plot (left, bottom, right, top)
# fig_rates + inset_element(fig_pars, 0.525, 0.015, 0.975, 0.35)
# fig_rates / (fig_pars_swe+fig_pars_usa)
fig_rates/fig_pars
## saving Figure
ggsave(file="figs/F3.pdf",width = 12,height=10)

## END