## --------------------------------------------------------- ##
##
##  FILE 25: plotting Figure 5
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
##  patchwork_1.2.0   viridis_0.6.5     viridisLite_0.4.2 lubridate_1.9.3   forcats_1.0.0    
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
library(viridis)
library(ggtext)
library(ggrepel)

## loading previously saved data 
load(file = "results/11_out10y.Rdata")

baseline.long.diagLC <- df.res %>%
  filter(model=="dLC") %>% 
  filter((cou=="FRATNP" & sex=="M")|(cou=="SWE" & sex=="F")) %>% 
  mutate(
    population = case_when(
      cou=="FRATNP" & sex=="M" ~ "France (M)",
      cou=="SWE" & sex=="F" ~ "Sweden (F)"
    )) %>% 
  pivot_longer(
    cols=c(rmse,cpd,dss),
    names_to="measure",
    values_to="value"
  ) %>% 
  mutate(measure=factor(measure,levels=c("rmse","cpd","dss"),
                        labels=c("RMSE","CPD 80%","DSS")))

## loading previously saved data 
load(file = "results/15_out10y_CCP_robustness.Rdata")

## preparing plot
df.plot <- df.res %>%
  mutate(
    population = case_when(
      cou=="FRATNP" & sex=="M" ~ "France (M)",
      cou=="SWE" & sex=="F" ~ "Sweden (F)"
    )
  ) %>%
  pivot_longer(
    cols=c(rmse,cpd,dss),
    names_to="measure",
    values_to="value"
  ) %>% 
  mutate(measure=factor(measure,levels=c("rmse","cpd","dss"),
                        labels=c("RMSE","CPD 80%","DSS")))

baseline.long <- df.res %>%
  mutate(
    population = case_when(
      cou=="FRATNP" & sex=="M" ~ "France (M)",
      cou=="SWE" & sex=="F" ~ "Sweden (F)"
    )) %>% 
  filter(
    lambda.com=="hat",
    delta.t==50,
    upper==TRUE,
  ) %>%
  pivot_longer(
    cols=c(rmse,cpd,dss),
    names_to="measure",
    values_to="value"
  ) %>% 
  mutate(measure=factor(measure,levels=c("rmse","cpd","dss"),
                        labels=c("RMSE","CPD 80%","DSS")))

lambda.labels <- c(
  hat = "(\u03BB<sub>x</sub>, \u03BB<sub>c</sub>)",
  HH  = "(10\u03BB<sub>x</sub>, 10\u03BB<sub>c</sub>)",
  HL  = "(10\u03BB<sub>x</sub>, 0.1\u03BB<sub>c</sub>)",
  LH  = "(0.1\u03BB<sub>x</sub>, 10\u03BB<sub>c</sub>)",
  LL  = "(0.1\u03BB<sub>x</sub>, 0.1\u03BB<sub>c</sub>)"
)

## plotting
df.plot %>%
  mutate(
    lambda.plot = lambda.labels[lambda.com],
    delta.t=factor(delta.t)) %>% 
  ggplot(aes(x=population,y=value)) +
  geom_boxplot(
    width=0.6,
    alpha=0,
    color="grey60",
    outlier.shape = NA
  ) +
  geom_point(
    data = baseline.long.diagLC,
    aes(x = population, y = value),
    size = 4,
    shape = 23,
    fill = "white",
    color = "grey20",
    stroke = 1.5
  )+
  # individual sensitivity runs
  geom_jitter(
    aes(color=lambda.com,shape=delta.t),
    width=0.12,size=2.5,alpha=0.75) +
  facet_wrap(~measure,nrow=1,scales="free_y") +
  geom_text(
    data = baseline.long.diagLC,
    aes(x = population, y = value, label = "dLC"),
    color = "grey20",
    fontface = "bold",
    hjust = -0.3,
    vjust = 0.1,
    size = 5
  )+
  scale_color_discrete(
    labels = parse(text = c(
    "(lambda[x] * ',' ~ lambda[c])",
    "(10*lambda[x] * ',' ~ 10*lambda[c])",
    "(10*lambda[x] * ',' ~ 0.1*lambda[c])",
    "(0.1*lambda[x] * ',' ~ 10*lambda[c])",
    "(0.1*lambda[x] * ',' ~ 0.1*lambda[c])"
  )))+
  theme_bw(base_size = 18) +
  labs(x=NULL,y=NULL,color="Smoothing pars",shape=expression(delta[c])) +
  theme(strip.text = element_text(face="bold"),
    legend.text=element_text(),
    panel.grid.minor = element_blank()
  )


## saving Figure
ggsave(file="figs/F5.pdf",width = 12,height=6)

## END
