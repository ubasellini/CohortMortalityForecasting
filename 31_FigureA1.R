## --------------------------------------------------------- ##
##
##  FILE 31: plotting Figure A1
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
##  viridis_0.6.5     viridisLite_0.4.2 lubridate_1.9.3   forcats_1.0.0    
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

## loading previously saved data 
load(file = "results/14_sensitivity.Rdata")

## generate figure for paper
df.res.long <- df.res %>% 
  mutate(par=case_when(
    par=="Alpha"~"alpha[x]",
    par=="Beta"~"beta[x]",
    par=="Kappa"~"kappa[c]"),
    cou=case_when(
      cou=="SWE_F" ~ "Sweden,females",
      cou=="FRATNP_M" ~ "France,males")) %>% 
  pivot_longer(-c(par,cou,scale),names_to = "Year")

## select cohorts to highlight
my.coh <- seq(1880,1950,5)

df.res.long %>% 
  ggplot(aes(x=scale,y=value,group = Year))+
  # geom_line(color="grey70",linewidth=0.75)+
  geom_line(data=df.res.long %>% filter(Year%in%my.coh),aes(color=Year),linewidth=1.15)+
  scale_color_viridis_d()+
  facet_wrap(cou~par,scales = "free",labeller = labeller(.cols = label_parsed, .multi_line = FALSE))+
  theme_bw(base_size = 22)+
  labs(x=NULL,y="Cohort LC parameters",color="First \nYear")

## saving Figure
ggsave(file="figs/FA1.pdf",width = 12,height=10)

## for all years
df.res.long %>% mutate(Year=as.numeric(Year)) %>% 
  ggplot(aes(x=scale,y=value,group = Year,color=Year))+
  geom_line(linewidth=0.85)+
  scale_color_viridis()+
  facet_wrap(cou~par,scales = "free",labeller = labeller(.cols = label_parsed, .multi_line = FALSE))+
  theme_bw(base_size = 22)+
  labs(x=NULL,y="Cohort LC parameters",color="First \nYear")

## saving Figure
# ggsave(file="figs/FA1_2.pdf",width = 12,height=10)

## END
