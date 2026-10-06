## --------------------------------------------------------- ##
##
##  FILE 44: plotting Figure S4
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

## seed(1) CCP
baseline.long.CCP <- df.res %>%
  filter(model=="CCP") %>% 
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

## second best model
baseline.long.second <- df.res %>%
  filter(model!="CCP") %>% 
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
  group_by(population,measure) %>% 
  summarise(value=min(value)) %>% 
  mutate(measure=factor(measure,levels=c("rmse","cpd","dss"),
                        labels=c("RMSE","CPD 80%","DSS")))



## loading previously saved data 
load(file = "results/16_out10y_CCP_seed.Rdata")

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



## plotting
df.plot %>%
  ggplot(aes(x=population,y=value)) +
  geom_boxplot(
    width=0.6,
    alpha=0,
    color="grey60",
    outlier.shape = NA
  ) +
  # individual sensitivity runs
  geom_jitter(
    width=0.12,size=2.5,alpha=0.75) +
  ## second best
  geom_point(
    data = baseline.long.second,
    aes(x = population, y = value),
    size = 4,
    shape = 23,
    # fill = "white",
    color = "#ff7f00",
    stroke = 1.5
  )+
  geom_text(
    data = baseline.long.second,
    aes(x = population, y = value, label = "2nd best"),
    color = "#ff7f00",
    fontface = "bold",
    hjust = -0.25,
    vjust = 0.1,
    size = 4
  )+
  ## seed (1)
  geom_point(
    data = baseline.long.CCP,
    aes(x = population, y = value),
    size = 4,
    shape = 23,
    # fill = "white",
    color = "#1b9e77",
    stroke = 1.5
  )+
  geom_text(
    data = baseline.long.CCP,
    aes(x = population, y = value, label = "CCP_s(1)"),
    color = "#1b9e77",
    fontface = "bold",
    hjust = -0.2,
    vjust = 0.1,
    size = 4
  )+
  facet_wrap(~measure,nrow=1,scales="free_y") +
  theme_bw(base_size = 18) +
  labs(x=NULL,y=NULL)

## saving Figure
ggsave(file="figs/FS4.pdf",width = 12,height=6)

## END
