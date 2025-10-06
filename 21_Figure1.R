## --------------------------------------------------------- ##
##
##  FILE 21: plotting Figure 1
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

## loading data
load(file=paste0("data/input/all_countries",".Rdata"))

## manually choosing the breaks
mybrek <- c(0.0003,0.0005,0.0008,0.0015,0.0020,0.0030,0.0060,0.0150,0.0525)
brek <- c(-20, log(mybrek), 20)
labe <- c("<0.0003",
          "0.0003-0.0005",
          "0.0005-0.0008",
          "0.0008-0.0015",
          "0.0015-0.0020",
          "0.0020-0.0030",
          "0.0030-0.0060",
          "0.0060-0.0150",
          "0.0150-0.0525",
          ">0.0525")

## visualizing data
Cohort.all %>%
  filter(Sex=="F") %>% 
  filter(Country=="SWE"|Country=="USA") %>% 
  mutate(Country=case_when(
    Country=="SWE"~"Sweden",
    TRUE ~ Country
  )) %>% 
  mutate(brks=cut(log(Rate), breaks=brek, labels=labe)) %>% 
  ggplot(aes(Year,Age)) + 
  geom_tile(aes(fill=brks))+
  scale_fill_viridis("", discrete = T,direction = -1,
                     guide = guide_legend(reverse=TRUE),
                     na.value = "grey90") +
  geom_vline(xintercept = seq(min(Cohort.all$Year)+10, max(Cohort.all$Year), 10), col="grey50", lty=3, lwd=0.2)+
  geom_hline(yintercept = seq(10, 100, 10), col="grey50", lty=3, lwd=0.2)+
  geom_abline(slope = -1, intercept = seq(1810, 2100,10), 
              col="grey60", lty=3, lwd=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(seq(min(Cohort.all$Year), 2000,20),2005,max(Cohort.all$Year))) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 100, 10)), labels = c(seq(0, 90, 10),"100+"))+
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.text.x=element_text(angle=45,hjust = 1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        panel.spacing.x=unit(1, "lines"))+
  labs(y="Age",x="Cohort") +
  facet_wrap(.~Country)

ggsave(file="figs/F1.pdf",width = 12,height=7)

## END