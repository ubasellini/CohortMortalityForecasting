## --------------------------------------------------------- ##
##
##  FILE 43: plotting Figure S3
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
##  lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
##  readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0      
##
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading packages
library(tidyverse)

## loading previously saved data 
load(file = "results/11_out10y.Rdata")


df.rmse.age %>% 
  filter((cou=="SWE" & sex=="F" )|(cou=="FRATNP" & sex=="M" )) %>% 
  mutate(label=case_when(
    cou=="SWE" ~ "Sweden, females",
    cou=="FRATNP" ~ "France, males"
  )) %>% 
  ggplot(aes(x=cohorts,y=ages))+
  geom_tile(aes(fill=diff))+
  facet_grid(.~label)+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0)+
  scale_x_continuous(limits = c(1910,2010), expand=c(0,0), breaks=c(seq(min(df.rmse.age$cohorts),1985,25),max(df.rmse.age$cohorts))) +
  scale_y_continuous(limits = c(0, 100),expand=c(0,0), breaks=c(seq(0, 100, 10)), labels = c(seq(0, 90, 10),"100+"))+
  theme_bw(base_size=22) +
  labs(x="Cohort",y="Age",fill=expression(Delta*SE))+
  theme(
    axis.text.x=element_text(angle=45,hjust = 1),
    panel.background = element_rect(fill = "grey85", colour = NA),
    panel.spacing.x=unit(1, "lines")
  ) 

## saving Figure
ggsave(file="figs/FS3.pdf",width = 12,height=8)


## END

