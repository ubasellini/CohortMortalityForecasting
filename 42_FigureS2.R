## --------------------------------------------------------- ##
##
##  FILE 42: plotting Figure S2
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

## age and cohorts
x <- 0:100
c <- 1850:2019
m <- length(x)
nc <- length(c)
AGES <- matrix(x,nrow=m,ncol=nc)
COHORTS <- matrix(rep(c,each=m),nrow=m)
PERIODS <- COHORTS + AGES
last.t <- max(c)

##---- first out-of-sample scenario (10y)
h.out <- 10
DF <- expand.grid(ages=x, cohorts=c) %>% 
  mutate(periods=c(PERIODS),
         dummy=case_when(
           periods > last.t ~ NA,
           periods > last.t - h.out & periods <= last.t ~ 0,
           TRUE           ~ 1),
         dummy=ifelse(cohorts>last.t - h.out,NA,dummy),
         dummy=as.factor(dummy),
         year = h.out)

f1 <- DF %>%
  ggplot(aes(cohorts,ages)) + 
  geom_tile(aes(fill=dummy))+
  scale_fill_viridis("", discrete = T,direction = -1,
                    guide = guide_legend(reverse=TRUE),
                    begin = 0.2, end = 0.8,
                    na.value = "grey90",na.translate = F,
                    labels = c("Test", "Training")) +
  geom_vline(xintercept = c(seq(min(c)+10, 2000, 10),last.t - h.out), col="grey50", lty=3, lwd=0.2)+
  geom_hline(yintercept = seq(10, 100, 10), col="grey50", lty=3, lwd=0.2)+
  geom_abline(slope = -1, intercept = seq(1810, 2100,10), 
              col="grey60", lty=3, lwd=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(seq(min(c), 2000,20),last.t - h.out,max(c))) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 100, 10)), labels = c(seq(0, 90, 10),"100+"))+
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.text.x=element_text(angle=45,hjust = 1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        panel.spacing.x=unit(1, "lines"))+
  labs(y="Age",x="Cohort") 

## second exercise
h.out <- 20
DF1 <- expand.grid(ages=x, cohorts=c) %>% 
  mutate(periods=c(PERIODS),
         dummy=case_when(
           periods > last.t ~ NA,
           periods > last.t - h.out & periods <= last.t ~ 0,
           TRUE           ~ 1),
         dummy=ifelse(cohorts>last.t - h.out,NA,dummy),
         dummy=as.factor(dummy),
         year = h.out)

## merging
DF0 <- DF %>% 
  bind_rows(DF1) %>% 
  mutate(year=as.factor(year))

year_names <- c(
  `10` = "10y validation",
  `20` = "20y validation"
)

DF0 %>%
  ggplot(aes(cohorts,ages)) + 
  geom_tile(aes(fill=dummy))+
  scale_fill_viridis("", discrete = T,direction = -1,
                     guide = guide_legend(reverse=TRUE),
                     begin = 0.2, end = 0.8,
                     na.value = "grey90",na.translate = F,
                     labels = c("Test", "Training")) +
  facet_grid(.~year,labeller = as_labeller(year_names))+
  geom_vline(xintercept = c(seq(min(c)+10, 1990, 10),1999,2009,2019), col="grey50", lty=3, lwd=0.2)+
  geom_hline(yintercept = seq(10, 100, 10), col="grey50", lty=3, lwd=0.2)+
  geom_abline(slope = -1, intercept = c(seq(1810, 1990,10),seq(1999,2159,10)), 
              col="grey60", lty=3, lwd=0.2)+
  scale_x_continuous(expand=c(0,0), breaks=c(seq(min(c), 1990,20),max(c)-20,max(c)-10,max(c))) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 100, 10)), labels = c(seq(0, 90, 10),"100+"))+
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.text.x=element_text(angle=45,hjust = 1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        panel.spacing.x=unit(1, "lines"), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  labs(y="Age",x="Cohort")

## saving Figure
ggsave("figs/FS2.pdf",width = 9, height = 4)

## counting the number of green cells in each exercise
DF %>% mutate(oos=case_when(
  dummy == 1 ~ 0,
  dummy == 0 ~ 1,
  is.na(dummy) ~ NA)) %>% 
  select(oos) %>% 
  pull() %>% as.numeric(.) %>% sum(.,na.rm = T)

DF1 %>% mutate(oos=case_when(
  dummy == 1 ~ 0,
  dummy == 0 ~ 1,
  is.na(dummy) ~ NA)) %>% 
  select(oos) %>% 
  pull() %>% as.numeric(.) %>% sum(.,na.rm = T)


## END
