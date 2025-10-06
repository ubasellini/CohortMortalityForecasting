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

## loading main results
load(file="results/13_main_results.Rdata")

## ---- plotting  -----

## select cohorts
my.coh <- c(1950,2019)

## transform everything in rates scale
df.mx <- df.lmx %>% 
  mutate(obs=exp(obs),
         dLC_med=exp(dLC_med),dLC_upp=exp(dLC_upp),dLC_low=exp(dLC_low),
         CCP_med=exp(CCP_med),CCP_upp=exp(CCP_upp),CCP_low=exp(CCP_low))

## for plotting, set to NAs the fitted value of dLC 
## (as they are not the fitted ones but rather the observed ones)
df.plot <- df.mx %>% 
  filter(cohort%in%my.coh) %>% 
  mutate(cou=case_when(
           cou == "AUS" ~ "Australia",
           cou == "FRATNP" ~ "France",
           cou == "SWE" ~ "Sweden",
           cou == "USA" ~ "USA"),
         sex=case_when(
           sex == "F" ~ "Females",
           sex == "M" ~ "Males"),
         cohort=case_when(
           cohort == 1950 ~ "Cohort 1950",
           cohort == 2019 ~ "Cohort 2019"),
         cohort=factor(cohort),sex=factor(sex),cou=factor(cou))

## df for models only
df.model.med <- df.plot %>% 
  dplyr::select(age,cohort,sex,cou,dLC=dLC_med,CCP=CCP_med) %>% 
  pivot_longer(-c(age,cohort,sex,cou),names_to = "model") %>% 
  mutate(model=factor(model))
df.model.ribbon <- df.plot %>%
  dplyr::select(age, cohort, sex, cou,
         CCP_low, CCP_upp,
         dLC_low, dLC_upp) %>%
  pivot_longer(cols = everything()[-c(1:4)], 
               names_to = c("model", ".value"), 
               names_pattern = "(.*)_(low|upp)") %>%
  rename(lower = low, upper = upp) %>%
  mutate(model = factor(model))


my.cols <- c("#984ea3","#1b9e77")

df.plot %>% 
  ggplot(aes(x=age,group=cohort))+
  geom_point(aes(y=obs,shape=cohort),size=1.75)+
  facet_grid(sex~cou)+
  geom_line(data = df.model.med,aes(y=value,color=model,group = interaction(model,cohort)),
            linewidth=0.8)+
  geom_ribbon(data = df.model.ribbon,
              aes(ymin = lower, ymax = upper,
                  fill = model, group = interaction(model, cohort)),
              alpha = 0.3)+
  theme_bw(base_size = 22)+
  labs(shape="Observed",y="death rates")+
  scale_color_manual(values = c("CCP" = my.cols[1], "dLC" = my.cols[2])) +
  scale_fill_manual(values = c("CCP" = my.cols[1], "dLC" = my.cols[2]))+
  scale_y_log10()+
  theme(axis.text.x= element_text(size=16))+
  guides(
    shape = guide_legend(order = 1, title = "Observed"),
    color = guide_legend(order = 2, title = "Model"),
    fill = guide_legend(order = 2, title = "Model")
  )

## saving Figure
ggsave(file="figs/F5.pdf",width = 12,height=8)

## END