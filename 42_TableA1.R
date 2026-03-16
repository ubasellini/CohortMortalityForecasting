## --------------------------------------------------------- ##
##
##  FILE 42: tabulating Table A1
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
load(file = "results/12_out20y.Rdata")

## printing values for manuscript table
for (i in 1:nrow(res)){
  print(paste(res[i,], collapse = " & "))
}

## finding the minimum
df.res %>% 
  pivot_longer(-c(cou,sex,model)) %>% 
  group_by(cou,sex,name) %>% 
  summarise(min   = min(value,na.rm = T),
            model = model[which.min(value)]) %>% print(n=100)

## END

