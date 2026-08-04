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

## observations
print(paste(df.obs$obs, collapse = " & "))

## effective dimensions
df.ed %>% 
  filter(model=="LC") %>% pull(ed) %>% paste(., collapse = " & ")
df.ed %>% 
  filter(model=="LLC") %>% pull(ed) %>% paste(., collapse = " & ")
df.ed %>% 
  filter(model=="CCP") %>% pull(ed) %>% paste(., collapse = " & ")

## here for RMSE, CPD and DSS
tab <- cbind(
  Australia_F = res[1,],
  Australia_M = res[2,],
  France_F    = res[3,],
  France_M    = res[4,],
  Sweden_F    = res[5,],
  Sweden_M    = res[6,],
  USA_F       = res[7,],
  USA_M       = res[8,]
)

for (i in 1:nrow(tab)){
  print(paste(tab[i,], collapse = " & "))
}

## compare with Table in manuscript
row.names <- c(
  "RMSE dLC",
  "RMSE dCP",
  "RMSE LC",
  "RMSE LLC",
  "RMSE CCP",
  "CPD dLC",
  "CPD dCP",
  "CPD LC",
  "CPD LLC",
  "CPD CCP",
  "DSS dLC",
  "DSS dCP",
  "DSS LC",
  "DSS LLC",
  "DSS CCP"
)

tab <- cbind(
  Row = row.names,
  Australia_F = res[1,],
  Australia_M = res[2,],
  France_F    = res[3,],
  France_M    = res[4,],
  Sweden_F    = res[5,],
  Sweden_M    = res[6,],
  USA_F       = res[7,],
  USA_M       = res[8,]
)

tab

## finding the minimum
df.res %>% 
  pivot_longer(-c(cou,sex,model)) %>% 
  group_by(cou,sex,name) %>% 
  summarise(min   = min(value,na.rm = T),
            model = model[which.min(value)]) %>% print(n=100)

## END




