## --------------------------------------------------------- ##
##
##  FILE 46: tabulating Table S2
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
load(file = "results/15_out10y_CCP_robustness.Rdata")

## printing values for manuscript table
df.res <- df.res %>%
  mutate(
    lambda.label = case_when(
      lambda.com == "hat" ~ "$(\\hat{\\lambda}_1,\\hat{\\lambda}_2)$",
      lambda.com == "HH"  ~ "$(10\\hat{\\lambda}_1,10\\hat{\\lambda}_2)$",
      lambda.com == "HL"  ~ "$(10\\hat{\\lambda}_1,0.1\\hat{\\lambda}_2)$",
      lambda.com == "LH"  ~ "$(0.1\\hat{\\lambda}_1,10\\hat{\\lambda}_2)$",
      lambda.com == "LL"  ~ "$(0.1\\hat{\\lambda}_1,0.1\\hat{\\lambda}_2)$"
    )
  )
names(df.res)

fr <- subset(df.res, cou=="FRATNP" & sex=="M")
sw <- subset(df.res, cou=="SWE" & sex=="F")
tab.out <- merge(
  fr,
  sw,
  by = c("lambda.com","lambda.label","delta.t","upper"),
  suffixes = c(".FR",".SW")
)

tab.out <- tab.out %>%
  mutate(
    lambda.com = factor(
      lambda.com,
      levels = c("hat","HH","HL","LH","LL")
    ),
    upper = factor(
      upper,
      levels = c(TRUE,FALSE)
    ),
    main.case = ifelse(
      lambda.com=="hat" &
        delta.t==50 &
        upper==TRUE,
      0,
      1
    )
  ) %>%
  arrange(main.case, lambda.com, delta.t, desc(upper)) %>%
  dplyr::select(-main.case)


for(i in seq_len(nrow(tab.out))){
  
  cat(
    tab.out$lambda.label[i], " & ",
    tab.out$delta.t[i], " & ",
    as.character(tab.out$upper[i]), " & ",
    
    sprintf("%.2f",tab.out$rmse.FR[i]), " & ",
    sprintf("%.2f",tab.out$cpd.FR[i]), " & ",
    sprintf("%.2f",tab.out$dss.FR[i]), " & ",
    
    sprintf("%.2f",tab.out$rmse.SW[i]), " & ",
    sprintf("%.2f",tab.out$cpd.SW[i]), " & ",
    sprintf("%.2f",tab.out$dss.SW[i]),
    
    "\\\\\n",
    sep=""
  )
  
}

## checking
df.res

## END

