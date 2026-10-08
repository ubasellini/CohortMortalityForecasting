
## --------------------------------------------------------- ##
##
##  FILE 00: install useful packages from CRAN and archive
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
## --------------------------------------------------------- ##


## available on CRAN
install.packages("tidyverse")
install.packages("viridis")
install.packages("HMDHFDplus")
install.packages("forecast")
install.packages("MASS")
install.packages("ggplot2")
install.packages("ggtext")
install.packages("patchwork")
install.packages("ggrepel")
install.packages("devtools")

## install archived packages
devtools::install_version('svcm', version = '0.1.2')
devtools::install_version('MortalitySmooth', version = '2.3.4')