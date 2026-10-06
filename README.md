# Forecasting Cohort Mortality: Lee-Carter Methods and _CCP_-Splines

This repository contains materials to fully reproduce the results of the conditionally accepted paper:

Basellini U. and C.G. Camarda. Forecasting Cohort Mortality: Lee—Carter Methods and _CCP_-Splines. _International Journal of Forecasting_


### Repository structure

- _data_ : input data 
- _figs_ : output figures presented in the manuscript
- _funs_ : functions used by the codes (including functions to fit the mortality models presented in the manuscript) 
- _results_ : output of the analysis (to be used to generate the figures and tables in the manuscript)
- codes to replicate analysis and visualisations:
  - codes starting with 0: packages installation and data download (need HMD credentials)
  - codes starting with 1: analysis
  - codes starting with 2 and 3: visualizations and tables of main manuscript
  - codes starting with 4: visualizations and tables of Supplementary Materials
  
### Computing environment  

All analyses are performed using the `R` statistical software (R Core Team, 2025). The packages needed to run the codes can be installed from the file `00_install-packages.R`. sessionInfo() details are provided in each `R` file. Analyses were carried using on a portable personal computer with an Intel i5-1245U processor, 1.6 GHz, and 16 GB of RAM (please check the manuscript for computational times). 

R Core Team (2025). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>


### Data 

All data can be freely downloaded (upon registration) from the Human Mortality Database (<https://www.mortality.org/>). Here, the file `01_downolad-data.R` downloads the data and stores them in the _data_ folder.

  
### Authors' contact information
  
- Ugofilippo Basellini: <basellini@demogr.mpg.de>
- Carlo G. Camarda: <carlo-giovanni.camarda@ined.fr>



This repository was last updated on October 6, 2026.
