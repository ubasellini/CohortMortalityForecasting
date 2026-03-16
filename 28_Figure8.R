## --------------------------------------------------------- ##
##
##  FILE 28: plotting Figure 8
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
library(patchwork)

## start by finding number of available ages for each cohort and population
## all countries
all.cou <- c("AUS","FRATNP","SWE","USA")
c <- 1850:2019
df.full.years <- expand.grid(country=all.cou,cohort=c,KEEP.OUT.ATTRS = F)
df.full.years$ages <- NA
i <- 1
for (i in 1:length(all.cou)){
  ## select country
  cou <- all.cou[i]
  ## loading data
  load(file=paste0("data/input/",cou,"_F",".Rdata"))
  j <- 1
  for (j in 1:nc){
    sum(!is.na(cE[,j]))
    df.full.years <- df.full.years %>% 
      mutate(ages=ifelse(country==cou & cohort==c[j], sum(!is.na(cE[,j])),ages))
  }
}

## modify tibble names and column type
df.full.years <- df.full.years %>% 
  rename(cou=country) %>% 
  mutate(cou=as.character(cou))

## remove everything expect new tibble
rm(list=setdiff(ls(), "df.full.years"))
  
## loading main results
load(file="results/13_main_results.Rdata")


## ---- plotting  -----

## for plotting, set to NAs the fitted value of dLC 
## (as they are not the fitted ones but rather the observed ones)
df.plot <- df.e0 %>% left_join(df.full.years) %>% 
  # mutate(period=100+cohort,
  #        forecast=ifelse(period > 2019, 1, 0),
  #        dLC_med=ifelse(forecast==1,dLC_med,NA),
  #        dLC_upp=ifelse(forecast==1,dLC_upp,NA),
  #        dLC_low=ifelse(forecast==1,dLC_low,NA)) %>% 
  mutate(cou=case_when(
    cou == "AUS" ~ "Australia",
    cou == "FRATNP" ~ "France",
    cou == "SWE" ~ "Sweden",
    cou == "USA" ~ "USA"),
    sex=case_when(
      sex == "F" ~ "Females",
      sex == "M" ~ "Males"),
    sex=factor(sex),cou=factor(cou))

## df for models only
df.model.med <- df.plot %>% 
  dplyr::select(cohort,sex,cou,dLC=dLC_med,CCP=CCP_med) %>% 
  ## remove backcast for Australia and USA
  mutate(CCP=case_when(
    cou=="Australia" & cohort< 1921 ~ NA,
    cou=="USA" & cohort< 1933 ~ NA,
    TRUE ~ CCP)) %>% 
  pivot_longer(-c(cohort,sex,cou),names_to = "model") %>% 
  mutate(model=factor(model))
df.model.ribbon <- df.plot %>%
  dplyr::select(cohort, sex, cou,
         CCP_low, CCP_upp,
         dLC_low, dLC_upp) %>%
  pivot_longer(cols = everything()[-c(1:3)], 
               names_to = c("model", ".value"), 
               names_pattern = "(.*)_(low|upp)") %>%
  rename(lower = low, upper = upp) %>%
  mutate(model = factor(model))

## colors
my.cols <- c("#984ea3","#1b9e77")
my.cols <- c("#e41a1c", "#ff7f00")
my.cols <- c("#1b9e77", "#ff7f00")

## test with raster

# Define y-range of your plot (expand as needed)
y_range <- range(c(df.plot$CCP_low,df.plot$CCP_upp), na.rm = TRUE) 

# Build a grid of (cohort, y) so the raster fills entire plot
c <- unique(df.plot$cohort)
df.bg <- expand_grid(
  cohort = c,
  y = seq(floor(y_range[1]), ceiling(y_range[2]), by = 1)
) %>%
  left_join(df.plot %>% dplyr::select(cohort, ages,cou,sex) %>% distinct(),
            by = "cohort", relationship = "many-to-many")

df.plot %>% 
  ggplot(aes(x=cohort,group=sex,shape=sex))+
  geom_raster(data = df.bg,
              aes(x = cohort, y = y, fill = ages),
              inherit.aes = FALSE) +
  geom_point(aes(y=obs),size=1.5,stroke=1.05)+
  facet_wrap(.~cou)+
  theme_bw(base_size = 22)+
  geom_line(data = df.model.med,
            aes(y=value,color=model,group = interaction(model,sex)),
            linewidth=0.9)+
  geom_ribbon(aes(ymin=CCP_low, ymax=CCP_upp),
              alpha=0.6, fill=my.cols[1]) +
  geom_ribbon(aes(ymin=dLC_low, ymax=dLC_upp),
              alpha=0.6, fill=my.cols[2]) +
  theme_bw(base_size = 22)+
  scale_shape_manual(values=c(15,16))+
  coord_cartesian(expand = FALSE) +   # removes margins
  scale_fill_viridis_c(name = "ages \navailable",option="A", direction = -1, alpha = 0.75) +
  scale_color_manual(values = c("CCP" = my.cols[1], "dLC" = my.cols[2])) +
  labs(y="Life expectancy") 

## find range of le for new plot
le.min <- min(df.plot$obs,na.rm=T)-1.15
le.max <- max(max(df.plot$CCP_upp,na.rm=T),max(df.plot$dLC_upp,na.rm=T))+1.15

f1 <- df.plot %>% 
  filter(cou%in%c("France","Sweden")) %>% 
  ggplot(aes(x=cohort,group=sex,shape=sex))+
  geom_raster(data = df.bg %>% 
                filter(cou%in%c("France","Sweden")),
              aes(x = cohort, y = y, fill = ages),
              inherit.aes = FALSE) +
  geom_point(aes(y=obs),size=1.45,stroke=0.6)+
  facet_wrap(.~cou,ncol = 1)+
  theme_bw(base_size = 22)+
  geom_line(data = df.model.med %>% 
              filter(cou%in%c("France","Sweden")),
            aes(y=value,color=model,group = interaction(model,sex)),
            linewidth=0.9)+
  geom_ribbon(aes(ymin=CCP_low, ymax=CCP_upp),
              alpha=0.6, fill=my.cols[1]) +
  geom_ribbon(aes(ymin=dLC_low, ymax=dLC_upp),
              alpha=0.6, fill=my.cols[2]) +
  theme_bw(base_size = 22)+
  scale_shape_manual(values=c(15,16))+
  coord_cartesian(expand = FALSE) +   # removes margins
  scale_fill_viridis_c(option="A", direction = -1, alpha = 0.75) +
  scale_color_manual(values = c("CCP" = my.cols[1], "dLC" = my.cols[2])) +
  # labs(y="Life expectancy")+
  labs(y=NULL) +
  ylim(c(le.min,le.max))
f1


f2 <- df.plot %>% 
  filter(!cou%in%c("France","Sweden")) %>% 
  ## remove backcast for Australia and USA
  mutate(CCP_low=case_when(
    cou=="Australia" & cohort< 1921 ~ NA,
    cou=="USA" & cohort< 1933 ~ NA,
    TRUE ~ CCP_low),
    CCP_upp=case_when(
      cou=="Australia" & cohort< 1921 ~ NA,
      cou=="USA" & cohort< 1933 ~ NA,
      TRUE ~ CCP_upp)) %>% 
  ggplot(aes(x=cohort,group=sex,shape=sex))+
  geom_raster(data = df.bg %>% 
                filter(!cou%in%c("France","Sweden")),
              aes(x = cohort, y = y, fill = ages),
              inherit.aes = FALSE) +
  geom_point(aes(y=obs),size=1.45,stroke=0.6)+
  facet_wrap(.~cou,ncol = 1)+
  theme_bw(base_size = 22)+
  geom_line(data = df.model.med %>% 
              filter(!cou%in%c("France","Sweden")),
            aes(y=value,color=model,group = interaction(model,sex)),
            linewidth=0.9)+
  geom_ribbon(aes(ymin=CCP_low, ymax=CCP_upp),
              alpha=0.6, fill=my.cols[1]) +
  geom_ribbon(aes(ymin=dLC_low, ymax=dLC_upp),
              alpha=0.6, fill=my.cols[2]) +
  theme_bw(base_size = 22)+
  scale_shape_manual(values=c(15,16))+
  coord_cartesian(expand = FALSE) +   # removes margins
  scale_fill_viridis_c(option="A", direction = -1, alpha = 0.75) +
  scale_color_manual(values = c("CCP" = my.cols[1], "dLC" = my.cols[2])) +
  # labs(y=NULL) +
  labs(y="Life expectancy")+
  xlim(c(1921,2019))+ylim(c(le.min,le.max))
f2
f2 <- f2 + theme(legend.position = "none")

## putting plots together
# f1+f2+patchwork::plot_layout(guides = 'collect', axes = 'collect', widths = c(1.5, 1))
f2+f1+patchwork::plot_layout(axes = 'collect', widths = c(1, 1.5))

## saving Figure
ggsave(file="figs/F8.pdf",width = 12,height=8)


## END
