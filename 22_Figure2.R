## --------------------------------------------------------- ##
##
##  FILE 22: plotting Figure 2
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
library(ggtext)

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

## colors for diagonals
my.cols <- c("#D55E00","#2979FF")

## 
shade_df <- data.frame(
  Year = c(2019,
           2100,
           2100,
           2019),
  Age  = c(-Inf,
           2100-2019,
           -Inf,
           -Inf)
)

## visualizing data
Period.all %>%
  mutate(Rate=Deaths/Exposures) %>% 
  filter(Sex=="F") %>% 
  filter(Country=="USA") %>% 
  mutate(Country=case_when(
    Country=="SWE"~"Sweden",
    TRUE ~ Country
  )) %>% 
  mutate(brks=cut(log(Rate), breaks=brek, labels=labe)) %>% 
  ggplot(aes(Year,Age)) + 
  geom_polygon(
    data = shade_df,aes(x = Year, y = Age),
    fill = "grey90",alpha = 0.9,
    inherit.aes = FALSE
  ) +
  geom_tile(aes(fill=brks))+
  scale_fill_viridis("Death Rates", discrete = T,direction = -1,
                     guide = guide_legend(reverse=TRUE),
                     na.value = "grey90") +
  geom_vline(xintercept = seq(min(Period.all$Year)+10, 2100, 10), col="grey50", lty=3, lwd=0.2)+
  geom_hline(yintercept = seq(10, 100, 10), col="grey50", lty=3, lwd=0.2)+
  geom_abline(slope = 1, intercept = seq(-4000, 2000,10), 
              col="grey60", lty=3, lwd=0.2)+
  geom_abline(slope = 1, intercept = c(-1950,-2019),
              col=my.cols, lty=1, lwd=1)+
  geom_abline(slope = 1, intercept = c(-1950.5,-2019.5),
              col=my.cols, lty=1, lwd=2)+
  geom_abline(slope = 1, intercept = c(-1951,-2020),
              col=my.cols, lty=1, lwd=1)+
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(seq(min(Period.all$Year), 1990, 20),
               2005,
               max(Period.all$Year),
               seq(2040, 2100, 20)),
    labels = function(x) {
      ifelse(x == 1950,
             paste0("<span style='color:", my.cols[1], "'>", x, "</span>"),
             ifelse(x == 2019,
                    paste0("<span style='color:", my.cols[2], "'>", x, "</span>"),
                    x))
    }
  ) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(0, 100, 10)), labels = c(seq(0, 90, 10),"100+"))+
  ## text
  ## text range
  geom_text(data = data.frame(Year = 2052,Age = 70,lab = "Text"),
            label = "Useful \nforecasts",size=10) +
  geom_text(data = data.frame(Year = 2075,Age = 22,lab = "Text"),
            label = "Unnecessary \nforecasts",size=10) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        panel.background = element_rect(fill = "grey70", colour = NA),
        axis.text.x = ggtext::element_markdown(),  # allows colored labels
        legend.title = element_text(size=16),
        panel.spacing.x=unit(1, "lines"))+
  labs(y="Age",x="Year") 

ggsave(file="figs/F2.pdf",width = 12,height=7)

## END





