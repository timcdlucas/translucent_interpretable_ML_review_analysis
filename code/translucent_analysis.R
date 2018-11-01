#'---
#'output:
#'  pdf_document: default
#'title: "Translucent box"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#' ## Setup

#+ knitrsetup, echo = FALSE

knitr::opts_chunk$set(cache = TRUE, fig.width = 7, fig.height = 5)


#+ libs,

library(dplyr)
library(ggplot2)
library(caret)


#'## The data

#' First we need to read in the data.

#+ data_read
p <- read.table(file = '../data/PanTHERIA_1-0_WR05_Aug2008.txt',
  header = TRUE, sep = "\t", na.strings = c("-999", "-999.00"))

names(p)
sapply(p, function(x) mean(is.na(x))) %>% sort
sapply(p, class)


#' Now we need to choose a variable of interest and make some basic exploratory plots.


#+ Choose 

# Maybe social group size.

ggplot(p, aes(X10.2_SocialGrpSize)) + geom_histogram()

p$X10.2_SocialGrpSize %>% summary

large_orders <- 
  p %>% 
    filter(!is.na(X10.2_SocialGrpSize)) %>% 
    group_by(MSW05_Order) %>% 
    count() %>%
    arrange(desc(n)) %>% 
    filter(n > 20) %>% 
    pull(MSW05_Order)

p %>% 
  filter(MSW05_Order %in% large_orders) %>% 
  ggplot(aes(X10.2_SocialGrpSize)) + 
    geom_histogram() + 
    facet_wrap(~ MSW05_Order)


ggplot(p, aes(x = MSW05_Order, y = X10.2_SocialGrpSize)) + geom_boxplot()

## Don't wish to do too many bivariate plots at this point. Going to do a priori variable selection and p values later.
ggplot(p, aes(x = X21.1_PopulationDensity_n.km2, y = X10.2_SocialGrpSize)) + 
  geom_point() + 
  geom_smooth() + 
  geom_smooth(method = 'lm', se = FALSE, colour = 'red') + 
  scale_x_log10()





