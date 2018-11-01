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

# Want something with quite a lot of data. Litter size?

sum(!is.na(p$X15.1_LitterSize))

ggplot(p, aes(X15.1_LitterSize)) + geom_histogram()

p$X15.1_LitterSize %>% summary

large_orders <- 
  p %>% 
    filter(!is.na(sum(!is.na(p$X15.1_LitterSize)))) %>% 
    group_by(MSW05_Order) %>% 
    count() %>%
    arrange(desc(n)) %>% 
    filter(n > 40) %>% 
    pull(MSW05_Order)

p %>% 
  filter(MSW05_Order %in% large_orders) %>% 
  ggplot(aes(X15.1_LitterSize)) + 
    geom_histogram() + 
    facet_wrap(~ MSW05_Order)


ggplot(p, aes(x = MSW05_Order, y = X15.1_LitterSize)) + geom_boxplot()

## Don't wish to do too many bivariate plots at this point. Going to do a priori variable selection and p values later.
ggplot(p, aes(x = X21.1_PopulationDensity_n.km2, y = X15.1_LitterSize)) + 
  geom_point() + 
  geom_smooth() + 
  geom_smooth(method = 'lm', se = FALSE, colour = 'red') + 
  scale_x_log10()


#' ## Some data cleaning


#' # Now fit 4 models.

#+ a_priori_var_selection



#+ elastic_net



#+ 


3 models
typical ecological modelling
regularised regression
  - maxent is regularised regression
  - linear can still have interactions,
 nonlinear, etc.
  - relationship to Bayes?
spline or gp
random forest/xgboost


2 or 3 questions.
  - generate hypotheses
    - var imp
    - interaction importance
    - ice etc.
  - gain some understanding of a system
     - predictability
     - complexity
     - r2
     - random effects?
  - understand individual points
    - lime




