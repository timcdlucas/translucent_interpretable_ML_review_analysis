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
library(lime)
library(doParallel)

source('helpers.R')

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

# Remove NAs in response and response where litter size is less than one (doesn't make sense).
p <- p %>% 
       filter(!is.na(X15.1_LitterSize)) %>% 
       filter(X15.1_LitterSize >= 1) %>% 
       mutate(y = log1p(X15.1_LitterSize)) %>% 
       select(-X15.1_LitterSize, -References, -X24.1_TeatNumber)

p_notaxa <- p %>% 
              select(-contains('MSW05'))

preprocesses <- preProcess(p, method = 'medianImpute')
p_impute <- predict(preprocesses, p_notaxa)


#' # Now fit 4 models.

#+ caret_Setup

folds <- createFolds(p$y, k = 5, returnTrain = TRUE)
trcntrl <- trainControl(index = folds, savePredictions = TRUE, search = 'random')


#+ paralell_setup


cl <- makePSOCKcluster(6)
registerDoParallel(cl)



#+ a_priori_var_selection
apriori_formula <- y ~ X5.1_AdultBodyMass_g + X3.1_AgeatFirstBirth_d + X18.1_BasalMetRate_mLO2hr + X9.1_GestationLen_d + X14.1_InterbirthInterval_d + X16.1_LittersPerYear + X17.1_MaxLongevity_m
m0_lm <- train(apriori_formula, data = p_impute, method = 'lm', trControl = trcntrl, na.action = na.omit)

plotCV(m0_lm)


summary(m0_lm$finalModel)



m0_pglm <- NULL

#+ elastic_net

enet_gr <- expand.grid(lambda = 10 ^ seq(0, -4, length.out = 20), fraction = c(seq(0.01, 1, length.out = 25)))
m1_enet <- train(y ~ ., data = p_impute, method = 'enet', tuneGrid = enet_gr, trControl = trcntrl, na.action = na.omit)

ggplot(m1_enet)

plotCV(m1_enet)


#+ gp?

gp_gr <- data.frame(sigma = c(0.01, 0.02, 0.04, 0.08, 0.16))
m2_gp <- train(y ~ ., data = p_impute, method = 'gaussprRadial', tuneGrid = gp_gr, trControl = trcntrl, na.action = na.omit)

plot(m2_gp)

plotCV(m2_gp)


#+ ranger

m3_rf <- train(y ~ ., data = p_impute, method = 'ranger', tuneLength = 15, trControl = trcntrl, na.action = na.omit)

plot(m3_rf)

plotCV(m3_rf)



#+ any_extras

compare_models(m2_gp, m3_rf)
compare_models(m1_enet, m3_rf)
compare_models(m0_lm, m3_rf)


#' # Global analysis.

#'  - gain some understanding of a system
#'     - predictability
#'     - complexity
#'     - r2

#' # Generate hypotheses to test more formally.


#'  - generate hypotheses (variable level)
#'    - var imp
#'    - interaction importance
#'    - ice etc.


#+ varimp

varImp(m1_enet)
varImp(m2_gp)
varImp(m3_rf)

#' # Examine correlation structure

#'  - examine correlation structure (variable level)
#'     - random effects
#'     - mean zero
#'     - unseen values
#'     - control Vs use in prediction
#'     - shared power
#'     - feature engineering (eg t-1 value)
#'     - stacked generalisation
#'     - priors
#'     - random slopes are regularised interactions. dealt with fairly natively by RF.
#'     - RF on distance to points could be applied to phylogeny.


#' # Point level understanding

#'  - understand individual points
#'    - lime


#+ lime

# Explain Tenrec ecaudatus. Largest litter size. And other large litters
l <- order(p_impute$y, decreasing = TRUE)[1:5]



m1_lime <- lime(p_impute, m1_enet)
m1_explain <- explain(p_impute[l, ], m1_lime, n_features = 6)

plot_features(m1_explain)


m2_lime <- lime(p_impute, m2_gp)
m2_explain <- explain(p_impute[l, ], m2_lime, n_features = 6)

plot_features(m2_explain)


m3_lime <- lime(p_impute, m3_rf)
m3_explain <- explain(p_impute[l, ], m3_lime, n_features = 6)

plot_features(m1_explain)





