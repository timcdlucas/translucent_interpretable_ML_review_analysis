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


#+ libs, cache = FALSE

library(dplyr)
library(ggplot2)
library(caret)
library(lime)
library(doParallel)
library(pdp)
library(ICEbox)
library(iml)

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


#+ paralell_setup, cache = FALSE


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

plot(m1_enet)

plotCV(m1_enet)


#+ gp?

gp_gr <- data.frame(sigma = c(0.01, 0.02, 0.04, 0.08, 0.16))
m2_gp <- train(y ~ ., data = p_impute, method = 'gaussprRadial', tuneGrid = gp_gr, trControl = trcntrl, na.action = na.omit)

plot(m2_gp)

plotCV(m2_gp)


#+ ranger, eval = TRUE

rf_gr <- expand.grid(mtry = c(2, 5, 10, 20, 30), splitrule = 'variance', min.node.size = c(5, 10, 20, 50))
m3_rf <- train(y ~ ., data = p_impute, method = 'ranger', tuneGrid = rf_gr, trControl = trcntrl, na.action = na.omit, importance = 'impurity')

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

#' # Find variable importance for each model


#+ surrogate_model


#+ varimp

varImp(m1_enet)
varImp(m2_gp)
varImp(m3_rf)



#' ## Now plot some functional forms

#+ pdp_gest

partial(m2_gp, 
        pred.var = c('X9.1_GestationLen_d'),
        parallel = TRUE, plot = TRUE)

partial(m3_rf, 
        pred.var = c('X9.1_GestationLen_d'),
        parallel = TRUE, plot = TRUE)

#+ pdp_lat

partial(m2_gp, 
        pred.var = c('X26.4_GR_MidRangeLat_dd'),
        parallel = TRUE, plot = TRUE)

partial(m3_rf, 
        pred.var = c('X26.4_GR_MidRangeLat_dd'),
        parallel = TRUE, plot = TRUE)

#+ pdp_pet

partial(m2_gp, 
        pred.var = c('X30.2_PET_Mean_mm'),
        parallel = TRUE, plot = TRUE)

partial(m3_rf, 
        pred.var = c('X30.2_PET_Mean_mm'),
        parallel = TRUE, plot = TRUE)


#+ pdp_temp

partial(m2_gp, 
        pred.var = c('X28.2_Temp_Mean_01degC'),
        parallel = TRUE, plot = TRUE)


partial(m3_rf, 
        pred.var = c('X28.2_Temp_Mean_01degC'),
        parallel = TRUE, plot = TRUE)


#+ pdp_gest_lat

partial(m2_gp, 
        pred.var = c('X9.1_GestationLen_d', 'X26.4_GR_MidRangeLat_dd'),
        parallel = TRUE, plot = TRUE)

partial(m3_rf, 
        pred.var = c('X9.1_GestationLen_d', 'X26.4_GR_MidRangeLat_dd'),
        parallel = TRUE, plot = TRUE)


#+ pdp_gest_pet

partial(m2_gp, 
        pred.var = c('X9.1_GestationLen_d', 'X30.2_PET_Mean_mm'),
        parallel = TRUE, plot = TRUE)

partial(m3_rf, 
        pred.var = c('X9.1_GestationLen_d', 'X30.2_PET_Mean_mm'),
        parallel = TRUE, plot = TRUE)


#' Test for interactions

#+ all_inter

predictor_gp = Predictor$new(m2_gp, data = select(p_impute, -y), y = p_impute$y)
predictor_rf = Predictor$new(m3_rf, data = select(p_impute, -y), y = p_impute$y)


interact_gp = Interaction$new(predictor_gp)
interact_gp
plot(interact_gp)


interact_rf = Interaction$new(predictor_rf)
interact_rf
plot(interact_rf)


#+ which_inter

interact_gp = Interaction$new(predictor_gp, feature = "X9.1_GestationLen_d")
plot(interact_gp)


interact_rf = Interaction$new(predictor_rf, feature = "X9.1_GestationLen_d")
plot(interact_rf)

interact_gp2 = Interaction$new(predictor_gp, feature = "X26.4_GR_MidRangeLat_dd")
plot(interact_gp2)


interact_rf2 = Interaction$new(predictor_rf, feature = "X26.4_GR_MidRangeLat_dd")
plot(interact_rf2)

#' # ICE plots

#+ ice

m2_gp_ice <- ice(m2_gp, p_impute, p_impute$y, 'X9.1_GestationLen_d', frac_to_build = 0.1)
plot(m2_gp_ice)

m3_rf_ice <- ice(m3_rf, p_impute, p_impute$y, 'X9.1_GestationLen_d', frac_to_build = 0.1)
plot(m3_rf_ice)


#+ ice_lat

m2_gp_ice_lat <- ice(m2_gp, p_impute, p_impute$y, 'X26.4_GR_MidRangeLat_dd', frac_to_build = 0.1)
plot(m2_gp_ice_lat)

m3_rf_ice_lat <- ice(m3_rf, p_impute, p_impute$y, 'X26.4_GR_MidRangeLat_dd', frac_to_build = 0.1)
plot(m3_rf_ice_lat)


#+ ice_temp

m2_gp_ice_temp <- ice(m2_gp, p_impute, p_impute$y, 'X28.2_Temp_Mean_01degC', frac_to_build = 0.1)
plot(m2_gp_ice_temp)

m3_rf_ice_temp <- ice(m3_rf, p_impute, p_impute$y, 'X28.2_Temp_Mean_01degC', frac_to_build = 0.1)
plot(m3_rf_ice_temp)


#+ clustered_ice


m2_gp_ice_c <- ice(m2_gp, p_impute, p_impute$y, 'X9.1_GestationLen_d')
clusterICE(m2_gp_ice_c, nClusters = 20, centered = TRUE)
clusterICE(m2_gp_ice_c, nClusters = 20, centered = FALSE)


m3_rf_ice_c <- ice(m3_rf, p_impute, p_impute$y, 'X9.1_GestationLen_d')
clusterICE(m3_rf_ice_c, nClusters = 20, centered = FALSE)



#+ clustered_ice_lat


m2_gp_ice_lat_c <- ice(m2_gp, p_impute, p_impute$y, 'X26.4_GR_MidRangeLat_dd')
clusterICE(m2_gp_ice_lat_c, nClusters = 20, centered = FALSE)


m3_rf_ice_lat_c <- ice(m3_rf, p_impute, p_impute$y, 'X26.4_GR_MidRangeLat_dd')
clusterICE(m3_rf_ice_lat_c, nClusters = 20, centered = FALSE)




#+ clustered_ice_temp


m2_gp_ice_temp_c <- ice(m2_gp, p_impute, p_impute$y, 'X28.2_Temp_Mean_01degC')
clusterICE(m2_gp_ice_temp_c, nClusters = 20, centered = FALSE)


m3_rf_ice_temp_c <- ice(m3_rf, p_impute, p_impute$y, 'X28.2_Temp_Mean_01degC')
clusterICE(m3_rf_ice_temp_c, nClusters = 20, centered = FALSE)



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


#+ lime, fig.width = 10, fig.height = 8

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











#+ session_info

sessionInfo()

