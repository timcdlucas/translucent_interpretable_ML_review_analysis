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

# General
library(dplyr)
library(ggplot2)
library(doParallel)

# Modelling libraries
library(caret)
library(lime)
library(pdp)
library(ICEbox)
library(iml)
library(elasticnet)
library(INLA)


# Phylogenetic libraries.
library(ape)
library(caper)

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
       dplyr::select(-X15.1_LitterSize, -References, -X24.1_TeatNumber)


p_notaxa <- p %>% 
              dplyr::select(-contains('MSW05'))

preprocesses <- preProcess(p, method = 'medianImpute')
p_impute <- predict(preprocesses, p_notaxa)



#' ## Get phylogeny data.

#' ### read in phylogeny data.

# Read in trees
tree <- read.nexus('https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.1461-0248.2009.01307.x&file=ELE_1307_sm_SA1.tre')

# Select best supported tree
tree <- tree[[1]]
tree$tip.label <- gsub('_', ' ', tree$tip.label)

# Check if species are available.
mean(p$MSW05_Binomial %in% tree$tip.label)
in_phylo <- p$MSW05_Binomial %in% tree$tip.label

# Remove data that is not in the phylogeny.

p <- p %>% filter(in_phylo) 
p_impute <- p_impute %>% filter(in_phylo)


#' ## Data summary

#+data_summary

dim(p_impute)
names(p_impute)

par(mfrow = c(2, 2))

for(i in 0:11){
  for( j in 1:4){

    if(j + 4 * i <= ncol(p_impute)){
      hist(p_impute[, j + 4 * i], breaks = 100, ylim = c(0, 80), main = j + 4 * i)
    }

  }
  print(i)
  par(mfrow = c(2, 2))
}



#' # Log Transforms
#+ logtrans

log_cols <- c(2, 4, 7, 8, 
              10, 11, 13, 14, 15, 17, 18, 19, 
              20, 21, 22, 23, 26, 27, 28, 29,  
              31, 32, 33, 
              40, 41, 42)

p_impute[, log_cols] <- log1p(p_impute[, log_cols])


par(mfrow = c(2, 2))

for(i in 0:11){
  for( j in 1:4){

    if(j + 4 * i <= ncol(p_impute)){
      hist(p_impute[, j + 4 * i], breaks = 100, ylim = c(0, 80), main = j + 4 * i)
    }

  }
  print(i)
  par(mfrow = c(2, 2))
}


#' # Now fit 4 models.

#+ caret_Setup

folds <- createFolds(p$y, k = 5, returnTrain = TRUE)
trcntrl <- trainControl(index = folds, savePredictions = TRUE, search = 'random')


#+ paralell_setup, cache = FALSE


cl <- makeForkCluster(6)
registerDoParallel(cl)



#+ a_priori_var_selection
apriori_formula <- y ~ X5.1_AdultBodyMass_g + X3.1_AgeatFirstBirth_d + X18.1_BasalMetRate_mLO2hr + X9.1_GestationLen_d + X16.1_LittersPerYear + X17.1_MaxLongevity_m
m0_lm <- train(apriori_formula, data = p_impute, method = 'lm', trControl = trcntrl, na.action = na.omit)

plotCV(m0_lm)

m0_lm


summary(m0_lm$finalModel)



m0_pglm <- NULL

#+ elastic_net

enet_gr <- expand.grid(lambda = 10 ^ seq(0, -4, length.out = 20), fraction = c(seq(0.01, 1, length.out = 25)))
m1_enet <- train(y ~ ., data = p_impute, method = 'enet', tuneGrid = enet_gr, trControl = trcntrl, na.action = na.omit)

plot(m1_enet)

plotCV(m1_enet)



# Find the final parameters 
#   https://stackoverflow.com/questions/40088228/how-to-retrieve-elastic-net-coefficients
#final_pars <- predict.enet(m1_enet$finalModel, s=m1_enet$bestTune[1, "fraction"], type="coef", mode="fraction")$coefficients
#final_pars
#sum(final_pars != 0)

#+ gp?

gp_gr <- data.frame(sigma = c(0.01, 0.02, 0.04, 0.08, 0.16))
m2_gp <- train(y ~ ., data = p_impute, method = 'gaussprRadial', tuneGrid = gp_gr, trControl = trcntrl, na.action = na.omit)

plot(m2_gp)

plotCV(m2_gp)

m2_gp

#+ ranger, eval = TRUE

rf_gr <- expand.grid(mtry = c(2, 5, 10, 20, 30), splitrule = 'variance', min.node.size = c(5, 10, 20, 50))
m3_rf <- train(y ~ ., data = p_impute, method = 'ranger', tuneGrid = rf_gr, trControl = trcntrl, na.action = na.omit, importance = 'impurity', num.trees = 1000)

plot(m3_rf)

plotCV(m3_rf)

m3_rf


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

predictor_gp = Predictor$new(m2_gp, data = dplyr::select(p_impute, -y), y = p_impute$y)
predictor_rf = Predictor$new(m3_rf, data = dplyr::select(p_impute, -y), y = p_impute$y)


interact_gp = Interaction$new(predictor_gp)
interact_gp$results %>% arrange(desc(.interaction)) %>% head
plot(interact_gp)


interact_rf = Interaction$new(predictor_rf)
interact_rf$results %>% arrange(desc(.interaction)) %>% head
plot(interact_rf)


#+ which_inter

interact_ges_gp = Interaction$new(predictor_gp, feature = "X9.1_GestationLen_d")
interact_ges_gp$results %>% arrange(desc(.interaction)) %>% head
plot(interact_ges_gp)


interact_ges_rf = Interaction$new(predictor_rf, feature = "X9.1_GestationLen_d")
interact_ges_rf$results %>% arrange(desc(.interaction)) %>% head
plot(interact_ges_rf)



interact_lat_gp = Interaction$new(predictor_gp, feature = "X26.4_GR_MidRangeLat_dd")
interact_lat_gp$results %>% arrange(desc(.interaction)) %>% head
plot(interact_lat_gp)


interact_lat_rf = Interaction$new(predictor_rf, feature = "X26.4_GR_MidRangeLat_dd")
interact_lat_rf$results %>% arrange(desc(.interaction)) %>% head
plot(interact_lat_gp)

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





#' ## Using genus as categorical variable.

#+ ranger_gen


genus_dummy <- dummyVars(~ MSW05_Genus, data = p)
genus_dummy_data <- data.frame(predict(genus_dummy, newdata = p))

genus_data <- cbind(p_impute, genus_dummy_data)



rf_gr_gen <- expand.grid(mtry = c(5, 20, 30, 100, 300, 500), splitrule = 'variance', min.node.size = c(5, 10, 20, 50))
m3_rf_gen <- train(y ~ ., data = genus_data, method = 'ranger', tuneGrid = rf_gr_gen, 
                   trControl = trcntrl, na.action = na.omit, importance = 'impurity', 
                   save.memory = TRUE, num.trees = 1000)

plot(m3_rf_gen)

plotCV(m3_rf_gen)

m3_rf_gen


#' ## PGLS phylogenetic regression
#' Fit a standard pylogenetic regression with caper so we can check the later INLA models are reasonable.

#+ apriori_phyloreg



comp_data <- comparative.data(tree, cbind(p_impute, MSW05_Binomial = p$MSW05_Binomial)[1:100, ], 'MSW05_Binomial')

# Model will not run without log transforming.
apriori_formula_phylo <- y ~ X5.1_AdultBodyMass_g + X3.1_AgeatFirstBirth_d + X18.1_BasalMetRate_mLO2hr + 
                             X9.1_GestationLen_d + X16.1_LittersPerYear + X17.1_MaxLongevity_m


m0_pglm <- pgls(apriori_formula_phylo, data = comp_data)

m0_pglm

phylo_profile <- pgls.profile(m0_pglm)

plot(phylo_profile$logLik ~ phylo_profile$x, type = 'l')


#+ INLA_phyloreg_compare

# Code broadly copied from https://github.com/daijiang/phyr/blob/master/R/pglmm-utils.R#L54
nspp <- length(comp_data$phy$tip.label)
Vphy <- ape::vcv(comp_data$phy)
Vphy <- Vphy / max(Vphy)
Vphy <- Vphy / exp(determinant(Vphy)$modulus[1]/nspp)

order <- match(p$MSW05_Binomial[1:100], colnames(Vphy))
Vphy <- Vphy[order, order] # same order as species levels

apriori_form_inla <- y ~ X5.1_AdultBodyMass_g + X3.1_AgeatFirstBirth_d + X18.1_BasalMetRate_mLO2hr + 
                             X9.1_GestationLen_d + X16.1_LittersPerYear + X17.1_MaxLongevity_m + 
                             f(ii, model = 'generic0', Cmatrix = Vphy)


m0_inla_comp <- inla(apriori_form_inla, data = cbind(p_impute[1:100, ], ii = 1:100), control.predictor = list(compute = TRUE))
m0_inla_comp$summary.fixed
m0_inla_comp$summary.hyperpar




#+ INLA_phyloreg

comp_data_full <- comparative.data(tree, cbind(p_impute, MSW05_Binomial = p$MSW05_Binomial), 'MSW05_Binomial')

# Code broadly copied from https://github.com/daijiang/phyr/blob/master/R/pglmm-utils.R#L54
nspp <- length(comp_data_full$phy$tip.label)
Vphy <- ape::vcv(comp_data_full$phy)
Vphy <- Vphy / max(Vphy)
Vphy <- Vphy / exp(determinant(Vphy)$modulus[1]/nspp)

order <- match(p$MSW05_Binomial, colnames(Vphy))
Vphy <- Vphy[order, order] # same order as species levels

apriori_form_inla2 <- y ~ X5.1_AdultBodyMass_g + X3.1_AgeatFirstBirth_d + X18.1_BasalMetRate_mLO2hr + 
                             X9.1_GestationLen_d + X16.1_LittersPerYear + X17.1_MaxLongevity_m + 
                             f(ii, model = 'generic0', Cmatrix = Vphy)

# fit full model

m0_inla_comp_full <- inla(apriori_form_inla2, data = cbind(p_impute, ii = 1:nrow(p_impute)), 
                      control.predictor = list(compute = TRUE),
                      control.inla= list(strategy = "gaussian", int.strategy = "eb"))




m0_inla_comp_full$summary.fixed
m0_inla_comp_full$summary.hyperpar


#+ INLA_phyloreg_cv
# cross validation
m0_inla_comp2 <- list()

for(f in 1:5){

  cv_data <- cbind(p_impute, ii = 1:nrow(p_impute))

  # remove hold out data
  cv_data$y[folds[[f]]] <- NA

  m0_inla_comp2[[f]] <- inla(apriori_form_inla2, data = cv_data, 
                      control.predictor = list(compute = TRUE),
                      control.inla= list(strategy = "gaussian", int.strategy = "eb"))

}

#+ INLA_phyloreg_cv_analysis
# calc rmse and scatter plot

inla_apriori_cvdat <- data.frame(y = p_impute$y, predy = NA)

for(f in 1:5){
  inla_apriori_cvdat$predy[folds[[f]]] <- m0_inla_comp2[[f]]$summary.fitted.values$mean[folds[[f]]]
}


sqrt(mean((inla_apriori_cvdat$predy - inla_apriori_cvdat$y) ^ 2))

ggplot(inla_apriori_cvdat, aes(y, predy)) +
  geom_point() +
  #scale_y_continuous(limits = c(NA, 3)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth() 


#+ stacked_gen_setup, eval = FALSE

mlist <- list(m1_enet, m2_gp, m3_rf)

cv_preds <- lapply(mlist, best_tune_preds)
cv_preds <- lapply(seq_along(cv_preds), function (x) cv_preds[[x]]$predicted)
cv_preds <- do.call(cbind, cv_preds)

names(cv_preds) <- c('enet', 'gp', 'rf')

cv_preds$y <- p_impute$y # check order
cv_preds$ii <- seq_len(nrow(cv_preds))

stacked_form <- y ~  
  f(enet, model='clinear', range=c(0, Inf), initial=0) +
  f(gp, model='clinear', range=c(0, Inf), initial=0) +
  f(rf, model='clinear', range=c(0, Inf), initial=0) +
  f(ii, model = 'generic0', Cmatrix = Vphy)


#+ stacked_full, eval = FALSE

m0_inla_comp_full <- inla(stacked_form, data = cv_preds, 
                      control.predictor = list(compute = TRUE),
                      control.inla= list(strategy = "gaussian", int.strategy = "eb"))




m0_inla_comp_full$summary.fixed
m0_inla_comp_full$summary.hyperpar



#' # Point level understanding

#'  - understand individual points
#'    - lime


#+ lime, fig.width = 10, fig.height = 12

# Explain Tenrec ecaudatus. Largest litter size. And other large litters
l1 <- pred_head(m1_enet)
l2 <- pred_head(m2_gp)
l3 <- pred_head(m3_rf)


m1_lime <- lime(p_impute, m1_enet)
m1_explain <- explain(p_impute[l1, ], m1_lime, n_features = 20, feature_select = 'auto')

plot_features(m1_explain)


m2_lime <- lime(p_impute, m2_gp)
m2_explain <- explain(p_impute[l2, ], m2_lime, n_features = 20, feature_select = 'auto')

plot_features(m2_explain)


m3_lime <- lime(p_impute, m3_rf)
m3_explain <- explain(p_impute[l3, ], m3_lime, n_features = 20, feature_select = 'auto')

plot_features(m3_explain)











#+ session_info

sessionInfo()


#+ close_cluster, cache = FALSE

stopCluster(cl)


