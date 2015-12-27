#' Ensemble predictions of potentially P,K and/or S deficient soils
#' EthioSIS Mehlich-3 extractable P,K & S measurements from 255 Woredas
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","glmnet","dismo","raster")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(glmnet)
require(dismo)
require(raster)

# Data setup --------------------------------------------------------------
# Look at & run this 1rst: https://github.com/mgwalsh/Ethiopia/blob/master/ETM3_ensemble_setup.R
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Ethiopia/blob/master/ETM3_ensemble_setup.R"
# source_url(SourceURL)

# Specify critical levels in ppm 
# P-test
Pcut <- 30
Pcal <- as.factor(ifelse(etm3_cal$P < Pcut, "Y", "N")) ## calibration data from 204 Woredas
Pval <- as.factor(ifelse(etm3_val$P < Pcut, "Y", "N")) ## validation data for 51 randomly selected Woredas

# K-test
Kcut <- 200
Kcal <- as.factor(ifelse(etm3_cal$K < Kcut, "Y", "N")) ## calibration data from 204 Woredas
Kval <- as.factor(ifelse(etm3_val$K < Kcut, "Y", "N")) ## validation data for 51 randomly selected Woredas

# S-test
Scut <- 20
Scal <- as.factor(ifelse(etm3_cal$S < Scut, "Y", "N")) ## calibration data from 204 Woredas
Sval <- as.factor(ifelse(etm3_val$S < Scut, "Y", "N")) ## validation data for 51 randomly selected Woredas

# Gridded covariates
GRIDSc <- etm3_cal[c(18:37)] ## gridded covariates for model calibration from 204 Woredas
GRIDSv <- etm3_val[c(18:37)] ## same for 51 randomly selected validation Woredas

# Model calibration -------------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 5, allowParallel = TRUE)

# P-test
P.glm <- train(GRIDSc, Pcal, method = "glmnet", 
               preProc = c("center", "scale"),
               trControl = tc)

# K-test
K.glm <- train(GRIDSc, Kcal, method = "glmnet", 
               preProc = c("center", "scale"),
               trControl = tc)

# S-test
S.glm <- train(GRIDSc, Scal, method = "glmnet", 
               preProc = c("center", "scale"),
               trControl = tc)

# Stop parallelize
stopCluster(mc)

# Model validation --------------------------------------------------------
# P-test
P_pred <- predict(P.glm, GRIDSv, type = 'prob')
P_test <- cbind(Pval, P_pred)
P_tpos <- subset(P_test, Pval=="Y", select=c(Y))
P_tneg <- subset(P_test, Pval=="N", select=c(Y))
P_eval <- evaluate(p=P_tpos[,1], a=P_tneg[,1]) ## calculate ROC's on test set <dismo>
plot(P_eval, 'ROC') ## plot ROC curve
P_thld <- threshold(P_eval, 'spec_sens') ## TPR+TNR threshold for classification
P_thld

# K-test
K_pred <- predict(K.glm, GRIDSv, type = 'prob')
K_test <- cbind(Kval, K_pred)
K_tpos <- subset(K_test, Kval=="Y", select=c(Y))
K_tneg <- subset(K_test, Kval=="N", select=c(Y))
K_eval <- evaluate(p=K_tpos[,1], a=K_tneg[,1]) ## calculate ROC's on test set <dismo>
plot(K_eval, 'ROC') ## plot ROC curve
K_thld <- threshold(K_eval, 'spec_sens') ## TPR+TNR threshold for classification
K_thld

# S-test
S_pred <- predict(S.glm, GRIDSv, type = 'prob')
S_test <- cbind(Sval, S_pred)
S_tpos <- subset(S_test, Sval=="Y", select=c(Y))
S_tneg <- subset(S_test, Sval=="N", select=c(Y))
S_eval <- evaluate(p=S_tpos[,1], a=S_tneg[,1]) ## calculate ROC's on test set <dismo>
plot(S_eval, 'ROC') ## plot ROC curve
S_thld <- threshold(S_eval, 'spec_sens') ## TPR+TNR threshold for classification
S_thld

# Spatial predictions -----------------------------------------------------
P_prob <- predict(grids, P.glm, type = 'prob')
P_mask <- 1-P_prob > P_thld
K_prob <- predict(grids, K.glm, type = 'prob')
K_mask <- 1-K_prob > K_thld
S_prob <- predict(grids, S.glm, type = 'prob')
S_mask <- 1-S_prob > S_thld
crit_pred <- stack(1-P_prob, 1-K_prob, 1-S_prob, P_mask, K_mask, S_mask)
names(crit_pred) <- c("P_prob","K_prob","S_prob","P_mask","K_mask","S_mask")
plot(crit_pred)

# Export Gtif's -----------------------------------------------------------
# Create a "Results" folder in current working directory
dir.create("ETM3_results", showWarnings=F)

# Export Gtif's to "./ETM3_results"
writeRaster(crit_ens, filename="./ETM3_results/crit_preds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
