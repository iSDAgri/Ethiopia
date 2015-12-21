#' Regularized regression ensembles of nutrient mass balances with spatial predictions
#' EthioSIS Mehlich-3 extractable P,K,S,Ca & Mg data from 255 Woredas
#' Spatial predictions derived from BART, DNN, GBM & RF models w spatial covariates
#' M. Walsh & J. Chen, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","glmnet","raster")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(glmnet)
require(raster)

# Data setup --------------------------------------------------------------
# Look at & run this first!
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Ethiopia/blob/master/ETM3_ensemble_setup.R"
# source_url(SourceURL)

# GLMNET stacking models --------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", number=10)

# V0 = ilr [P,K,S,Ca,Mg | Fv]
V0.ens <- train(V0 ~ V0_bart+V0_dnn+V0_gbm+V0_rf, data = etm3_cal,
                method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V0.ens)
v0.imp <- varImp(V0.ens)
plot(v0.imp, top=4)

# V3 = ilr [P,K | K,Ca,Mg]
V3.ens <- train(V3 ~ V3_bart+V3_dnn+V3_gbm+V3_rf, data = etm3_cal,
                method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V3.ens)
v3.imp <- varImp(V3.ens)
plot(v3.imp, top=4)

# V4 = ilr [K | Ca,Mg]
V4.ens <- train(V4 ~ V4_bart+V4_dnn+V4_gbm+V4_rf, data = etm3_cal,
                method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V4.ens)
v4.imp <- varImp(V4.ens)
plot(v4.imp, top=4)

# V5 = ilr [P | S]
V5.ens <- train(V5 ~ V5_bart+V5_dnn+V5_gbm+V5_rf, data = etm3_cal,
                method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V5.ens)
v5.imp <- varImp(V5.ens)
plot(v5.imp, top=4)

# V6 = ilr [Ca | Mg]
V6.ens <- train(V6 ~ V6_bart+V6_dnn+V6_gbm+V6_rf, data = etm3_cal,
                method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V6.ens)
v6.imp <- varImp(V6.ens)
plot(v6.imp, top=4)

# Stop parallelize
stopCluster(mc)

# Test set predictions ----------------------------------------------------
V0_ens <- predict(V0.ens, GRIDSv)
V3_ens <- predict(V3.ens, GRIDSv)
V4_ens <- predict(V4.ens, GRIDSv)
V5_ens <- predict(V5.ens, GRIDSv)
V6_ens <- predict(V6.ens, GRIDSv)
pred <- cbind.data.frame(V0_ens,V3_ens,V4_ens,V5_ens,V6_ens)
test <- etm3_val[c("PID","V0","V3","V4","V5","V6")]
ens_eval <- cbind(test, pred)

# Gridded ensemble predictions --------------------------------------------
V0_ens <- predict(grids, V0.ens)
V3_ens <- predict(grids, V3.ens)
V4_ens <- predict(grids, V4.ens)
V5_ens <- predict(grids, V5.ens)
V6_ens <- predict(grids, V6.ens)
ens_pred <- stack(V0_ens,V3_ens,V4_ens,V5_ens,V6_ens)
names(ens_pred) <- c("V0_ens","V3_ens","V4_ens","V5_ens","V6_ens")
plot(ens_pred)

# Export Gtif's -----------------------------------------------------------
# Create a "Results" folder in current working directory
dir.create("ETM3_results", showWarnings=F)

# Export Gtif's to "./ETM3_results"
writeRaster(ens_pred, filename="./ETM3_results/ens_preds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
