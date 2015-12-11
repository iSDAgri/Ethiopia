#' Generalized boosted regression predictions of nutrient mass balances with spatial covariates
#' EthioSIS Mehlich-3 extractable P,K,S,Ca & Mg data from 255 Woredas
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","gbm")), dependencies=TRUE)
require(devtools)
require(caret)
require(gbm)
require(raster)

# Data setup --------------------------------------------------------------
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Ethiopia/blob/master/ETM3_setup.R"
# source_url(SourceURL)

# Mehlich-3 nutrient mass balance variables
V0 <- etm3_cal$V0
V3 <- etm3_cal$V3
V4 <- etm3_cal$V4
V5 <- etm3_cal$V5
V6 <- etm3_cal$V6

# Gridded covariates
GRIDSc <- etm3_cal[c(18:45)] ## gridded covariates for model calibration from 204 Woredas
GRIDSv <- etm3_val[c(18:45)] ## same for 51 randomly selected validation Woredas

# GBM models --------------------------------------------------------------
set.seed(1385321)

# Cross-validation setup
tc <- trainControl(method = "cv", number=10)

# V0 = ilr [P,K,S,Ca,Mg | Fv]
V0.gbm <- train(GRIDSc, V0, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 3,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 100))
print(V0.gbm)
v0.imp <- varImp(V0.gbm)
plot(v0.imp, top=28)

# V3 = ilr [P,K | K,Ca,Mg]
V3.gbm <- train(GRIDSc, V3, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 3,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 100))
print(V3.gbm)
v3.imp <- varImp(V3.gbm)
plot(v3.imp, top=28)

# V4 = ilr [K | Ca,Mg]
V4.gbm <- train(GRIDSc, V4, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 3,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 100))
print(V4.gbm)
v4.imp <- varImp(V4.gbm)
plot(v4.imp, top=28)

# V5 = ilr [P | S]
V5.gbm <- train(GRIDSc, V5, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 3,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 100))
print(V5.gbm)
v5.imp <- varImp(V5.gbm)
plot(v5.imp, top=28)

# V6 = ilr [Ca | Mg]
V6.gbm <- train(GRIDSc, V6, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 3,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 100))
print(V6.gbm)
v6.imp <- varImp(V6.gbm)
plot(v6.imp, top=28)

# Test set predictions ----------------------------------------------------
V0_gbm <- predict(V0.gbm, GRIDSv)
V3_gbm <- predict(V3.gbm, GRIDSv)
V4_gbm <- predict(V4.gbm, GRIDSv)
V5_gbm <- predict(V5.gbm, GRIDSv)
V6_gbm <- predict(V6.gbm, GRIDSv)
pred <- cbind.data.frame(V0_gbm,V3_gbm,V4_gbm,V5_gbm,V6_gbm)
test <- etm3_val[c("PID","V0","V3","V4","V5","V6")]
eval <- cbind(test, pred)

# Gridded predictions -----------------------------------------------------
V0_gbm <- predict(grids, V0.gbm)
V3_gbm <- predict(grids, V3.gbm)
V4_gbm <- predict(grids, V4.gbm)
V5_gbm <- predict(grids, V5.gbm)
V6_gbm <- predict(grids, V6.gbm)
gbm_pred <- stack(V0_gbm,V3_gbm,V4_gbm,V5_gbm,V6_gbm)
names(gbm_pred) <- c("V0_gbm","V3_gbm","V4_gbm","V5_gbm","V6_gbm")
plot(gbm_pred)

# Export Gtif's -----------------------------------------------------------
writeRaster(gbm_pred, filename="gbm_pred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
