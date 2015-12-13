#' Deep-Learning Neural Network predictions of nutrient mass balance variables with spatial covariates
#' EthioSIS Mehlich-3 extractable P,K,S,Ca & Mg, from 255 Woredas
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","deepnet","raster")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(deepnet)
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

# Random Forest models ----------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", number = 10)

# V0 = ilr [P,K,S,Ca,Mg | Fv]
V0.dnn <- train(GRIDSc, V0, 
                method = "dnn", 
                preProc = c("center", "scale"), 
                trControl = tc,
                tuneGrid = expand.grid(layer1 = 2:6,
                                       layer2 = 0:3,
                                       layer3 = 0:3,
                                       hidden_dropout = 0,
                                       visible_dropout = 0))
print(V0.dnn)
v0.imp <- varImp(V0.dnn, useModel = FALSE)
plot(v0.imp, top=28)

# V3 = ilr [P,K | K,Ca,Mg]
V3.dnn <- train(GRIDSc, V3, 
                method = "dnn", 
                preProc = c("center", "scale"), 
                trControl = tc,
                tuneGrid = expand.grid(layer1 = 2:6,
                                       layer2 = 0:3,
                                       layer3 = 0:3,
                                       hidden_dropout = 0,
                                       visible_dropout = 0))
print(V3.dnn) 
v3.imp <- varImp(V3.dnn, useModel = FALSE)
plot(v3.imp, top=28)

# V4 = ilr [K | Ca,Mg]
V4.dnn <- train(GRIDSc, V4, 
                method = "dnn", 
                preProc = c("center", "scale"), 
                trControl = tc,
                tuneGrid = expand.grid(layer1 = 2:6,
                                       layer2 = 0:3,
                                       layer3 = 0:3,
                                       hidden_dropout = 0,
                                       visible_dropout = 0))
print(V4.dnn)
v4.imp <- varImp(V4.dnn, useModel = FALSE)
plot(v4.imp, top=28)

# V5 = ilr [P | S]
V5.dnn <- train(GRIDSc, V5, 
                method = "dnn", 
                preProc = c("center", "scale"), 
                trControl = tc,
                tuneGrid = expand.grid(layer1 = 2:6,
                                       layer2 = 0:3,
                                       layer3 = 0:3,
                                       hidden_dropout = 0,
                                       visible_dropout = 0))
print(V5.dnn)
v5.imp <- varImp(V5.dnn, useModel = FALSE)
plot(v5.imp, top=28)

# V6 = ilr [Ca | Mg]
V6.dnn <- train(GRIDSc, V6, 
                method = "dnn", 
                preProc = c("center", "scale"), 
                trControl = tc,
                tuneGrid = expand.grid(layer1 = 2:6,
                                       layer2 = 0:3,
                                       layer3 = 0:3,
                                       hidden_dropout = 0,
                                       visible_dropout = 0))
print(V6.dnn)
v6.imp <- varImp(V6.dnn, useModel = FALSE)
plot(v6.imp, top=28)

stopCluster(mc)

# Test set predictions ----------------------------------------------------
V0_dnn <- predict(V0.dnn, GRIDSv)
V3_dnn <- predict(V3.dnn, GRIDSv)
V4_dnn <- predict(V4.dnn, GRIDSv)
V5_dnn <- predict(V5.dnn, GRIDSv)
V6_dnn <- predict(V6.dnn, GRIDSv)
pred <- cbind.data.frame(V0_dnn,V3_dnn,V4_dnn,V5_dnn,V6_dnn)
test <- etm3_val[c("PID","V0","V3","V4","V5","V6")]
dnn_eval <- cbind(test, pred)

# Gridded predictions -----------------------------------------------------
V0_dnn <- predict(grids, V0.dnn)
V3_dnn <- predict(grids, V3.dnn)
V4_dnn <- predict(grids, V4.dnn)
V5_dnn <- predict(grids, V5.dnn)
V6_dnn <- predict(grids, V6.dnn)
dnn_pred <- stack(V0_dnn,V3_dnn,V4_dnn,V5_dnn,V6_dnn)
names(dnn_pred) <- c("V0_dnn","V3_dnn","V4_dnn","V5_dnn","V6_dnn")
plot(dnn_pred)

# Export Gtif's -----------------------------------------------------------
# Create a "Results" folder in current working directory
dir.create("ETM3_results", showWarnings=F)

# Export Gtif's to "./ETM3_results"
writeRaster(dnn_pred, filename="./ETM3_results/dnn_preds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
