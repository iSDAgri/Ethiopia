#' PLS predictions of nutrient mass balance variables with spatial covariates
#' EthioSIS Mehlich-3 extractable P,K,S,Ca & Mg, from 255 Woredas
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","pls","raster")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(pls)
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

# Ensemble-PLS models -----------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", number = 10)

# V0 = ilr [P,K,S,Ca,Mg | Fv]
V0.pls <- train(GRIDSc, V0, 
                method = "pls", 
                preProc = c("center", "scale"),
                tuneGrid = expand.grid(ncomp = 1:10),
                trControl = tc)
print(V0.pls)
v0.imp <- varImp(V0.pls)
plot(v0.imp, top=28)

# V3 = ilr [P,K | K,Ca,Mg]
V3.pls <- train(GRIDSc, V3, 
                method = "pls", 
                preProc = c("center", "scale"),
                tuneGrid = expand.grid(ncomp = 1:10),
                trControl = tc)
print(V3.pls) 
v3.imp <- varImp(V3.pls)
plot(v3.imp, top=28)

# V4 = ilr [K | Ca,Mg]
V4.pls <- train(GRIDSc, V4, 
                method = "pls", 
                preProc = c("center", "scale"),
                tuneGrid = expand.grid(ncomp = 1:10),
                trControl = tc)
print(V4.pls)
v4.imp <- varImp(V4.pls)
plot(v4.imp, top=28)

# V5 = ilr [P | S]
V5.pls <- train(GRIDSc, V5, 
                method = "pls", 
                preProc = c("center", "scale"),
                tuneGrid = expand.grid(ncomp = 1:10),
                trControl = tc)
print(V5.pls)
v5.imp <- varImp(V5.pls)
plot(v5.imp, top=28)

# V6 = ilr [Ca | Mg]
V6.pls <- train(GRIDSc, V6, 
                method = "pls", 
                preProc = c("center", "scale"),
                tuneGrid = expand.grid(ncomp = 1:10),
                trControl = tc)
print(V6.pls)
v6.imp <- varImp(V6.pls)
plot(v6.imp, top=28)

stopCluster(mc)

# Test set predictions ----------------------------------------------------
V0_pls <- predict(V0.pls, GRIDSv)
V3_pls <- predict(V3.pls, GRIDSv)
V4_pls <- predict(V4.pls, GRIDSv)
V5_pls <- predict(V5.pls, GRIDSv)
V6_pls <- predict(V6.pls, GRIDSv)
pred <- cbind.data.frame(V0_pls,V3_pls,V4_pls,V5_pls,V6_pls)
test <- etm3_val[c("PID","V0","V3","V4","V5","V6")]
pls_eval <- cbind(test, pred)

# Gridded predictions -----------------------------------------------------
V0_pls <- predict(grids, V0.pls)
V3_pls <- predict(grids, V3.pls)
V4_pls <- predict(grids, V4.pls)
V5_pls <- predict(grids, V5.pls)
V6_pls <- predict(grids, V6.pls)
pls_pred <- stack(V0_pls,V3_pls,V4_pls,V5_pls,V6_pls)
names(pls_pred) <- c("V0_pls","V3_pls","V4_pls","V5_pls","V6_pls")
plot(pls_pred)

# Export Gtif's -----------------------------------------------------------
# Create a "Results" folder in current working directory
dir.create("ETM3_results", showWarnings=F)

# Export Gtif's to "./ETM3_results"
writeRaster(pls_pred, filename="./ETM3_results/pls_preds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
