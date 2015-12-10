#' Regularized regression predictions of nutrient mass balances with spatial covariates
#' EthioSIS Mehlich-3 extractable P,K,S,Ca & Mg data from 255 Woredas
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","randomForest")), dependencies=TRUE)
require(devtools)
require(caret)
require(glmnet)
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

# GLMNET models -----------------------------------------------------------
set.seed(1385321)

# Cross-validation setup
tc <- trainControl(method = "cv", number=10)

# V0 = ilr [P,K,S,Ca,Mg | Fv]
V0.rr <- train(GRIDSc, V0, method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V0.rr)
v0.imp <- varImp(V0.rr)
plot(v0.imp, top=28)

# V3 = ilr [P,K | K,Ca,Mg]
V3.rr <- train(GRIDSc, V3, method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V3.rr)
v3.imp <- varImp(V3.rr)
plot(v3.imp, top=28)

# V4 = ilr [K | Ca,Mg]
V4.rr <- train(GRIDSc, V4, method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V4.rr)
v4.imp <- varImp(V4.rr)
plot(v4.imp, top=28)

# V5 = ilr [P | S]
V5.rr <- train(GRIDSc, V5, method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V5.rr)
v5.imp <- varImp(V5.rr)
plot(v5.imp, top=28)

# V6 = ilr [Ca | Mg]
V6.rr <- train(GRIDSc, V6, method = "glmnet", preProc = c("center", "scale"), trControl = tc)
print(V6.rr)
v6.imp <- varImp(V6.rr)
plot(v6.imp, top=28)

# Test set predictions ----------------------------------------------------
V0_rr <- predict(V0.rr, GRIDSv)
V3_rr <- predict(V3.rr, GRIDSv)
V4_rr <- predict(V4.rr, GRIDSv)
V5_rr <- predict(V5.rr, GRIDSv)
V6_rr <- predict(V6.rr, GRIDSv)
pred <- cbind.data.frame(V0_rr,V3_rr,V4_rr,V5_rr,V6_rr)
test <- etm3_val[c("PID","V0","V3","V4","V5","V6")]
eval <- cbind(test, pred)

# Gridded predictions -----------------------------------------------------
V0_pred <- predict(grids, V0.rr)
V3_pred <- predict(grids, V3.rr)
V4_pred <- predict(grids, V4.rr)
V5_pred <- predict(grids, V5.rr)
V6_pred <- predict(grids, V6.rr)
